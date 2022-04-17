# Copyright (C) 2021 by University of Cambridge

# This software and algorithm was developed as part of the Cambridge Alliance
# for Medicines Safety (CAMS) initiative, funded by AstraZeneca and
# GlaxoSmithKline

# This program is made available under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the
# License, or at your option, any later version.
import logging
import os
from multiprocessing import Pool
from collections.abc import Callable
from tqdm.auto import tqdm
import pandas as pd
from joblib import Parallel, delayed
from cdwave import fnc
from cdwave.data import WaveformFull, Dataset
from cdwave.fnc import Waveform, BloodPressure

logger = logging.getLogger(__name__)


def calc_parameter(waveform: WaveformFull) -> dict:
    """Calculate parameters of a waveform"""
    series = waveform.get_signal_series()
    wave = Waveform(series)
    if not wave.get_peaks():
        return None
    wave.analyse()
    try:
        r = wave.get_parameters()
    except Exception as e:
        logger.error('Cannot get parameters of %s', waveform)
        raise e
    return r


pbar_wrapper = {'x': None}


def default_process_fnc(status=-1, total=0):
    pbar = pbar_wrapper['x']
    if status == 0 and total != 0:
        pbar_wrapper['x'] = tqdm(total=total)
    elif status == 1:
        pbar.update()
    elif status == 2:
        pbar.close()


def calc_parameters_for_waveforms(dataset: Dataset,
                                  process_fnc: Callable = None,
                                  batch: int = 200,
                                  processes: int = None,
                                  custom_calculator: Callable = None):
    """Calculate parameter for waveforms

    Args:
        dataset: The waveform dataset.
        process_fnc: A processing function used to send out the progress,
            see `default_process_fnc`, which uses tqdm
        batch: Batch size for multi-processing.
        processes: Number of processors.
        custom_calculator: A custom calculator which can setup the custom thresholds.
            If `None`, `calc_parameter` will be used by default.
    """
    if processes is None:
        try:
            processes = int(os.environ['NUMEXPR_MAX_THREADS'])
        except KeyError:
            processes = 10
        except ValueError:
            processes = 10
    if custom_calculator is None:
        custom_calculator = calc_parameter
    n = dataset.size
    if process_fnc:
        process_fnc(status=0, total=n)
    reses = []
    batch_waveforms = []
    for i, waveform in enumerate(dataset.waveforms):
        batch_waveforms.append(waveform)
        if i % batch == batch - 1 or i == n - 1:
            with Pool(processes=processes) as pool:
                for _waveform in batch_waveforms:
                    reses.append(pool.apply_async(
                        custom_calculator, (_waveform,)))
                for j, res in enumerate(reses):
                    r = res.get()
                    if r is not None:
                        batch_waveforms[j].parameters = r
            reses = []
            batch_waveforms = []
        if process_fnc:
            process_fnc(status=1)
    if process_fnc:
        process_fnc(status=2)


def calc_parameters_with_threshold(dataset: Dataset,
                                   threshold: int = 100,
                                   processes: int = None,
                                   custom_calculator: Callable = None):
    """Calculate parameter for waveforms

    Args:
        dataset: The waveform dataset
        threshold: The threshold
        processes: Number of processors
        custom_calculator: A custom calculator which can setup the custom
            thresholds. If `None`, :func:`~calc_parameter_with_threshold` will
            be used by default.
    """
    if processes is None:
        try:
            processes = int(os.environ['NUMEXPR_MAX_THREADS'])
        except (KeyError, ValueError):
            processes = 4
    if custom_calculator is None:
        custom_calculator = calc_parameter_with_threshold
    logger.debug('Calculating parameters for %d waveforms with %d preceossors', len(dataset), processes)
    parameters = Parallel(processes)(
        delayed(custom_calculator)(_waveform, threshold) for _waveform in tqdm(dataset.waveforms))
    for waveform, p in zip(dataset.waveforms, parameters):
        if p is not None:
            waveform.parameters = p


def calc_parameter_with_threshold(waveform: WaveformFull, threshold, method='prominence') -> dict:
    """Calculate parameters of a waveform"""
    series = waveform.get_signal_series()
    wave = Waveform(series)
    if method == 'height':
        res = wave.get_peaks(height=threshold)
    elif method == 'prominence':
        res = wave.get_peaks(prominence=threshold)
    else:
        raise ValueError(f"Method {method} is not supported. Please use 'height' or 'prominence'")
    if not res:
        return None
    wave.analyse()
    try:
        r = wave.get_parameters()
    except Exception as e:
        logging.error('Cannot get parameters of %s', waveform)
        raise e
    return r


def derive_bp_parameters(bp: BloodPressure, start = 0, end = None, processes=4):
    if end is None:
        end = bp.max_time
    series = bp.series.loc[start: end]
    batch_generator = bp.get_batch_series(series)
    # Calculate total batches
    total = (bp.max_time - bp.window_size) // bp.batch_window
    sub_dfs = Parallel(processes)(
        delayed(derive_batch_bp_parameters)(batch_bp) for batch_bp in tqdm(batch_generator, total=total))
    df = pd.concat(sub_dfs)
    df['time(h)'] = df['start_time'] / 1000 / 3600
    df['time(m)'] = df['start_time'] / 1000 / 60
    return df


def derive_batch_bp_parameters(batch_bp: BloodPressure):
    batch_bp.run_filter()
    total, generator = batch_bp.get_windows_generator()
    if total == 0:
        return pd.DataFrame()
    rows = []
    for start_time, window in zip(batch_bp.get_start_times(), generator()):
        window = window[window['time_diff']>=50].set_index('group_id')
        if len(window) == 0:
            continue
        hrv = batch_bp.calc_hrv(window)
        tao = int(hrv['RR_mean'] // 3)
        tao = tao // batch_bp.sample_interval * batch_bp.sample_interval
        angle = batch_bp.calc_angle(start_time, tao)
        row = {
            'PP_min': window['PP'].min(),
            'PP_max': window['PP'].max(),
            'PP_mean': window['PP'].mean(),
            'PP_25': window['PP'].quantile(.25),
            'PP_75': window['PP'].quantile(.75),
            'start_time': window.time.values[0],
            'angle': angle
        }
        row.update(hrv)
        rows.append(row)
    return pd.DataFrame(rows)
