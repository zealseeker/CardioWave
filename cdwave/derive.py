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
from tqdm import tqdm
from joblib import Parallel, delayed
from cdwave.data import WaveformFull, Dataset
from cdwave.fnc import Waveform

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
        dataset: The waveform dataset
        process_fnc: A processing function used to send out the progress,
            see `default_process_fnc`, which uses tqdm
        batch: Batch size for multi-processing
        processes: Number of processors
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
        process_fnc: A processing function used to send out the progress,
            see `default_process_fnc`, which uses tqdm
        batch: Batch size for multi-processing
        processes: Number of processors
        custom_calculator: A custom calculator which can setup the custom thresholds.
            If `None`, `calc_parameter` will be used by default.
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


def calc_parameter_with_threshold(waveform: WaveformFull, threshold) -> dict:
    """Calculate parameters of a waveform"""
    series = waveform.get_signal_series()
    wave = Waveform(series)
    if not wave.get_peaks(height=threshold):
        return None
    wave.analyse()
    try:
        r = wave.get_parameters()
    except Exception as e:
        logging.error('Cannot get parameters of %s', waveform)
        raise e
    return r
