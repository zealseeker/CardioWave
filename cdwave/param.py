# Copyright (C) 2021 by University of Cambridge

# This software and algorithm was developed as part of the Cambridge Alliance
# for Medicines Safety (CAMS) initiative, funded by AstraZeneca and
# GlaxoSmithKline

# This program is made available under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the
# License, or at your option, any later version.
import logging
from typing import List
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.decomposition import FactorAnalysis as FA
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import median_absolute_deviation
import pandas as pd
import numpy as np
from tqdm import tqdm
from cdwave import data
from cdwave import hillcurve
logger = logging.getLogger(__name__)


def parameter_projection(df: pd.DataFrame, parameters: list = None, method='tsne', n_components=2):
    if parameters is None:
        parameters = data.WaveformFull.parameter_names
    X = df[parameters].astype('float').values  # Transfer bool to float
    if method == 'pca':
        X = preprocessing.scale(X)
    models = {
        'tsne': TSNE,
        'pca': PCA,
        'fa': FA
    }
    model = models[method](n_components)
    newX = model.fit_transform(X)
    for i in range(n_components):
        df['X%d' % (i+1)] = newX[:, i]
    return df, model


def aggrate_parameters(df: pd.DataFrame,
                       parameters: list = None, method="median",
                       compound_column='uniname',
                       plates: dict = None):
    """Aggregate the parameters of the same compound under the same concentration
    by methods such as median or mean

    Args:
        df: Dataframe of the whole parameters
        parameters: The parameter list to process
        method: The method to aggrate the parameters
        plates: A dictionary with the key of compounds and values of list of
            plates to use.

    Returns:
        DataFrame: Dataframe with aggregated parameters.
    """
    res = []
    if parameters is None:
        parameters = data.WaveformFull.parameter_names
    for (concentration, compound), gdf in df.groupby(['concentration', compound_column]):
        if plates and compound in plates:
            gdf = gdf[gdf.plate.isin(plates[compound])]
        if len(gdf) == 0:
            continue
        s = gdf[parameters].agg(method)
        s['compound'] = compound
        s['concentration'] = concentration
        res.append(s)
    return pd.DataFrame(res)


def parameter_correlation(df: pd.DataFrame, parameters: list = None):
    """Calculate the correlation between the parameters"""
    if parameters is None:
        parameters = data.WaveformFull.parameter_names
    X = df[parameters].astype(float).values
    corr = np.corrcoef(X, rowvar=False)
    return corr


def remove_well_by_baseline(pdf: pd.DataFrame) -> list:
    """Remove wells by the quality of baseline (pre-measurement)

    When the quality of baseline is low under the following critera, the well
    should be removed.

    1. There is at least one multi-peak
    2. standard deviation of peak space is higher than 1
    3. maximum amplitude is lower than 100
    4. Some key time point (such as decay point) cannot be recognised

    Args:
        pdf: DataFrame of one plate

    Returns:
        list: A list of removing wells
    """
    df = pdf[pdf['state'] == 'prior']
    wells = df[(df['max_combo_peaks'] > 1) |
               (df['std_lambda'] > 1) |
               (df['maximum'] < 100) |
               (df['fail_analysis'] == 'FALSE')]['well'].tolist()
    return wells


def calc_rcv(x: np.ndarray) -> float:
    """Robust coefficient of variation

    Using the second approch in this paper
    https://arxiv.org/pdf/1907.01110.pdf

    .. math::

       MAD = med | x_i - m |

       RCV_M = 1.4826 * \\frac{MAD}{m}

    Args:
        x: A 1-d array of parameters

    Returns:
        float: robust coefficient of variation
    """
    m = x.median()
    xm = x - m
    mad = xm.median()
    return 1.4826 * mad / m


def remove_low_quality(df: pd.DataFrame):
    """Remove waveforms with low quality

    Wells of a plate will be removed if:

    1. Double peak in negative control
    2. High standard deviation of peak space in negative control
    3. Low quality in baseline. See :func:`~remove_well_by_baseline`
    
    The whole plate will be removed if RCV is higher than 0.
    See :func:`~calc_rcv`

    Args:
        df: A dataframe of all the samples with parameters

    Returns:
        tuple: A tuple containing a filtered dataframe and a dictionary of
            removed wells
    """
    pdfs = []
    removed_wells = {}
    for plate, pdf in df.groupby('plate'):
        removed_wells[plate] = {}
        removing_wells = remove_well_by_baseline(pdf)
        removed_wells[plate].update(
            {x: 'low baseline' for x in removing_wells})
        pdf = pdf[~pdf['well'].isin(removing_wells)]
        removing_wells = pdf[pdf['compound'].isin(['NegCtrl', 'DMSO']) &
                             ((pdf['max_combo_peaks'] > 1) | (pdf['std_lambda'] > 1))]['well'].tolist()
        removed_wells[plate].update({x: 'double peak' for x in removing_wells})
        pdf = pdf[~pdf['well'].isin(removing_wells)]
        rcv = calc_rcv(pdf[pdf['compound'].isin(['NegCtrl', 'DMSO'])]['freq'])
        if rcv >= 0.2:
            logger.debug(
                "Plate %s is removed due to the rcv of %s", plate, rcv)
            continue
        pdfs.append(pdf)
    return pd.concat(pdfs), removed_wells


def normalise_by_negctrl(df: pd.DataFrame,
                         standardiser: str = 'sdm',
                         standardisers: dict = None,
                         parameters: list = None,
                         control_compound: str = 'DMSO') -> pd.DataFrame:
    """Normalise the parameters by negative control of the plate
    Due to the fact that there will be one negative control in a plate, we
    use mean to aggragate the parameters.

    Args:
        df: DataFrame of parameters got from CardioWave
        standariser: Method to standardise the datak, including  
            sdm: Subtract and divide by median of negative control  
            sm: Subtract by median of negative control  
            smdmad: Subtract median and divide by median absolute deviation
        standardisers: A dictionary of which the keys are standarise methods
            and the values are parameters implementing the standardisers.
        control_compound: The name of control samples in the compound column

    Return:
        DataFrame: Normalised parameters
    """
    if parameters is None:
        parameters = data.WaveformFull.parameter_names
    # TODO remove specific names such as DMSO, NegCtrl
    control_dict = {'GSK': 'DMSO', 'AstraZeneca': 'NegCtrl'}
    agg_df_list = []
    for (plate, vendor), pdf in df.groupby(['plate', 'vendor']):
        if control_compound is None:
            control_name = control_dict[vendor]
        else:
            control_name = control_compound
        if standardisers is None:
            standardisers = {standardiser: parameters}
        control_samples = pdf[pdf['compound'] == control_name]
        assert pdf.well.duplicated().sum() == 0
        df_list = []
        for method in standardisers:
            samples = []
            ps = standardisers[method]
            if method == 'smdmad':
                mad = median_absolute_deviation(control_samples[ps])
            median = control_samples[ps].agg('median')
            for idx, item in pdf.iterrows():
                if method == 'sm':
                    sample = item[ps] - median
                elif method == 'sdm':
                    sample = (item[ps] - median).divide(median)
                elif method == 'smdmad':
                    sample = (item[ps] - median).divide(mad)
                else:
                    raise ValueError(f'Method "{method}" is invalid')
                samples.append(sample)
            df_list.append(pd.DataFrame(samples, index=pdf.index))
        tdf = pd.concat(df_list, axis=1)
        for col in ['compound', 'concentration', 'well', 'vendor', 'plate']:
            tdf[col] = pdf[col]
        agg_df_list.append(tdf)
    agg_df = pd.concat(agg_df_list)
    return agg_df


def normalise_by_baseline(df: pd.DataFrame,
                          subtract_params: list,
                          divide_params: list,
                          divide_only_params: list = None,
                          std_params: dict = None) -> pd.DataFrame:
    """Normalise the parameters by baseline of the well

    Args:
        subtract_params: Parameter list to be subtracted only.
        divide_params: Parameter list to be subtracted and divided.
        divide_only_params: Parameter list to be divided only.
        std_params: A dictionary mapping standard deviation parameters to its
            average parameters, such as `{'std_amplitude': 'avg_amplitude'}`.
            Parameters in this dictionary will be processed via following
            equation. :math:`std(A)= \\frac{std(A)}{\\overline{A}}`

    Return:
        DataFrame: Normalised parameters
    """
    if divide_only_params is None:
        divide_only_params = []
    if std_params is None:
        std_params = {}
    if 'index' not in df.columns:
        df = df.reset_index()
    if len(set(df.state) & set(['prior', 'treat'])) != 2:
        raise ValueError("states of treat and prior are required")
    agg_df_list = []
    for (plate, well), cdf in df.groupby(['plate', 'well']):
        baseline = cdf[cdf['state'] == 'prior']
        treat = cdf[cdf['state'] == 'treat']
        if len(baseline) != 1 or len(treat) != 1:
            continue
        baseline = baseline.iloc[0]
        treat = treat.iloc[0]
        sub_sample = treat[subtract_params] - baseline[subtract_params]
        div_sample = (treat[divide_params] - baseline[divide_params]).divide(
            baseline[divide_params])
        div_only_sample = treat[divide_only_params].divide(
            baseline[divide_only_params])
        std_sample = treat[std_params.keys()].divide(
            baseline[list(std_params.values())].values)
        sample_dict = {'compound': baseline['compound'], 'plate': plate,
                       'concentration': baseline['concentration'],
                       'well': well,
                       'vendor': baseline['vendor'],
                       'index': treat['index']}
        sample_dict.update(sub_sample.to_dict())
        sample_dict.update(div_sample.to_dict())
        sample_dict.update(div_only_sample.to_dict())
        sample_dict.update(std_sample.to_dict())
        agg_df_list.append(sample_dict)
    return pd.DataFrame(agg_df_list)


def select_concentration(df: pd.DataFrame, c: float, lim=(0.01, 10)):
    ids = []
    for (plate, compound), pdf in df.groupby(['plate', 'cmp']):
        tmp = pdf[pdf['concentration'] >= c].sort_values(by='concentration')
        if len(tmp) == 0:
            tmp = pdf.sort_values(by='concentration', ascending=False)
        if len(tmp) == 0:
            logger.warning('No sample in %s - %s', plate, compound)
        concentration = tmp.iloc[0]['concentration']
        if concentration / c > lim[1] or concentration / c < lim[0]:
            continue
        ids.append(tmp.index[0])
    return df.loc[ids]


def select_concentration_by_log(df: pd.DataFrame, c: float):
    ids = []
    for _, pdf in df.groupby(['plate', 'cmp']):
        tmp = pdf[['concentration']].copy()
        tmp['logc'] = np.log10(tmp['concentration'])
        tmp['diff_logc'] = np.abs(tmp['logc'] - np.log10(c))
        ids.append(tmp.sort_values(by='diff_logc').index[0])
    return df.loc[ids]


def calc_inflection(df: pd.DataFrame, parameters=None, min_n=4):
    if parameters is None:
        parameters = ['freq', 'rms', 'PW10_mean']
    items = []
    for (plate, cmp), cdf in tqdm(df.groupby(['plate', 'cmp'])):
        item = {'plate': plate, 'cmp': cmp}
        if len(cdf) <= min_n:
            continue
        for parameter in parameters:
            cs = cdf['concentration']
            fs = cdf[parameter]
            try:
                curve = hillcurve.HillCurve(cs, fs)
            except ValueError:
                print(plate, cmp)
                continue
            perr = curve.calc_perr()
            if perr[2] > np.abs(curve.popt[2]):
                # A magic number, to be improved
                item[f'{parameter}_EC50'] = -4.3
            else:
                item[f'{parameter}_EC50'] = curve.EC50
            # item[f'{parameter}_EC50_std'] = min(perr[1], 10)
            item[f'{parameter}_DIS'] = curve.curve_diff
            # item[f'{parameter}_DIS_std'] = min(perr[2], 10)
        items.append(item)
    return pd.DataFrame(items)


def select_plates(df: pd.DataFrame, t=0.2):
    """If the amplitude and freq of the lowest concentration is out of +- 0.2
    then remove the plate

    Args:
        df: input dataframe of the parameters
        t: threshold of the good quality

    Returns:
        dict: A dictionary with key of compounds and value of a list of available
        plates
    """
    res = {}
    for (plate, mol), cdf in df.groupby(['plate', 'uniname']):
        c_min = cdf['concentration'].min()
        tmp = cdf[(cdf['concentration'] == c_min) &
                  ((np.abs(cdf['avg_amplitude']) >= t) |
                   (np.abs(cdf['freq']) >= t))]
        if len(tmp) > 0:
            continue
        res.setdefault(mol, []).append(plate)
    return res


def linear_regression_with_logc(df: pd.DataFrame,
                                parameter: str,
                                remove_beatstop=True,
                                error='raise'):
    """Use linear regression to derive slope and intercept for a parameter

    Args:
        df: A dataframe that contains the samples needed.
        parameter: The name of the parameter that are included in the dataframe
        remove_beatstop: Whether to remove the samples of which 'beat_stop' is True
        error: The return k and b if all the samples are beat_stop or there is no sample.
            if `raise` is used, it will raise a ValueError from LinearRegression

    Returns:
        Tuple: slope (k) and intercept (b)
    """
    if parameter not in df or 'logc' not in df:
        raise KeyError("Do not have the parameter or logc")
    if remove_beatstop and 'beat_stop' not in df:
        raise KeyError("Do not have beat_stop column")
    if remove_beatstop:
        df = df[~df['beat_stop']]
    if error != 'raise' and len(df) == 0:
        return error, error
    reg = LinearRegression().fit(
        df['logc'].values.reshape(-1, 1), df[parameter])
    k = reg.coef_[0]
    b = reg.intercept_
    return k, b


def npoint_descriptor(df: pd.DataFrame, parameter: str, n: int):
    "This function only works for 8 concentrations"
    df = df.sort_values(by='concentration')
    res = []
    if n == 3:
        res.append(df.iloc[0][parameter])
        # res.append(df[parameter].median())
        res.append((df.iloc[3][parameter]+df.iloc[4][parameter])/2)
        res.append(df.iloc[len(df)-1][parameter])
    elif n == 4:
        curve = hillcurve.TCPLHill(df['concentration'], df[parameter])
        res.append(curve.curve_min)
        res.append(curve.curve_max)
        res.append(curve.E50)
        res.append(curve.EC50)
    elif n == 8:
        res = df[parameter].values
        if len(res) != 8:
            logger.warning(
                'This compound does not have 8 points - %d', len(res))
    else:
        raise ValueError("Only 3 or 8 points are supported.")
    return res


def calc_4_descriptors(df: pd.DataFrame, parameters: List[str], compounds: List[str]) -> pd.DataFrame:
    """Calculate four descriptors for each parameter, including minimum concentration,
    maximum concentration, median concentration and slope of the concentration-response"""
    suffixes = ['_min', '_median', '_max']
    res = {'compound': compounds}
    for p in parameters:
        for suffix in suffixes:
            res[p+suffix] = []
        res[p+'_slope'] = []
        for compound in compounds:
            tdf = df[df['compound'] == compound]
            k = npoint_descriptor(tdf, p, 3)
            for i, suffix in enumerate(suffixes):
                res[p+suffix].append(k[i])
            try:
                k, _ = linear_regression_with_logc(tdf, p, error=0)
            except (KeyError, ValueError) as e:
                logger.error("Cannot calculate the slope of %s - %s", compound, p)
                raise e
            res[p+'_slope'].append(k)
    kdf = pd.DataFrame(res)
    return kdf


def calc_grit(df: pd.DataFrame, parameter: str):
    control_dist = df.loc[df['compound']=='DMSO', [parameter]].values
    target_dist = df.loc[df['compound']!='DMSO', [parameter]].values
    assert len(df['compound'].unique()) == 2, 'Only two compounds can be included'
    control_std = control_dist.std()
    assert control_std !=0, 'DMSO has not variance.'
    scaler = StandardScaler()
    scaler.fit(control_dist)
    scores = scaler.transform(target_dist)
    return pd.Series(scores, index=target_dist.index)
