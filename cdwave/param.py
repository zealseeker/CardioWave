import logging
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.decomposition import FactorAnalysis as FA
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
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


def remove_well_by_baseline(pdf: pd.DataFrame):
    """Remove wells by baseline (prior)

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


def calc_rcv(x: np.ndarray):
    """Robust coefficient of variation

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
    The low quality includes:

    1. Double peak in negative control
    2. Low quality in baseline

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


def normalise_by_negctrl(df: pd.DataFrame, subtract_params: list, divide_params: list, control_compound=None):
    """Normalise the parameters by negative control of the plate
    Due to the fact that there will be one negative control in a plate, we
    use mean to aggragate the parameters.

    Args:
        subtract_params: Parameter list to be subtracted
        divide_params: Parameter list to be subtracted and then be divided
        control_compound: The name of control samples in the compound column

    Return:
        DataFrame: Normalised parameters
    """
    parameters = subtract_params + divide_params
    # TODO remove specific names such as DMSO, NegCtrl
    control_dict = {'GSK': 'DMSO', 'AstraZeneca': 'NegCtrl'}

    def aggregate_negative_control(pdf, name):
        samples = pdf[pdf['compound'] == name]
        return samples[parameters].agg('median')
    agg_df_list = []
    for (plate, vendor), pdf in df.groupby(['plate', 'vendor']):
        if control_compound is None:
            control_name = control_dict[vendor]
        else:
            control_name = control_compound
        control_param = aggregate_negative_control(pdf, control_name)
        for (compound, _), cdf in pdf.groupby(['compound', 'concentration']):
            # Sometimes  there is more than one compound-concentration in one plate
            for _, item in cdf.iterrows():
                sub_sample = item[subtract_params] - \
                    control_param[subtract_params]
                if (control_param[divide_params] == 0).any():
                    zero = control_param[divide_params]
                    zero = zero[zero == 0].index
                    logger.warning("Find zero in %s - %s", plate, zero)
                div_sample = (item[divide_params] - control_param[divide_params]).divide(
                    control_param[divide_params])
                tmp0 = item
                sample_dict = {'cmp': compound, 'plate': plate,
                               'concentration': tmp0['concentration'],
                               'vendor': vendor,
                               'well': tmp0['well'],
                               'index': tmp0['index']}
                sample_dict.update(sub_sample.to_dict())
                sample_dict.update(div_sample.to_dict())
                agg_df_list.append(sample_dict)
    agg_df = pd.DataFrame(agg_df_list)
    return agg_df


def normalise_by_baseline(df: pd.DataFrame, subtract_params: list, divide_params: list):
    """Normalise the parameters by baseline of the well

    Args:
        subtract_params: Parameter list to be subtracted only
        divide_params: Parameter list to be subtracted and divided

    Return:
        DataFrame: Normalised parameters
    """
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
        sample_dict = {'compound': baseline['compound'], 'plate': plate,
                       'concentration': baseline['concentration'],
                       'well': well,
                       'vendor': baseline['vendor'],
                       'index': treat['index']}
        sample_dict.update(sub_sample.to_dict())
        sample_dict.update(div_sample.to_dict())
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
    elif n == 8:
        res = df[parameter].values
        if len(res) != 8:
            logger.warning(
                'This compound does not have 8 points - %d', len(res))
    else:
        raise ValueError("Only 3 or 8 points are supported.")
    return res
