# Copyright (C) 2021 by University of Cambridge

# This software and algorithm was developed as part of the Cambridge Alliance
# for Medicines Safety (CAMS) initiative, funded by AstraZeneca and
# GlaxoSmithKline

# This program is made available under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the
# License, or at your option, any later version.
from itertools import product
import pandas as pd
from cdwave import param

model_parameters = ['up_length', 'down_length', 'full_down', 'rd_ratio', 'freq',
                    'avg_amplitude', 'avg_lambda', 'tail_proportion', 'avg_tail', 'std_tail',
                    'avg_shoulder', 'avg_shoulder_tail', 'std_shoulder',
                    'std_shoulder_tail', 'std_amplitude', 'std_lambda', 'PW10_mean', 'PW10_std',
                    'PW25_mean', 'PW25_std', 'PW50_mean', 'PW50_std', 'PW80_mean',
                    'PW80_std', 'PW90_mean', 'PW90_std']


def four_point_parameter_generator(parameters, suffixes):
    param_matrix = {'parameter': [], 'suffix': [], 'product': []}
    suffixes = ['_slope', '_median', '_max']
    for p, s in product(parameters, suffixes):
        param_matrix['parameter'].append(p)
        param_matrix['suffix'].append(s)
        param_matrix['product'].append(p + s)
    param_factor_df = pd.DataFrame(param_matrix).set_index('product')
    return param_factor_df


def prepare_four_point_model(agg_df: pd.DataFrame, endpoint: pd.Series, parameters, suffixes):
    res = {}
    cmps = agg_df['cmp'].unique()
    for p in parameters:
        for suffix in ['_min', '_median', '_max']:
            res[p+suffix] = []
        res[p+'_slope'] = []
        for cmp in cmps:
            df = agg_df[agg_df['compound'] == cmp]
            k = param.npoint_descriptor(df, p, 3)
            for i, suffix in enumerate(suffixes):
                res[p+suffix].append(k[i])
            k, _ = param.linear_regression_with_logc(df, p, error=0)
            res[p+'_slope'].append(k)
    kdf = pd.DataFrame(res)
    kdf['endpoint'] = kdf['compound'].map(endpoint)
    features = four_point_parameter_generator(
        parameters, suffixes).index.values
    X = kdf[features]
    y = kdf['ednpoint']
    return X, y, kdf
