# Copyright (C) 2021 by University of Cambridge

# This software and algorithm was developed as part of the Cambridge Alliance
# for Medicines Safety (CAMS) initiative, funded by AstraZeneca and
# GlaxoSmithKline

# This program is made available under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the
# License, or at your option, any later version.
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


def fit_parameter(df: pd.DataFrame, parameter):
    """Fit the S curve of concentration-response

    Args:
        df: Dataframe from the dataset
        ax: Axes object from the matplotlib

    Returns:
        popt: parameters of the S curve
        perr: RMSE of the fitted curve

    Raise:
        RuntimeError: When the curve cannot be fitted
    """
    df = df.copy()
    freqs = df[parameter]
    concentrations = df.concentration
    hillcurve = HillCurve(concentrations, freqs)
    popt = hillcurve.popt
    perr = hillcurve.calc_perr()
    return popt, perr, hillcurve


def fsigmoid(x, a, x0, k, b):
    return k / (1.0 + np.exp(-a*(x-x0))) + b


def gain_loss(x, gw, ga, tp, lw, la, s, b):
    return tp * (1/(1+np.exp(gw*(ga-x)))+s) * (1/(1+np.exp(lw*(x-la)))) + b


def plain(x, b):
    return b

class TCPL:
    def __init__(self,
                 concentrations: np.ndarray,
                 responses: np.ndarray,
                 concentration_unit=-6,
                 boundary='auto',
                 logit=True):
        if (concentrations<=0).sum() > 0:
            raise ValueError("Concentration of 0 or less is invalid.")
        if logit:
            concentrations = concentrations.apply(
                np.log10) + concentration_unit
        if boundary == 'auto':
            c_min = concentrations.min()
            c_max = concentrations.max()
        else:
            c_max, c_min = boundary
        p_min = responses.min()
        p = responses - p_min
        p_diff = p.max()
        p = p / p_diff
        self.concentrations = concentrations
        self.p_response = p
        self.responses = responses
        self.p_min = p_min
        self.p_diff = p_diff
        self.c_max = c_max
        self.c_min = c_min
        self.rmsd = None
        self.bound = self.get_bound(c_max, c_min)
        if self.bound and self.fnc:
            popt, pcov = self.fit(self.fnc)
            self.popt, pcov = popt, pcov
    
    def get_bound(self, c_max, c_min):
        raise NotImplementedError()
        # return None

    def fit(self, fnc):
        popt, pcov = curve_fit(
            fnc, self.concentrations, self.p_response, method='dogbox',
            bounds=self.bound)
        return popt, pcov
    
    @property
    def EC50(self):
        if len(self.popt)>1:
            return self.popt[1]
        else:
            return None
    
    @property
    def E50(self):
        if self.EC50:
            return self.predict(self.EC50)
        else:
            return None

    @property
    def curve_diff(self):
        if len(self.popt>1):
            return np.abs(self.popt[2]) * self.p_diff
        else:
            return None
    
    @property
    def curve_min(self):
        if len(self.popt>1):
            return self.popt[-1] * self.p_diff + self.p_min
        else:
            return None
    
    @property
    def curve_max(self):
        if len(self.popt>1):
            return np.abs(self. popt[2]) * self.p_diff + self.p_min
        else:
            return None

    @property
    def RMSD(self):
        if self.rmsd is None:
            y_true = self.responses
            y_pred = np.vectorize(self.predict)(self.concentrations)
            self.rmsd = np.linalg.norm(y_true-y_pred)
        return self.rmsd

    def calc_perr(self):
        perr = np.sqrt(np.diag(self.pcov))
        return perr

    def predict(self, x):
        y = self.fnc(x, *self.popt)
        y = y * self.p_diff + self.p_min
        return y

class TCPLHill(TCPL):
    fnc = staticmethod(fsigmoid)
    def get_bound(self, c_max, c_min):
        return [0.3, c_min-1, -1, 0], [8, c_max+0.5, 1, 1]

    @property
    def hill(self):
        return self.popt[0] * np.sign(self.popt[3])


class TCPLGainLoss(TCPL):
    fnc = staticmethod(gain_loss)
    name='TCPL-GainLoss'
    def get_bound(self, c_max, c_min):
        return [0.3, c_min-1, 0.01, 0.3, c_min-1, 0, 0], [8, c_max+0.5, 1, 18, c_max+0.5, 1, 1]
        

class TCPLPlain(TCPL):
    fnc = staticmethod(plain)
    def get_bound(self, *args):
        return [-1], [1]

class HillCurve:
    # Learn from https://github.com/jbloomlab/neutcurve/blob/master/neutcurve/hillcurve.py
    # Deprecated
    def __init__(self,
                 concentrations: np.ndarray,
                 responses: np.ndarray,
                 concentration_unit=-6,
                 boundary='auto',
                 logc=True):
        if logc:
            concentrations = concentrations.apply(
                np.log10) + concentration_unit
        if boundary == 'auto':
            c_min = concentrations.min()
            c_max = concentrations.max()
        else:
            c_max, c_min = boundary
        p_min = responses.min()
        p = responses - p_min
        p_diff = p.max()
        p = p / p_diff
        popt, pcov = curve_fit(fsigmoid, concentrations, p,
                               method='dogbox', bounds=([1., c_min, -1, 0], [10, c_max, 1, 1]))
        self.popt = popt
        self.pcov = pcov
        self.p_min = p_min
        self.p_diff = p_diff
        self.c_max = c_max
        self.c_min = c_min
        self.concentrations = concentrations
        self.responses = responses

    @property
    def EC50(self):
        return self.popt[1]

    @property
    def curve_diff(self):
        return np.abs(self.popt[2]) * self.p_diff

    @property
    def hill(self):
        return self.popt[0] * np.sign(self.popt[3])

    def calc_perr(self):
        perr = np.sqrt(np.diag(self.pcov))
        return perr

    def predict(self, x):
        y = fsigmoid(x, *self.popt)
        y = y * self.p_diff + self.p_min
        return y
