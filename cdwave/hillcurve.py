import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

def fit_parameter(df: pd.DataFrame, parameter, ax=None):
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
    if ax:
        hillcurve.plot(ax, ylabel=parameter)
    popt = hillcurve.popt
    perr = hillcurve.calc_perr()
    return popt, perr

def fsigmoid(x, a, x0, k, b):
    return k / (1.0 + np.exp(-a*(x-x0))) + b

class HillCurve:
    # Learn from https://github.com/jbloomlab/neutcurve/blob/master/neutcurve/hillcurve.py
    def __init__(self,
                 concentrations: np.ndarray,
                 responses: np.ndarray,
                 concentration_unit=-6,
                 boundary='auto',
                 logc=True):
        if logc:
            concentrations = concentrations.apply(np.log10) + concentration_unit
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

    def plot(self, ax=None, density=0.1, ylabel=None):
        if ax == None:
            ax = plt.subplots()
        cc = np.arange(self.c_min, self.c_max, density)
        vsigmoid = np.vectorize(lambda x: fsigmoid(
            x, *self.popt) * self.p_diff + self.p_min)
        ax.plot(cc, vsigmoid(cc))
        ax.plot(self.concentrations, self.responses, 'o', label="normal")
        ax.text(0.7, 0.2, 'EC50: {:.2f}'.format(self.EC50), transform=ax.transAxes)
        if ylabel:
            ax.set_ylabel(ylabel)
        return ax
        
    def calc_perr(self):
        perr = np.sqrt(np.diag(self.pcov))
        return perr
