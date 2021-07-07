# Copyright (C) 2021 by University of Cambridge

# This software and algorithm was developed as part of the Cambridge Alliance
# for Medicines Safety (CAMS) initiative, funded by AstraZeneca and
# GlaxoSmithKline

# This program is made available under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the
# License, or at your option, any later version.
import logging
import pandas as pd
import numpy as np
from scipy.stats import kstest
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.interpolate import UnivariateSpline
from scipy.signal import butter
from scipy.signal import filtfilt
from scipy.signal import hilbert
from scipy.signal import savgol_filter
from scipy.signal import welch
import statsmodels.nonparametric.api as smnp
logger = logging.getLogger(__name__)


class Waveform:
    """Class for waveform analysis

    Args:
        series (pd.Series): A series of which the name is the well name, values
            are amplitudes.
        index_penalty (float): Add penalty to signals to prioritise former time
            point during peak detection.

    Attributes:
        df (pd.DataFrame): The main dataframe of the waveform, including several important
            columns, peak: 1: main peak, 2-x: double peaks; status: 0: normal,
            1: raising point, 2: peak, 3: down starting; category: 0-9: intensities
            are categorised into 10 levels in terms of the `span`.
        num_peak (int): Number of peaks
        n (int): Number of points
        maximum (int): Maximum of intensity
        minimum (int): Minimum of intensity
    """
    max_shoulder_tail_ratio = 2.5

    def __init__(self, series: pd.Series, index_penalty: float = 0.000001):
        # opt_series is to prioritise early points
        opt_series = series - series.index * index_penalty
        self.series = opt_series
        self.series.name = series.name
        self.name = series.name
        self.num_peak = 0
        self.df = pd.DataFrame(series).reset_index()
        self.df['status'] = 0
        self.n = len(self.df)
        self.maximum = self.series.max()
        self.minimum = self.series.min()
        self.span = self.maximum - self.minimum
        self.variance = self.span / 10
        self.group_valley = {}
        self.fail_analysis = False
        self._group = None

    @property
    def group(self):
        """Return a iterator of group (i, gdf), excluding the first and the last"""
        def iterator():
            for i, gdf in self._group:
                if i in (self.num_peak, 0):
                    continue
                yield i, gdf
        if self._group is None:
            if self.num_peak == 0:
                raise ValueError("Cannot get group before peak detection")
            self._group = self.df.groupby('group')
        return iterator()

    def regroup(self):
        """After the peaks are changed, the groups need to be recalculated

        This function calculate the number of peaks and re-define the group
        number of each point.
        """
        df = self.df
        last_point = 0
        peaks = df[df['peak'] == 1]['peak']
        self.num_peak = len(peaks)
        for i, point in enumerate(peaks.index):
            self.df.loc[last_point:point, 'group'] = i
            last_point = point
        self.df.loc[last_point:, 'group'] = self.num_peak
        group = df.groupby('group')
        self._group = group
        for i, gdf in group:
            cond = gdf['peak'].notna() & (gdf['peak'] > 1)
            n_subpeak = cond.sum()
            index = gdf.loc[cond].index
            self.df.loc[index, 'peak'] = list(range(2, n_subpeak+2))

    def analyse_normal_waves(self, diff_n=1, l=3):
        """Get status of points for a normal wave

        Args:
            series: A series of which the name is the well name, values are amplitudes
            diff_n: n times of continuous points higher than this point to be
                regarded as a starting point
            l: length of exploration to find a point after higher than
                upper_half, indicating it can be a starting point
        """
        series = self.series
        # 1 Starting point of raising, 2: Peak point, 0: Normal (down)
        last_status = 0
        status = []
        diff = series.diff(diff_n)
        diff.iloc[0] = 0
        diff.iloc[1] = series.iloc[1] - series.iloc[0]
        upper_half = series.max() / 2
        for i, row in enumerate(
                self.df[['category', 'peak', 'group', self.name]].itertuples(index=False)):
            if i >= self.n - l:  # For the last group, don't calculte status
                status.append(0)
                continue
            if last_status == 0:  # Checking if it starts to raise
                for j in range(max(1, diff_n - i), l+1):
                    # One point between i and i+diff_n is smaller than i, reject the starting point
                    if diff.iloc[i+j] < 0:
                        status.append(0)
                        break
                    # Once a point between i and i+l is higher than upper_half,
                    # accept the starting point.
                    if series.iloc[i+j] > upper_half:
                        status.append(1)
                        last_status = 1
                        break
                else:
                    status.append(0)
            elif last_status < 2 and row[1] == 1:
                last_status = 2
                status.append(2)
            elif last_status == 2 and row[2] in self.group_valley and \
                    row[3] - self.group_valley[row[2]] < self.variance:
                status.append(3)
                last_status = 0
            else:
                status.append(0)
        self.df['status'] = status

    def analyse(self):
        """Analyse the status of each point

        The status indicates the thrend of the point, such as risng and
        declining. For waveforms with a frequency higher than 10, it will
        identify double peaks according their prominences and tail duration.
        This function also calculates the valley possitions.
        """
        self.get_valleys()
        # Check if there is less then 10 points between two peaks
        df = self.df.reset_index()
        diff = df[df['prominence'] != 0]['index'].diff().dropna()
        high_freq = diff.max() < 10
        if high_freq:
            peaks = self.df[self.df['prominence'] != 0].index
            valley_pos = self.df[self.df['valley'] == 1].index
            self.df.loc[valley_pos, 'status'] = 1
            self.df.loc[peaks, 'status'] = 2
        else:
            self.fix_double_peak_by_prominence()
            self.analyse_normal_waves()
            self.fix_double_peak_by_tail()
            self.analyse_normal_waves()

    def check_status_points(self):
        """A debugging function checking whether any period has more than 1 or
        has not starting point or downing point.
        """
        Success = True
        for i, gdf in self.df.groupby('group'):
            if i in (self.num_peak, 0):
                continue
            starting_point = gdf[gdf['status'] == 1]
            downing_point = gdf[gdf['status'] == 3]
            if len(starting_point) > 1:
                logger.debug("More than one starting point in group %d", i)
                Success = False
            elif len(starting_point) == 0:
                logger.debug("No starting point in group %d", i)
                Success = False
            if len(downing_point) > 1:
                logger.debug("More than one downing_point in group %d", i)
                Success = False
            elif len(downing_point) == 0:
                logger.debug("No downing point in group %d", i)
                Success = False
        return Success

    def get_peaks(self, height=None, prominence=None, min_prominence=20, span_ratio=0.1):
        """Identify the peaks and group of the whole waveform

        Args:
            height: Number or ndarray or sequence, optional Required height of peaks.
                Either a number, None, an array matching x or a 2-element sequence of the former.
                The first element is always interpreted as the minimal and the second,
                if supplied, as the maximal required height.
            prominence: Number or ndarray or sequence, optional Required prominence of peaks.
                Either a number, None, an array matching x or a 2-element sequence of the former.
                The first element is always interpreted as the minimal and the second, if supplied,
                as the maximal required prominence. If None, the minimal prominence will be the
                `max(min_prominence, span_ratio*self.span)`
            min_prominence: The absolute minimal prominence
            span_ratio: The ratio of span as the prominence threshold
        """
        series = self.df[self.name]
        self.df['category'] = (series - self.minimum) // self.variance
        if prominence is None:
            prominence = max(min_prominence, span_ratio*self.span)
        peaks, properties = find_peaks(
            series, height=height, prominence=prominence)
        prominences = properties['prominences']
        self.df['peak'] = 0
        self.df['prominence'] = 0
        for i in range(len(peaks)):
            peak = peaks[i]
            p = prominences[i]
            if i == 0:
                p = series[peak] - series[properties['right_bases'][i]]
            elif i == len(peaks)-1:
                p = series[peak] - series[properties['left_bases'][i]]
            self.df.loc[peak, 'peak'] = 1
            self.df.loc[peak, 'prominence'] = p
        self.num_peak = len(self.df[self.df['peak'] == 1])
        if self.num_peak == 0:
            return False
        self.regroup()
        return True

    def get_valleys(self):
        """Add valley infomation into the main dataframe"""
        self.df['valley'] = 0
        self.group_valley = {}
        for group_id, gdf in self.group:
            s = gdf[self.name]
            pos = s.idxmin()
            self.group_valley[group_id] = s[pos]
            self.df.loc[pos, 'valley'] = 1

    def fix_double_peak_by_prominence(self):
        """Identify double peaks by their prominence

        Definition of a double peak: A peak of which the prominence is below the
        threshold (see below about the threshold) and the signal value is close
        to the last real peak. Close means the difference between the signal
        value is smaller than variance (10% of maximum value).
        For waves of which the maximum is higher than 250, the threshold of the
        prominence is 0.7 * maixmum of all prominence. For others less than 250,
        the threshold is 0.5 * maximum of all prominence.
        Returns:
            Always True
        """
        if self.span < 250:
            prominences_t = 0.5
        else:
            prominences_t = 0.7
        # The order of current peak (1 is main peak, 2 is the first double peak)
        series = self.df[self.name]
        real_peaks = []
        peaks = self.df[self.df['peak'] == 1].index
        last_peak = None
        last_double_peak = None
        t = prominences_t * self.df.prominence.max()
        for peak in peaks:
            p = self.df.loc[peak, 'prominence']
            if p >= t:
                real_peaks.append(peak)
                self.df.loc[peak, 'peak'] = 1
                if last_double_peak and series[peak] - series[last_double_peak] < self.variance:
                    self.df.loc[last_double_peak, 'peak'] = 2
                last_peak = peak
            elif p < 0.3 * self.df.prominence.max():
                self.df.loc[peak, 'peak'] = 0  # Remove the peak
                continue
            elif last_peak and series[last_peak] - series[peak] < self.variance:
                self.df.loc[peak, 'peak'] = 2
                last_double_peak = None
            else:
                last_double_peak = peak
                self.df.loc[peak, 'peak'] = 0
        self.regroup()
        self.get_valleys()
        return True

    def calc_frequency_parameters(self):
        df = self.df
        lambdas = df[df['peak'] == 1].diff().iloc[1:]['time']
        avg_lambda = lambdas.mean() if len(lambdas) > 0 else 0
        std_lambda = lambdas.std() if len(lambdas) > 1 else 0
        if len(lambdas) > 5:
            quartiles = np.percentile(lambdas, [25, 50, 75])
            iqr = quartiles[-1] - quartiles[0]
            top = 1.5 * iqr + quartiles[-1]
            bottom = quartiles[0] - 1.5 * iqr
            inner_lambdas = lambdas[(lambdas <= top) & (lambdas >= bottom)]
            avg_inner_lambda = inner_lambdas.mean()
            std_inner_lambda = inner_lambdas.std()
        else:
            avg_inner_lambda = avg_lambda
            std_inner_lambda = std_lambda
        freq = (1 / avg_lambda * 100) if avg_lambda != 0 else 0
        return {
            'avg_lambda': avg_lambda,
            'std_lambda': std_lambda,
            'avg_inner_lambda': avg_inner_lambda,
            'std_inner_lambda': std_inner_lambda,
            'freq': freq
        }

    def calc_status_parameter(self):
        ups = []
        downs = []
        full_downs = []
        rd_ratios = []  # Raise decay ratio
        tails = []
        t_proportion = []  # Tail proportion
        for _, gdf in self.group:
            peak_time = gdf.iloc[0]['time']
            try:
                starting_point = gdf[gdf['status'] == 1].iloc[0]['time']
                downing_point = gdf[gdf['status'] == 3].iloc[0]['time']
            except IndexError:
                self.fail_analysis = True
                continue
            # gdf.index[-1] + 1 is next peak index
            next_peak_time = self.df.iloc[gdf.index[-1] + 1]['time']
            up = next_peak_time - starting_point
            down = downing_point - peak_time
            full_down = starting_point - peak_time
            tail = starting_point - downing_point
            ups.append(up)
            downs.append(down)
            full_downs.append(full_down)
            tails.append(tail)
            rd_ratios.append(up / full_down)
            t_proportion.append(tail / (up + full_down))
        if tails:
            avg_tail = np.mean(tails)
            std_tail = np.std(tails)
            t_proportion = np.mean(t_proportion)
        else:
            avg_tail = std_tail = t_proportion = 0
        return {
            'up_length': np.mean(ups) if ups else 0,
            'down_length': np.mean(downs) if downs else 0,
            'full_down': np.mean(full_downs) if full_downs else 0,
            'rd_ratio': np.mean(rd_ratios) if rd_ratios else 0,
            'avg_tail': avg_tail, 'std_tail': std_tail, 'tail_proportion': t_proportion
        }

    def calc_shoulder_parameters(self):
        """Calculate shoulder related parameters

        Shoulder parameters include mean and standard deviation of
        shoulder position (in amplitude), and median and standard deviation of
        shoulder_tail ratio. Median is used because sometime the ratio is
        unstable. For example, when tail is missed in one period, the ratio
        will become 100, which is cover other normal values.
        """
        shoulder_amplitudes = []
        shoulder_tail_ratio = []
        shoulder_position = []
        for _, gdf in self.group:
            amplitude = gdf.iloc[0][self.name]
            wdist = self.get_width_distribution(gdf)
            if not wdist:
                continue
            grid, y = wdist
            wd_peaks, wd_properties = find_peaks(
                y, prominence=0)  # Width Distribution
            peak_amplitudes = grid[wd_peaks]
            prominences = wd_properties['prominences']
            wd_tail = 0
            wd_shoulder_prominence = 0
            wd_shoulder = 0
            for wd_peak, prominence in zip(peak_amplitudes, prominences):
                if wd_peak > 0.80 * amplitude:
                    wd_tail = max(wd_tail, prominence)
                elif wd_peak > 0.15 * amplitude:
                    if prominence > wd_shoulder_prominence:
                        wd_shoulder_prominence = prominence
                        wd_shoulder = wd_peak
            shoulder_amplitudes.append(wd_shoulder)
            shoulder_position.append(wd_shoulder/amplitude)
            if wd_tail == 0 and wd_shoulder != 0:
                shoulder_tail_ratio.append(self.max_shoulder_tail_ratio)
            elif wd_tail != 0:
                shoulder_tail_ratio.append(
                    min(self.max_shoulder_tail_ratio, wd_shoulder_prominence/wd_tail))
            else:
                shoulder_tail_ratio.append(0)
        return {
            'avg_shoulder': np.mean(shoulder_position) if shoulder_position else 0,
            'avg_shoulder_tail': np.median(shoulder_tail_ratio) if shoulder_tail_ratio else 0,
            'std_shoulder': np.std(shoulder_position) if shoulder_position else 0,
            'std_shoulder_tail': np.std(shoulder_tail_ratio) if shoulder_tail_ratio else 0,
            'avg_shoulder_amp': np.mean(shoulder_amplitudes) if shoulder_amplitudes else 0,
            'std_shoulder_amp': np.std(shoulder_amplitudes) if shoulder_amplitudes else 0,
        }

    def calc_peak_width(self):
        df = self.df
        peaks = df[df['peak'] == 1].index
        rel_height_dict = {
            'PW10': 0.1,
            'PW25': 0.25,
            'PW50': 0.5,
            'PW80': 0.8,
            'PW90': 0.9
        }
        parameters = {name + '_mean': 0 for name in rel_height_dict}
        parameters.update({name + '_std': 0 for name in rel_height_dict})
        for name, rel_height in rel_height_dict.items():
            widths = peak_widths(df[self.name], peaks, rel_height)[0]
            parameters[name+'_mean'] = widths.mean()
            parameters[name+'_std'] = widths.std()
        return parameters

    def calc_amplitudes(self):
        df = self.df
        intensities = df.loc[df['peak'] == 1, self.name].values
        valleys = list(self.group_valley.values())
        amplitudes = []
        rms = np.sqrt(np.mean(np.square(df[self.name].values)))
        for i, v in enumerate(valleys):
            amplitudes.append(intensities[i] - v)
            amplitudes.append(intensities[i+1] - v)
        amplitudes = np.array(amplitudes)
        parameters = {
            'avg_amplitude': amplitudes.mean() if valleys else 0,
            'std_amplitude': amplitudes.std() if valleys else 0,
            'maximum': amplitudes.max() if valleys else 0,
            'avg_intensity': intensities.mean(),
            'std_intensity': intensities.std(),
            'max_intensity': self.maximum,
            'min_intensity': self.minimum,
            'avg_valley': np.mean(valleys) if valleys else 0,
            'rms': rms
        }
        return parameters

    def calc_fft_freq_ratio(self):
        r"""Calculate the energy ratio between double frequency (minor) and major

        .. math::

            FFT Ratio = \frac{\sum_{x=f_{mi}-0.05}^{f_{mi}+0.05}p(x)}
                        {\sum_{x=f_{ma}-0.05}^{f_{ma}+0.05}p(x)}
        """
        signals = self.resample(100)
        frq, psd = wave_transform(signals)
        # Remove peaks of which frq < 0.1 Sometimes butter doesn't work well
        sorted_frq_psd = sorted(
            filter(lambda x: x[0] > 0.05, zip(frq, psd)),
            key=lambda x: x[1], reverse=True)
        # In some cases, double freq may have higher density than main freq
        # We should select the smaller one (ideally we should select the one
        # closer to real peak freqency)
        if sorted_frq_psd[0][1] / sorted_frq_psd[1][1] < 2 and \
                sorted_frq_psd[0][0] > sorted_frq_psd[1][0]:
            # Swap the first two
            sorted_frq_psd[0], sorted_frq_psd[1] = sorted_frq_psd[1], sorted_frq_psd[0]
        main_frq, main_psd = sorted_frq_psd[0]
        if main_frq == 0 or main_psd == 0:
            return 0
        for _frq, _psd in sorted_frq_psd[1:]:
            if 1.9 < _frq/main_frq < 2.1:
                minor_frq = _frq
                break
        else:
            return 0
        S_main_fft = psd[(main_frq-0.05 < frq) & (frq < main_frq+0.05)]
        P_main = np.sum(S_main_fft) * (frq[1] - frq[0])
        S_minor_fft = psd[(minor_frq-0.05 < frq) & (frq < minor_frq+0.05)]
        P_minor = np.sum(S_minor_fft) * (frq[1] - frq[0])
        if P_minor == 0:
            return 0
        return P_minor / P_main

    def get_parameters(self):
        """Get all parameters from the waveform"""
        df = self.df
        is_double_peak = (df['peak'] > 1).any()
        max_combo_peaks = df['peak'].max()
        if self.num_peak < 4:
            ks_p = 1
        elif self.num_peak < 6:
            ks_p = self.peak_uniform_test(interpolation=True)
        else:
            ks_p = self.peak_uniform_test()
        parameters = {'double_peak': is_double_peak,
                      'uniform': ks_p > 0.99, 'noise': self._noise_test(),
                      'n_peak': self.num_peak, 'max_combo_peaks': max_combo_peaks,
                      'fail_analysis': self.fail_analysis,
                      'fft_ratio': self.calc_fft_freq_ratio()
                      }
        parameters.update(self.calc_frequency_parameters())
        parameters.update(self.calc_status_parameter())
        parameters.update(self.calc_shoulder_parameters())
        parameters.update(self.calc_peak_width())
        parameters.update(self.calc_amplitudes())
        return parameters

    def draw_status(self, figure):
        """Plot the points into a figure colored by their status
        black indicates normal points; orange is starting point, the first point
        of raising wave; red is the peak and green is the first point of the
        tail.
        Args:
            figure: A Matplotlib figure object
        """
        figure.clear()
        colors = {0: 'black', 1: 'orange', 2: 'red', 3: 'green'}
        ax1 = figure.add_subplot(111)
        for color in colors:
            tdf = self.df[self.df['status'] == color]
            ax1.plot(tdf['time'], tdf[self.series.name],
                     'o', color=colors[color])
        tdf = self.df[(self.df['peak'] > 1)]
        ax1.plot(tdf['time'], tdf[self.series.name], 'mo')

    def draw_series(self, figure, series=None, style='-'):
        """Plot points of the waveform
        Args:
            figure: A Matplotlib figure object
            series: A pandas Series of the waveform
            style: A Maplotlib line style, `-` by default indicating a line plot
        """
        if series is None:
            series = self.series
        elif type(series) == str:
            series = self.df[series]
        figure.clear()
        ax1 = figure.add_subplot(111)
        ax1.plot(series.index, series, style)
        return ax1

    def export(self, filename):
        """Export the waveform to a csv file"""
        self.df.to_csv(filename)

    def peak_uniform_test(self, interpolation=False):
        """Test if the peak points are in uniform distribution

        Args:
            interpolation: True if interpolation is applied to the waveform.
                Mostly used when ther are only 3 peaks in the waveform so we
                need to interpolate the points between the three points.
                Otherwise KS test may underestimate the possibility

        Returns:
            float: The probability of of waveform being in uniform distribution
        """
        df = self.df
        x = df[df['peak'] == 1].index
        # First to check if the minimum and maximum is close to the boundary
        diff_max = np.diff(x).max()
        if x[0] > diff_max or x[-1] + diff_max < df.index.max():
            return 0
        x = x - min(x)
        x = x / max(x)
        if interpolation:
            new_x = [x[0]]
            for i in range(len(x)-1):
                new_x.append(x[i]+(x[i+1]-x[i])/2)
                new_x.append(x[i+1])
            x = new_x
        _, p = kstest(x, 'uniform')
        return p

    def _noise_test(self):
        '''(Deprecated) Test if this waveform is noise'''
        sum_error = 0
        for _, gdf in self.group:
            p1 = gdf.iloc[0][self.name]
            p2 = self.df.loc[gdf.index[-1] + 1, self.name]
            p_max = max(p1, p2)
            for x in gdf[self.name]:
                if x > p_max:
                    sum_error += 1
                if sum_error > 2:
                    return True
        return False

    def get_width_distribution(self, gdf: pd.DataFrame, normalise=True):
        """Calculate the distirbution of width between position of the peak point
        and the points.

        Args:
            gdf: subdataframe of a group

        Returns:
            tuple: A tuple of grid and y. `grid` is grid sampled in all widthes.
            `y` is the probability of the width
        """
        peak_amplitude = gdf.iloc[0][self.name]
        next_start = gdf[gdf['status'] == 1]
        if len(gdf) < 10 or len(next_start) == 0:
            return False
        next_start = next_start.index[0]
        prs = []
        for x in gdf.loc[:next_start, self.name]:
            pr = peak_amplitude - x
            prs.append(pr)
        kde = smnp.KDEUnivariate(prs)
        kde.fit('biw', bw=0.2*peak_amplitude, fft=False, gridsize=100,
                cut=3, clip=(-np.inf, np.inf))  # TODO hyper-parameter bw
        grid, y = kde.support, kde.density
        if normalise:
            y = y / y.max()  # Normalise the amplitude of density probability function
        return grid, y

    def fix_double_peak_by_tail(self, min_amplitude=100, std_threshold=0.5):
        """Identify double peaks by comparing the average signal values of tails

        In some situations, "real" peaks have a long tail but the subpeaks
        don't, so we can recognise the subpeaks by comparing their tail with
        maximum. But for some waveforms, the tails are too short to be used as a
        symbol to identify double peaks.

        Args:
            min_amplitude: Minimum of amplitude the wave should have as a
                prerequisite to find double peak.
            std_threshold: Number of times that the standard deviation of the
                tails is higher than mean of the tails which will indicate
                there is double peak.
        Returns:
            bool: If starting point or downing point were not found, return
                False and self.fail_analysis will be True. Otherwise always
                return True
        """
        if self.maximum <= min_amplitude:
            return False
        tails = []
        group = self.df.groupby('group')
        for i, gdf in group:
            if i in (0, self.num_peak):
                continue
            starting_point = gdf[gdf['status'] == 1]
            downing_point = gdf[gdf['status'] == 3]
            if len(downing_point) == 0:
                self.fail_analysis = True
                return False
            if len(starting_point) == 0:
                # This means starting_point is the same to downing point
                tails.append(0)
            else:
                tails.append(
                    starting_point.iloc[0]['time'] - downing_point.iloc[0]['time'])
        tails = np.array(tails)
        if len(tails) == 0:
            return False
        std_tail = tails.std()
        avg_tail = tails.mean()
        max_tail = tails.max()
        need_regroup = False
        if std_tail > std_threshold * avg_tail:
            for i, tail in enumerate(tails):
                if tail < max_tail - std_tail:
                    self.df.loc[group.get_group(i+1).index[0], 'peak'] = 2
                    need_regroup = True
        if need_regroup:
            self.regroup()
            self.get_valleys()  # TODO Maybe it can be merged to regroup
        return True

    def resample(self, sample_rate=100) -> np.ndarray:
        """Resample points in the waveform

        Args:
            sample_rate: How many points per second, default 100

        Returns:
            np.ndarray: An numpy array with `sample_rate` points per second
        """
        x = self.df['time']
        y = self.df[self.name]
        # new_y = savgol_filter(y, 5, 2)
        x_new = np.linspace(x.min(), x.max(), int(x.max()*sample_rate+1))
        interp = UnivariateSpline(x, y, k=3)
        return interp(x_new)

    def standardise_by_filter(self):
        signals = self.series.values
        new_signals = signal_filter(signals)
        series = pd.Series(
            new_signals, index=self.series.index, name=self.series.name)
        self.__init__(series)


def wave_transform(signals: np.ndarray, sample_rate=100, method='fft'):
    """Transform the wave using methods such as Fourier transformation

    Args:
        signals (np.ndarray): resampled signals from waveform
        sample_rate (int): How many points per second, default 100
        method (str): Transformation method, default `fft`, fast fourier transform

    Returns:
        tuple: A tuple containing:
            - **frq** (*np.ndarray*): Frequency points (Hz)
            - **psd** (*np.ndarray*): Power Spectral Density from FFT
    """
    freq_range = np.array([0.01, 1])  # frequency from 0.03 to 1 Hz
    # Smooth by
    b, a = butter(2, freq_range/0.5/sample_rate, 'bandpass')
    signals = filtfilt(b, a, signals)
    if method == 'fft':
        datalen = len(signals)
        frq = np.fft.fftfreq(datalen, d=((1/sample_rate)))
        frq = frq[range(int(datalen/2))]
        Y = np.fft.fft(signals)/datalen
        Y = Y[range(int(datalen/2))]
        psd = np.power(np.abs(Y), 2)
    elif method == 'hilbert':
        frq = (np.arange(len(signals)) / sample_rate)[1:]
        analytic_signal = hilbert(signals)
        # amplitude_envelope = np.abs(analytic_signal)
        instantaneous_phase = np.unwrap(np.angle(analytic_signal))
        instantaneous_frequency = (np.diff(instantaneous_phase) /
                                   (2.0*np.pi) * sample_rate)
        psd = instantaneous_frequency
    elif method == 'welch':
        frq, psd = welch(signals, fs=100, nperseg=len(signals)-1)
    else:
        raise ValueError('The method is not available')
    return frq, psd


def signal_filter(signals: np.ndarray, window_length=5) -> np.ndarray:
    """Apply a Savitzky-Golay filter to an array."""
    return savgol_filter(signals, window_length, 2)
