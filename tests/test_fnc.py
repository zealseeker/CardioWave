import unittest
import warnings
import numpy as np
import pandas as pd
from cdwave import fnc
from tests import Dataset


class TestFnc(unittest.TestCase):
    data = Dataset()

    def test_find_peak(self):
        df = self.data['D1']
        series = df.set_index('time')['B14']
        wave = fnc.Waveform(series)
        self.assertTrue(wave.get_peaks(), "Cannot get peaks")
        wave.analyse()
        self.assertEqual(wave.num_peak, 17)

    def test_double_peak(self):
        df = self.data['double_peak']
        series = df.set_index('time')['F6']
        wave = fnc.Waveform(series)
        self.assertTrue(wave.get_peaks(), "Cannot get peaks")
        wave.fix_double_peak_by_tail()
        wave.analyse()
        r = wave.get_parameters()
        self.assertEqual(r['n_peak'], 4)
        self.assertEqual(r['max_combo_peaks'], 3)

    def test_width_distribution(self):
        df = self.data['D1']
        series = df.set_index('time')['B14']
        wave = fnc.Waveform(series)
        self.assertTrue(wave.get_peaks(), "Cannot get peaks")
        wave.analyse()
        r = wave.get_parameters()
        self.assertGreater(r['avg_shoulder_tail'], 2)

    def test_high_frequency(self):
        wave = self.data.get_waveform('high_freq')
        self.assertTrue(wave.get_peaks(), "Cannot get peaks")
        wave.analyse()
        r = wave.get_parameters()
        self.assertEqual(r['max_combo_peaks'], 1)
        self.assertEqual(r['up_length'], 0)
        self.assertFalse(r['fail_analysis'])
        self.assertTrue(r['uniform'])

    def test_one_peak(self):
        wave = self.data.get_waveform('one_peak')
        self.assertTrue(wave.get_peaks(), "Cannot get peaks")
        wave.analyse()
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            r = wave.get_parameters()
            for w in ws:
                print('\n==========\n{}-{}\n{}'.format(w.filename, w.lineno, w.message))
            self.assertEqual(len(ws), 0, "Warnning has been raised {}".format(
                [x.message for x in ws]))
        self.assertEqual(r['n_peak'], 1)
        for p in r:
            self.assertFalse(np.isnan(r[p]))

    def test_uniform(self):
        wave = self.data.get_waveform('uniform')
        wave.get_peaks()
        wave.analyse()
        self.assertEqual(wave.peak_uniform_test(), 0)

    def test_unstable_waveform(self):
        """
        Unstable means that the valleys are fluctuating during the waveform,
        so the prominences will be "inaccurate" leading to the fault double peak
        This is not easy to fix but we can try to use "QC"
        This test may fail if slightly change the hyper-parameter (threshold) or
        algorithm.
        """
        wave = self.data.get_waveform('unstable')
        wave.standardise_by_filter()
        wave.get_peaks()
        wave.analyse()
        r = wave.get_parameters()
        self.assertEqual(r['max_combo_peaks'], 1)

    def test_fft_freq_ratio(self):
        wave = self.data.get_waveform('fft')
        wave.standardise_by_filter()
        wave.get_peaks()
        wave.analyse()
        r = wave.calc_fft_freq_ratio()
        self.assertAlmostEqual(r, 1.1, 1)


class TestBloodPressure(unittest.TestCase):
    times = np.arange(0, 100000, 2)
    signal = np.sin(times/10)*100 + np.random.random(len(times))*10
    data = pd.Series(signal, index=times)

    def test_batch_bp(self):
        bp = fnc.BloodPressure(self.data, 10000, 1000)
        bp.run_filter()
        batches = list(bp.get_batch_series())
        self.assertEqual(len(batches), 10)
    
    def test_hrv(self):
        bp = fnc.BloodPressure(self.data, 1000, 100)
        batch_bp = next(bp.get_batch_series())
        n, generator = batch_bp.get_windows_generator()
        for window in generator():
            res = batch_bp.calc_hrv(window)
            for each in res:
                self.assertNotEqual(each, 0)
            break
    
    def test_angle(self):
        bp = fnc.BloodPressure(self.data, 1000, 100)
        batch_bp = next(bp.get_batch_series())
        start_time = batch_bp.get_start_times()[0]
        angle = batch_bp.calc_angle(start_time, 10)
        self.assertLess(angle, 90)
