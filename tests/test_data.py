import unittest
from io import StringIO
import pandas as pd
from cdwave import data, fnc
from tests import Dataset


class TestData(unittest.TestCase):
    _data = Dataset()

    def get_mock_waveform(self):
        df = self._data['D1']
        input_item = {
            'plate': 'plate1',
            'compound': 'compound1',
            'concentration': 100,
            'cpid': 'x',
            'vendor': 'AstraZeneca',
            'well': 'B14',
            'state': 'prior'
        }
        input_item['signal'] = {
            'x': df['time'].values,
            'y': df['B14'].values
        }
        waveform = data.WaveformFull(input_item)
        return waveform

    def test_parameter_consistency(self):
        waveform = self.get_mock_waveform()
        series = waveform.get_signal_series()
        wave = fnc.Waveform(series)
        wave.get_peaks()
        wave.analyse()
        parameters = wave.get_parameters()
        parameters_more = set(parameters.keys()) - \
            set(data.WaveformFull.parameter_names)
        self.assertFalse(
            parameters_more, "These parameters are not included in default parameter names")
        default_more = set(data.WaveformFull.parameter_names) - \
            set(parameters.keys())
        self.assertFalse(default_more, "These parameters are not calculated")

    def test_csvloader(self):
        csv_string = """compound,concentration,well,plate,time,signal
        CP1,0.1,A1,P1,0,1000
        CP1,0.1,A1,P1,0.33,1001
        CP2,0.1,A2,P1,0,1000""".replace(' ', '')
        obj = StringIO(csv_string)
        df = pd.read_csv(obj)
        dataset = data.StandardCSVLoader(df).transfer()
        self.assertEqual(len(dataset), 2)
        obj = StringIO(csv_string+'\nCP1,0.1,A1,P1,0,100')
        df = pd.read_csv(obj)
        with self.assertRaises(ValueError):
            dataset = data.StandardCSVLoader(df).transfer()
