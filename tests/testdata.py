import unittest
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
        self.assertFalse(default_more, "These parameters are not calculated}")

def main():
    dataset = data.Dataset.loaddata(r'C:\Users\hy929891\OneDrive - University of Cambridge\Projects\Cardiotox\WaveformData\alldata.pickle.gz')
    wave = dataset.waveforms[1622]
    waveform = fnc.Waveform(wave.get_signal_series())
    waveform.get_peaks()
    waveform.analyse()
    waveform.get_parameters()
