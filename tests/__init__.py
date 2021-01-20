import pandas as pd
from cdwave import fnc


class Dataset:
    _data = None

    def __init__(self):
        if Dataset._data is None:
            Dataset._data = self.load_data()
        self.data = Dataset._data

    def get_series(self, dataname):
        df = self.data[dataname]
        series = df.set_index('time').iloc[:, 0]
        return series

    def get_waveform(self, dataname):
        series = self.get_series(dataname)
        wave = fnc.Waveform(series)
        return wave

    @staticmethod
    def load_data() -> dict:
        times = []
        dataname = ''
        waveforms = {}
        with open('tests/test_data.csv') as fp:
            for i, line in enumerate(fp):
                items = line.strip().split(',')
                if i % 2 == 0:
                    times = items[1:]
                    dataname = items[0]
                if i % 2 == 1:
                    waveforms[dataname] = pd.DataFrame({
                        'time': times,
                        items[0]: items[1:]
                    }, dtype=float)
        return waveforms

    def __getitem__(self, dataname):
        return self.data[dataname]
