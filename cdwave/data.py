import os
import pickle
import gzip
import logging
import warnings
import pandas as pd
from cdwave import fnc

logger = logging.getLogger(__name__)
crt_path = os.path.abspath(os.path.dirname(__file__))

class WaveformFull:
    """WaveformFull contains all the information such as compound name,
    concentration or vendor, and all calculated parameters of a waveform.

    Attributes:
        profile: A dictionary of The profile of the waveform, including:
            plate: The plate name form which the waveform is generated.
            compound: The compound name of the waveform.
            concentration: The concentration of the compound of the waveform.
            well: The well of the waveform.
            cpid: Compound id.
        signal (dict): A dictionary with with keys, x and y, which are are time
            and calcium transient of the waveform.
        state: The state of the waveform, such as before treatment and after
            treatment.
        parameters (dict): A dictionary containing all the parameters of the
            waveform.

    Args:
        item (dict): A dictionary contains all the information of the waveform
            required keys: signals, plate, well, concentration. item['singal']
            is a dictionary with times (x) and signals (y)
            `{'x': [0.1, 0.2], 'y':[100, 101]}`

    Example:
        >>> item = {'signal': {'x': [0.1, 0.2], 'y': [100, 101]},
        ...         'plate': 'P1',
        ...         'compound': 'cmp1',
        ...         'concentration': 0.01,
        ...         'well': 'A1'}
        >>> wave = WaveformFull(item)
    """
    parameter_df = pd.read_csv(os.path.join(crt_path, 'parameters.csv'), index_col=0)
    parameter_names = parameter_df.index.tolist()
    parameter_annotations = parameter_df['Annotation'].to_dict()

    def __init__(self, item):
        if 'compound' not in item:
            item['compound'] = 'undefined'
        if 'concentration' not in item:
            item['concentration'] = 0
        if 'cpid' not in item:
            item['cpid'] = item['compound']
        self.signal = self.standardise_signal(item['signal'])
        self.parameters = {x: 0 for x in self.parameter_names}
        self.profile = {
            'plate': str(item['plate']),
            'compound': item['compound'],
            'concentration': item['concentration'],
            'well': item['well'],
            'vendor': item.get('vendor', 'None'),
            'cpid': str(item['cpid']),
            'state': item.get('state', 'treat')
        }

    def get_dict(self):
        d = {'signal': self.signal}
        d.update(self.profile)
        d['waveform'] = self
        return d

    @staticmethod
    def standardise_signal(signal: dict):
        signals = signal['y']
        if len(signals) < 5:
            warnings.warn(
                "Number of signal is less than 5 and cannot be smoothed")
            return signals
        signals = fnc.signal_filter(signals)
        minimum = signals.min()
        signal['y'] = signals - minimum
        return signal

    def get_signal_series(self):
        se = pd.Series(self.signal['y'], index=self.signal['x'])
        se.name = self.profile['well']
        se.index.name = 'time'
        return se

    def get_parameters(self, fillna=0):
        """Return a dictionary with all parameters
        Args:
            fillna: If fill na is 'raise', then an exception will be raised when
            a parameter has not been calculated. For other values, it will be used
            to fill the empty parameter.

        Returns:
            dict: A dictionary with all parameters listed in `parameter_names`
        """
        if fillna == 'raise':
            return {x: self.parameters[x] for x in self.parameter_names}
        return {x: self.parameters[x] if x in self.parameters else fillna for x in self.parameter_names}

    def __str__(self):
        return "waveform: {}".format(self.profile)


class Dataset:
    """A dataset contains a list of WaveformFull objects and a meta table
    with the information and parameters of the waveforms

    Attributes:
        waveforms: A list of WaveformFull objects.
        dataframe: The meta table of the waveforms.
        filtered_df: A view of the meta table.
        filterable_columns: columns that can be used to filter the waveforms
        size: Number of waveforms in the dataset

    Args:
        waveforms: A list of WaveformFull objects
    """

    def __init__(self, waveforms: list = None):
        self.waveforms = waveforms
        self._dataframe = None
        self.filtered_df = None
        self.filterable_columns = [
            'plate', 'compound', 'state', 'well', 'vendor', 'cpid', 'concentration']
        self.size = len(waveforms)

    def __len__(self):
        return len(self.waveforms)

    @staticmethod
    def loaddata(filepath) -> 'Dataset':
        file_items = filepath.split('.')
        if len(file_items) > 1 and file_items[-1] == 'gz' and file_items[-2] == 'pickle':
            fp = gzip.open(filepath, 'rb')
        elif file_items[-1] == 'pickle':
            fp = open(filepath, 'rb')
        else:
            return None
        dataset = pickle.load(fp)
        dataset.dataframe = dataset.get_df()
        dataset.filtered_df = dataset.dataframe
        fp.close()
        return dataset

    @staticmethod
    def concat(datasets: list):
        merged_dataset = datasets[0].copy()
        for dataset in datasets[1:]:
            merged_dataset.merge(dataset)
        return merged_dataset

    def merge(self, dataset):
        self.waveforms += dataset.waveforms

    def copy(self):
        return Dataset(self.waveforms)

    def get_df(self) -> pd.DataFrame:
        return pd.DataFrame([x.get_dict() for x in self.waveforms])

    def get_parameter_df(self) -> pd.DataFrame:
        """Return a dataframe with all papameters for each wave"""
        df = self.df[self.filterable_columns+['waveform']]
        for parameter in WaveformFull.parameter_names:
            df[parameter] = [x.parameters[parameter] for x in df.waveform]
        for p in ['noise', 'double_peak', 'uniform']:
            df[p] = df[p].astype(bool)
        return df

    def save(self, filename, compress=True):
        if compress:
            fp = gzip.open(filename + '.gz', 'wb')
        else:
            fp = open(filename, 'wb')
        pickle.dump(self, fp)
        fp.close()

    @property
    def df(self):
        if self._dataframe is None:
            self._dataframe = self.get_df()
            self.filtered_df = self._dataframe
        return self._dataframe

    @property
    def dtypes(self):
        return self.df.dtypes

    def filter_by_filters(self, filters, replace=True):
        """
        Args:
            filters: a dictionary of {column: value}
        return:
            filtered dataframe
        """
        df = self.df
        condition = [True] * len(df)
        for column in filters:
            value = filters[column]
            condition &= df[column] == value
        if replace:
            self.filtered_df = df[condition]
            return self.filtered_df
        return df[condition]


class DataLoader:

    def __init__(self, filepath: str = None, log=None):
        self.filepath = filepath
        self.dataframe = pd.DataFrame()
        self.signals = []
        self.waveforms = []
        self.pbar = None
        if filepath and not os.path.exists(filepath):
            raise Exception(
                "Cannot open the file {}, doesn't exist!".format(filepath))
        if log is None:
            self.log = logger.debug
        else:
            self.log = log

    def transfer(self) -> Dataset:
        """Parse the waveform from tables and generate a dataset

        Returns:
            Dataset: The dataset containing all the waveforms and parameters
        """
        self.log('Start parsing....')
        self.parse()
        self.log('Parsing complete, {} waveforms in total'.format(
            len(self.waveforms)))
        dataset = Dataset(self.waveforms)
        return dataset

    def parse(self):
        """This function should parse waveform data from files into self.waveforms"""
        raise NotImplementedError


class StandardCSVLoader(DataLoader):
    """A loader parsing "standard csv file". The format of the table file is like
    ```
    compound,concentration,well,plate,time,signal
    CP1,0.1,A1,P1,0,1000
    CP1,0.1,A1,P1,0.33,1001
    CP2,0.1,A2,P1,0,1000
    ```
    """

    def __init__(self, filepath: str = None, data: pd.DataFrame = None, log=None):
        if isinstance(filepath, pd.DataFrame):
            data = filepath
            filepath = None
        super().__init__(filepath=filepath, log=log)
        if filepath is None:
            if data is None:
                raise ValueError('File path and data cannot both be None')
            self.data = data

    def parse(self):
        if self.filepath is None:
            df = self.data
        else:
            df = pd.read_csv(self.filepath)
        required_columns = ['plate', 'compound',
                            'concentration', 'well', 'time', 'signal']
        attribute_columns = ['plate', 'compound',
                             'concentration', 'well', 'cpid', 'vendor']
        miss_column = set(required_columns) - set(df.columns)
        if miss_column:
            raise KeyError('Column {} is missing'.format(miss_column))
        if df.duplicated(['plate', 'well', 'time']).any():
            raise ValueError('Find duplicacy in set')
        for _, gdf in df.groupby(['plate', 'well']):
            item = gdf.iloc[0]
            input_item = {x: item[x] for x in attribute_columns if x in item}
            signal = {'x': gdf['time'].values, 'y': gdf['signal'].values}
            input_item['signal'] = signal
            waveform = WaveformFull(input_item)
            self.waveforms.append(waveform)


if __name__ == "__main__":
    print('Number of parameters: {}'.format(len(WaveformFull.parameter_names)))
