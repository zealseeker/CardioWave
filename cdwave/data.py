# Copyright (C) 2021 by University of Cambridge

# This software and algorithm was developed as part of the Cambridge Alliance
# for Medicines Safety (CAMS) initiative, funded by AstraZeneca and
# GlaxoSmithKline

# This program is made available under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the
# License, or at your option, any later version.
import os
import pickle
import gzip
import logging
import warnings
from itertools import product
from typing import List
import pandas as pd
from cdwave import fnc

logger = logging.getLogger(__name__)
crt_path = os.path.abspath(os.path.dirname(__file__))


class WaveformFull:
    """WaveformFull contains all the information about a sample, such as
    compound name, concentration or vendor, and all calculated parameters of a
    waveform.

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
        scale (bool): Whether to scale the minimum to 0. Default is True
        window_length (int): The length of the filter window (i.e. the number of
            coefficients). window_length must be a positive odd integer. If set
            to 0, the waveform will not be smoothed.

    Example:
        >>> item = {'signal': {'x': [0.1, 0.2], 'y': [100, 101]},
        ...         'plate': 'P1',
        ...         'compound': 'cmp1',
        ...         'concentration': 0.01,
        ...         'well': 'A1'}
        >>> wave = WaveformFull(item)
    """
    parameter_df = pd.read_csv(os.path.join(
        crt_path, 'parameters.csv'), index_col=0)
    parameter_names = parameter_df.index.tolist()
    parameter_annotations = parameter_df['Annotation'].to_dict()

    def __init__(self, item, scale=True, window_length=5):
        if 'compound' not in item:
            item['compound'] = 'undefined'
        if 'concentration' not in item:
            item['concentration'] = 0
        if 'cpid' not in item:
            item['cpid'] = item['compound']
        self.signal = self.standardise_signal(
            item['signal'], scale, window_length)
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
    def standardise_signal(signal: dict, scale=True, window_length=5):
        signals = signal['y']
        minimum = signals.min()
        if len(signals) < 5:
            warnings.warn(
                "Number of signal is less than 5 and cannot be smoothed")
            return signals
        if window_length != 0:
            signals = fnc.signal_filter(signals, window_length)
        if scale:
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

    def __init__(self, waveforms: List[WaveformFull] = None):
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
        dataset: 'Dataset' = pickle.load(fp)
        dataset.dataframe = dataset.get_df()
        dataset.filtered_df = dataset.dataframe
        fp.close()
        return dataset

    @staticmethod
    def concat(datasets: list) -> 'Dataset':
        merged_dataset = datasets[0].copy()
        for dataset in datasets[1:]:
            merged_dataset.merge(dataset)
        return merged_dataset

    def merge(self, dataset):
        self.waveforms += dataset.waveforms
        self.size += len(dataset.waveforms)

    def copy(self):
        return Dataset(self.waveforms.copy())

    def get_df(self) -> pd.DataFrame:
        return pd.DataFrame([x.get_dict() for x in self.waveforms])

    def get_parameter_df(self) -> pd.DataFrame:
        """Return a dataframe with all papameters for each wave"""
        df = self.df[self.filterable_columns+['waveform']].copy()
        parameters = WaveformFull.parameter_names
        missing_parameters = set(parameters) - \
            set(df.waveform.iloc[0].parameters)
        if len(missing_parameters) > 0:
            warnings.warn("Parameters {} cannot be found in the dataset."
                          "This dataset may be generated by previous version".format(
                              ','.join(missing_parameters)))
            parameters = df.waveform.iloc[0].parameters
        for parameter in parameters:
            df[parameter] = [x.parameters[parameter] for x in df.waveform]
        for p in ['noise', 'double_peak', 'uniform']:
            df[p] = df[p].astype(bool)
        return df

    def save(self, filename, compress=True):
        # TODO remove _dataframe to save space
        if compress:
            if filename[-10:] != '.pickle.gz':
                filename += '.pickle.gz'
            fp = gzip.open(filename, 'wb')
        else:
            fp = open(filename, 'wb')
        pickle.dump(self, fp)
        fp.close()
    
    def export_raw(self, filename=None, compression='infer'):
        """Export the row data into a csv file with the columns of
        compound,concentration,well,plate,time,signal
        """
        rows = []
        for waveform in self.waveforms:
            profile = waveform.profile
            times = waveform.signal['x']
            signals = waveform.signal['y']
            for t, s in zip(times, signals):
                rows.append({
                    'compound': profile['compound'],
                    'well': profile['well'],
                    'plate': profile['plate'],
                    'concentration': profile['concentration'],
                    'state': profile['state'],
                    'time': t,
                    'signal': s
                })
        df = pd.DataFrame(rows)
        if filename:
            df.to_csv(filename, compression=compression)
            return None
        return df

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
        self.signals = []
        self.waveforms: list[WaveformFull] = []
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
    """A loader parsing "standard csv file".

    The format of the table file is like::
        
        compound,concentration,well,plate,time,signal
        CP1,0.1,A1,P1,0,1000
        CP1,0.1,A1,P1,0.33,1001
        CP2,0.1,A2,P1,0,1000
    
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
        attribute_columns = ['plate', 'compound', 'state',
                             'concentration', 'well', 'cpid', 'vendor']
        miss_column = set(required_columns) - set(df.columns)
        if miss_column:
            raise KeyError('Column {} is missing'.format(miss_column))
        if 'state' not in df.columns:
            df['state'] = 'treat'
        if df.duplicated(['plate', 'well', 'time', 'state']).any():
            raise ValueError('Find duplicacy in set')
        for _, gdf in df.groupby(['plate', 'well', 'state']):
            item = gdf.iloc[0]
            input_item = {x: item[x] for x in attribute_columns if x in item}
            signal = {'x': gdf['time'].values, 'y': gdf['signal'].values}
            input_item['signal'] = signal
            waveform = WaveformFull(input_item)
            self.waveforms.append(waveform)


class SeqLoader(DataLoader):
    def __init__(self, filepath: str, plate: str = None, state: str = None, opts=None, log=None):
        filename = os.path.basename(filepath)
        assert filename[-5:] == '.seq1'
        if plate is None and state is None:
            plate, state = filename.split('.')[0].split('_')
        elif plate is None or state is None:
            raise ValueError(
                "Plate and state should be both None or neither None")
        self.plate = plate
        self.state = state
        self.opts = {'cut': {'pre': 350, 'post': None}
                     } if opts is None else opts
        super().__init__(filepath, log)

    def parse(self):
        df = pd.read_csv(self.filepath, sep='\t')
        time = df.iloc[:, 4].values
        state = self.state
        waveforms = []
        for well in get_wells():
            input_item = {
                'plate': self.plate,
                'well': well,
                'state': state,
                'signal': {
                    'x': time,
                    'y': df[well].values
                }
            }
            cut_opts = self.opts['cut']
            if state in cut_opts and cut_opts[state] is not None:
                cut = cut_opts[state]
                input_item['signal'] = {
                    'x': time[:cut], 'y': df[well].values[:cut]}
            waveform = WaveformFull(input_item)
            waveforms.append(waveform)
        self.waveforms = waveforms

    def transfer_with_meta(self, df: pd.DataFrame, columns=None) -> Dataset:
        assert len(df) == len(df.well.unique())
        df = df.set_index('well')
        if columns is None:
            columns = df.columns.tolist()
        self.parse()
        waveforms = []
        for waveform in self.waveforms:
            well = waveform.profile['well']
            if well in df.index:
                for col in columns:
                    waveform.profile[col] = df.at[well, col]
                waveforms.append(waveform)
        dataset = Dataset(waveforms)
        return dataset


def get_wells():
    wells = product(list('ABCDEFGHIJKLMNOP'), range(1, 25))
    wells = ['{}{}'.format(x, y) for x, y in wells]
    return wells


if __name__ == "__main__":
    print('Number of parameters: {}'.format(len(WaveformFull.parameter_names)))
