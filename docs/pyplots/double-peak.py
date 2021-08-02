#%%
from matplotlib import pyplot as plt
import pandas as pd
from cdwave import fnc
import numpy as np

def get_tails(wave):
    group = wave.df.groupby('group')
    tails = []
    for i, gdf in group:
        if i in (0, wave.num_peak):
            continue
        starting_point = gdf[gdf['status'] == 1]
        downing_point = gdf[gdf['status'] == 3]
        if len(downing_point) == 0:
            return False
        if len(starting_point) == 0:
            tails.append(0)
        else:
            tails.append(
                starting_point.iloc[0]['time'] - downing_point.iloc[0]['time'])
    tails = np.array(tails)
    return tails
#%%
waveforms = {}
times = []
dataname = ''
with open('../../tests/test_data.csv') as fp:
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

f, axes = plt.subplots(1,2, figsize=(8, 4), sharey='all')
series = waveforms['double_peak'].set_index('time')['F6']
series.index.name = 'time'
wave = fnc.Waveform(series)
wave.get_peaks()
wave.get_valleys()
wave.analyse_normal_waves()
tails = get_tails(wave)
df = wave.df
ax = axes[0]
ax.plot(series.index, series, '-')
colors = { 1: 'orange', 2: 'red', 3: 'green'}
for color in colors:
    tdf = df[df['status'] == color]
    ax.plot(tdf['time'], tdf[series.name],
                'o', color=colors[color])
tdf = df[(df['peak'] > 1)]
ax.plot(tdf['time'], tdf[series.name], 'mo')

ax = axes[1]
wave.get_peaks()
wave.analyse()
ax.plot(series.index, series, '-')
colors = { 1: 'orange', 2: 'red', 3: 'green'}
for color in colors:
    tdf = df[df['status'] == color]
    ax.plot(tdf['time'], tdf[series.name],
                'o', color=colors[color])
tdf = df[(df['peak'] > 1)]
ax.plot(tdf['time'], tdf[series.name], 'mo')

f

# %%
