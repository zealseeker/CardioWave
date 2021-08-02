#%%
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
from scipy.signal import find_peaks
from cdwave import fnc

#%%
df = pd.read_csv('peak-example.csv', index_col=0)
f, axes = plt.subplots(1,1)

series = df.set_index('time1')['signal1']
series.index.name = 'time'
wave = fnc.Waveform(series)
wave.get_peaks()
ax = axes
ax.plot(wave.df['time'], wave.df[wave.name], '-')
s = wave.df[wave.df['peak']==1]
ax.scatter(s.time, s.signal1, c='green')
xlim = ax.get_xlim()
ylim = ax.get_ylim()
first_x = wave.df.loc[wave.df.peak==1,'time'].iloc[0]
last_x = wave.df.loc[wave.df.peak==1,'time'].iloc[-1]
ax.add_patch(Rectangle((xlim[0], ylim[0]), first_x-xlim[0], ylim[1]-ylim[0],
                      alpha=0.4, facecolor='grey'))
ax.add_patch(Rectangle((last_x, ylim[0]), xlim[1]-last_x, ylim[1]-ylim[0],
                      alpha=0.4, facecolor='grey'))

f
