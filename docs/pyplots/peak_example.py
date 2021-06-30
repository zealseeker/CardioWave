from matplotlib import pyplot as plt
import pandas as pd
from scipy.signal import find_peaks
from cdwave import fnc

df = pd.read_csv('peak-example.csv', index_col=0)
f, axes = plt.subplots(1,2, figsize=(8, 4), sharey='all')

# left plot
series = df.set_index('time1')['signal1']
opt_series = series - series.index * 0.0001
span = series.max() - series.min()
peaks, properties = find_peaks(opt_series, prominence=0)
prominences = properties['prominences']
real_peaks, false_peaks = [], []
times = series.index
for i in range(len(peaks)):
    peak = peaks[i]
    if prominences[i] > 0.1 * span:
        real_peaks.append((times[peak], series.iloc[peak]))
    else:
        false_peaks.append((times[peak], series.iloc[peak]))
ax = axes[0]
ax.plot(series.index, series.values, '-')
ax.scatter(*zip(*false_peaks), c='blue')
ax.scatter(*zip(*real_peaks), c='green')
ax.set_xlim(0, 10)
ax.set_xlabel('Time(s)')
ax.set_title('Example 1')
ax.set_ylabel('RFU')

# right plot
series = df.set_index('time2')['signal2']
series.index.name = 'time'
wave = fnc.Waveform(series)
wave.get_peaks()
ax = axes[1]
ax.plot(wave.df['time'], wave.df[wave.name], '-')
s = wave.df[wave.df['peak']==1]
ax.scatter(s.time, s.signal2, c='green')
ax.set_xlabel('Time(s)')
ax.set_title('Example 2')

plt.show()
