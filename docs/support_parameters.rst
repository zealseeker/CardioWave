*	Rise Time (`up_length`): Average rise time
*	Peak to End (`down_length`): Average duration from peak point to the tail, the first point lower than valley value plus variance (10% of max)
*	Decay Time (`full_down`): Average of duration between peak and next starting point
*	Rise/Decay (`rd_ratio`): Average rise time/decay time ratio
*	Peak Frequency (`freq`): Average number of peaks for the 100 seconds, 100/Peak Space
*	Maximum Amplitude (`maximum`): The highest amplitude
*	Average Peak Amplitude (`avg_amplitude`): Average peak amplitude
*	σ(Amplitude) (`std_amplitude`): Standard deviation of peak amplitude
*	Average Inner Peak Space (`avg_inner_lambda`): Average inner peak space, where inner means excluding cycles where the peak space is an outlier
*	Multi-Peak (`max_combo_peaks`): Maximum of multi-Peak, which is the number of peaks in each cycle
*	Average Peak Space (`avg_lambda`): Average peak space, the duration between two peaks
*	σ(Peak Space) (`std_lambda`): Standard deviation of peak space
*	σ(Inner Peak Space) (`std_inner_lambda`): Standard deviation of inner peak space
*	Average Tail Duration (`avg_tail`): Average tail duration, from the first point lower than valley value plus variance (10% of max) to the next starting point
*	σ(Tail Duration) (`std_tail`): Standard deviation of Tail Duration
*	Average Tail Proportion (`tail_proportion`): Average tail proportion, the ratio between tail druation and the peak space
*	Average Shoulder Position (`avg_shoulder`): Average shoulder position, a normalised (0-1) float indicating the position of shoulder from peak to valley
*	Average Shoulder/Tail (`avg_shoulder_tail`): Average shoulder tail ratio, which is based on kernel density estimation
*	σ(Shoulder Position) (`std_shoulder`): Standard deviation of shoulder position
*	FFT Ratio (`fft_ratio`): Fast Fourier transformation Ratio
*	σ(Shoulder/Tail) (`std_shoulder_tail`): Standard deviation of shoulder tail ratio
*	Average Intensity (`avg_intensity`): Average difference between the peak amplitude and valley amplitude in each cycle
*	σ(Average Intensity) (`std_intensity`): Standard deviation of average intensity
*	Valley (`avg_valley`): Average valley amplitude
*	Peak Count (`n_peak`): Number of peaks in the waveform
*	RMS (`rms`): Root mean square of the calcium transients
*	PW10 (`PW10_mean`): Average Peak width at 10% of prominence
*	σ(PW10) (`PW10_std`): Standard deviation of peak width at 10% of prominence
*	PW25 (`PW25_mean`): Average Peak width at 25% of prominence
*	σ(PW25) (`PW25_std`): Standard deviation of peak width at 25% of prominence
*	PW50 (`PW50_mean`): Average Peak width at 50% of prominence
*	σ(PW50) (`PW50_std`): Standard deviation of peak width at 50% of prominence
*	PW80 (`PW80_mean`): Average Peak width at 80% of prominence
*	σ(PW80) (`PW80_std`): Standard deviation of peak width at 80% of prominence
*	PW90 (`PW90_mean`): Average Peak width at 90% of prominence
*	σ(PW90) (`PW90_std`): Standard deviation of peak width at 90% of prominence
