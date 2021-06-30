*	Uniformity (`uniform`): Peak uniformity (p > 0.99 in KS test)
*	Noisy Waveform (`noise`): If the waveform is noise, which means it is probabily a bad waveform
*	Double Peak (`double_peak`): If there is double peak (max_combo_peaks > 1)
*	Fail (`fail_analysis`): If there is a missing key point (starting point, down point etc.)
*	σ(Amplitude) (`std_amplitude`): The standard deviation of amplitudes
*	σ(Intensity) (`std_intensity`): The standard deviation of intensity
*	σ(Peak Space) (`std_lambda`): The standard deviation of period
*	σ(Inner Peak Space) (`std_inner_lambda`): The standard deviation of period removing outliers
*	σ(Shoulder Position) (`std_shoulder`): The standard deviation of shoulder position
*	σ(Shoulder/Tail) (`std_shoulder_tail`): The standard deviation of shoulder tail ratio
*	σ(PW10) (`PW10_std`): Standard deviation of the peak widths at 10% of prominence
*	σ(PW25) (`PW25_std`): Standard deviation of the peak widths at 25% of prominence
*	σ(PW50) (`PW50_std`): Standard deviation of the peak widths at 50% of prominence
*	σ(PW80) (`PW80_std`): Standard deviation of the peak widths at 80% of prominence
*	σ(PW90) (`PW90_std`): Standard deviation of the peak widths at 90% of prominence
*	σ(Tail Duration) (`std_tail`): Standard deviation of Tail Duration
*	σ(Shoulder Amplitude) (`std_shoulder_amp`): The standard deviation of amplitude of shoulder
*	Average Rise/Decay (`rd_ratio`): Average rise time/decay time ratio
*	Average Tail Proportion (`tail_proportion`): The average proportion of tail in one period
*	Average Shoulder/Tail (`avg_shoulder_tail`): The average ratio between prominence of shoulder and tail in a kernel density distribution
*	FFT Ratio (`fft_ratio`): Ratio between second and first peak position in Fast Fourier Transformation
*	Average Rise Time (`up_length`): Average rise time
*	Average Decay Time (`down_length`): Average duration of decreasing waveform, from peak point to the first point lower than valley value plus variance (10% of max)
*	Peak To End (`full_down`): Average of duration between peak and next starting point
*	Number of Peaks (`n_peak`): Number of peaks
*	Peak Frequency (`freq`): Number of period for the 100 seconds, 100/Peak Space
*	Maximum Amplitude (`maximum`): The highest amplitude
*	Maximum Intensity (`max_intensity`): The highest signal in the waveform
*	Average Peak Amplitude (`avg_amplitude`): Average peak amplitude
*	Average Intensity (`avg_intensity`): The average intensity (highest signal in a period)
*	RMS (`rms`): Root mean square of peak amplitudes
*	Average Peak Space (`avg_lambda`): The average duration of period
*	Average Inner Peak Space (`avg_inner_lambda`): The average duration of periods removing outliers
*	Average Tail Duration (`avg_tail`): The average duration of tail
*	Average Valley (`avg_valley`): The average intensity of valley points
*	Average Shoulder Position (`avg_shoulder`): The average relative position of shoulder in a period (0-1)
*	Multi-Peaks (`max_combo_peaks`): Maximum of the number of peaks in one period
*	PW10 (`PW10_mean`): The average of peak width at 10% from peak to bottom.
*	PW25 (`PW25_mean`): The average of peak width at 25% from peak to bottom
*	PW50 (`PW50_mean`): The average of peak width at 50% from peak to bottom
*	PW80 (`PW80_mean`): The average of peak width at 80% from peak to bottom
*	PW90 (`PW90_mean`): The average of peak width at 90% from peak to bottom
*	Average Shoulder Amplitude (`avg_shoulder_amp`): The average of amplitude of shoulder
*	Minimum Intensity (`min_intensity`): The minimum signal of the waveform, should be 0 if scaled
