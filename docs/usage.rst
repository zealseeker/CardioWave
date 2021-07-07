Usage and Examples
==================

Getting Started
~~~~~~~~~~~~~~~


Prepare a CSV table like this format:

======== ============= ==== ===== ========== ==== ======
compound concentration well plate ..others.. time signal
======== ============= ==== ===== ========== ==== ======
CP1      0.1           A1   P1    ...        0    1000
CP1      0.1           A1   P1    ...        0.33 1001
...      ...           ...  ...   ...        ...  ...
CP2      0.1           A2   P1    ...        0    1000
...      ...           ...  ...   ...        ...  ...
======== ============= ==== ===== ========== ==== ======

The order of the rows and columns do not need to be fixed but the column names
must be exactly the same to the required (e.g. lowercase). The following columns
are compursory: 'plate', 'compound', 'concentration', 'well', 'time', 'signal'.
Optional columns include 'cpid' (compound ID) and 'vendor'. Other columns will
not be used.

.. code-block:: python

    import pandas as pd
    from cdwave import data
    from cdwave import derive

    # Load and convert
    df = pd.read_csv('data.csv')
    loader = data.StandardCSVLoader(data=df)
    dataset = loader.transfer()

    # Calculate parameters
    derive.calc_parameters_for_waveforms(dataset)
    dataset.save('data.pickle.gz')

    # Export parameters
    df = dataset.get_parameter_df()
    df.to_csv(os.path.join(data_path, 'parameters.csv'))

Advanced Usage
~~~~~~~~~~~~~~

There are two ways to derive parameters, with and without a threshold.
`calc_parameters_for_waveforms <source/cdwave.html#cdwave.derive.calc_parameters_for_waveforms>`_ and
`calc_parameters_with_threshold <source/cdwave.html#cdwave.derive.calc_parameters_with_threshold>`_.
In either method, you can define a customed function to derive parameters from the waveform data.


Derive parameters using a customed calculator
---------------------------------------------

The function should include these steps:

1. Get signal series from waveform data to instantiate `Waveform`

.. code-block:: python

    series = waveform.get_signal_series()
    wave = Waveform(series)

2. Detect peaks

.. code-block:: python
    
    wave.get_peaks()

3. Identify key time points

During this step, double-peaks or multi-peaks will be identified.

.. code-block:: python

    wave.analyse()

4. Derive parameters

By default, use `r = wave.get_parameters()` to calculate all parameters supported.
If only peak count and amplitudes are required, we can just use:

.. code-block:: python

    r = {'maximum': wave.maximum, 'n_peak': wave.num_peak}
    r.update(wave.calc_amplitudes())

The following function is the default calculator in CardioWave for derving
parameters with a threshold.

.. code-block:: python

    def calc_parameter_with_threshold(waveform: WaveformFull, threshold) -> dict:
        series = waveform.get_signal_series()
        wave = Waveform(series)
        if not wave.get_peaks(height=threshold):
            return None
        wave.analyse()
        try:
            r = wave.get_parameters()
        except Exception as e:
            logging.error('Cannot get parameters of %s', waveform)
            raise e
        return r


Selection of hyper-parameters
-----------------------------

The hyper-parameters will affect the detection of the peaks and subpeaks.
In peak detect `get_peaks`, either prominence and height can be set to ignore
"false peaks". If prominence is None, the prominence will be set to
`max(min_prominence, span_ratio * span)`. Span is the difference between the
highest intensity and the lowest intensity. The default span_ratio is 0.1 and
the default min_prominence is 20.

.. plot:: pyplots/peak_example.py

In the left example, the peaks before 6s are all false peaks (red points),
which may caused by noises. 
Their prominences are all lower than 10% of the highest itensity.
The green points are real points as their prominences are far higher than the
threshold.

In the right example, althrough these peaks are real in terms of the pre-set
rule of prominence higher than 0.1 * span, we know that they are actually false
peaks as their signals are very low compared to the left example. However,
CardioWave does not know that. So, one potential solution is to manually set the
threshold by `get_peak(prominence=100)`, where the threshold 100 depends on
the data we have. We can set the threshold according to positive control.
