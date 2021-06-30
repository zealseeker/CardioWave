Setup and Installation
==========================================================================================

Development and documentation occurs on GitHub_.

It also has the following dependencies:

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

- numpy>=1.16
- scipy>=1.2
- tqdm>=4.32
- pandas>=0.24
- statsmodels>=0.10.2

For GUI support and visualisation analysis

- matplotlib>=3.1
- pyqt5>=5.9


Install with Conda
~~~~~~~~~~~~~~~~~~

1. Download and install miniconda for Python 3.8 from https://docs.conda.io/en/latest/miniconda.html

2. Create envrionment (optional):

- Open terminal in Mac/Linux and run `conda create -n cardiowave`

- Run `conda activate cardiowave`

3. Install via Pypi:

- Run `pip install cardiowave`


.. include:: substitutions.rst

Install with source code
~~~~~~~~~~~~~~~~~~~~~~~~

1. Download source code from GitHub_
and change directory to the LICENSE files

.. code-block:: shell-session

    .
    ├── cdwave
    ├── tests
    ├── setup.py
    └── LICENSE

2. Install dependencies by running `pip install -r requirements`

3. Run `python setup.py install`
