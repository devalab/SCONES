# DataCollection

Scripts to collect predictions from methods and store them in a CSV file. Run `fetch_METHODNAME.py` to collect predictions CSV.

**NOTE:** The results CSV files shared may not be complete.

# Installation

1. `conda create --name scones_collection python=3.6`
2. `conda activate scones_collection`
3. `conda install -c conda-forge selenium`
4. Download geckodriver and specify the path to the binary in `DataFetch.DRIVER` (`DataFetch.py`)
5. `conda install numpy pandas`