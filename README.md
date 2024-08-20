# PyChooch

Python binding to calculate anomalous scattering factors from X-ray fluorescence spectra.

This repository adds a Python C extension to [chooch](https://github.com/fadnywg/chooch).

To compile in a conda environment:

```bash
conda create -n chooch c-compiler make cmake plplot gsl python=3.11
conda activate chooch

# To compile and install the Python binding
python -m pip install .
# To test
python test_pychooch.py CHOOCH_TEST.raw

# To build the chooch binary
cmake -S . -B build
cmake --build build
# To test
./build/chooch -e Se -a K CHOOCH_TEST.raw
```
