# -*- coding: utf-8 -*-

"""
## Installation

 1. Clone the repository https://github.com/ckald/pyBBN

        git clone --recursive https://github.com/ckald/pyBBN.git


 1. Install dependencies using favorite OS package manager and \
    [pip](https://pypi.python.org/pypi/pip) \
    (usually preinstalled on Linux and Mac OS X along with the standard CPython interpreter)

     * Ubuntu and aptitute: `sudo aptitude install libblas-dev liblapack-dev python-dev gfortran`
     * Mac OS X and [Homebrew](https://brew.sh/): `brew install scipy`

    Currently the code supports only Python3.5+. In case the default `python3` command points to an older version, \
    install Python3.5+ and use the exact version (e.g., `python3.5` or `python3.6`) everywhere instead of `python3`.

```
    python3 -m pip install -r requirements.txt
```
    
    (Installation of NumPy/SciPy scientific libraries might require `sudo`)

 1. Compile modified KAWANO nucleosynthesis code and run unit tests

        cd KAWANO/
        make all

 1. Compile code extensions and run unit tests

        PYTHON=python3 make all


## Usage

To run some script in the context of the code, one should temporarily alter the `PYTHONPATH`\
environment variable:

    PYTHONPATH=. python3 some_script.py

For example, to run the cosmic neutrino background temperature test:

    PYTHONPATH=. python3 tests/cosmic_neutrino_temperature

"""
