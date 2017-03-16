# -*- coding: utf-8 -*-

"""
## Installation

 1. Download the project source or clone the repository https://github.com/ckald/pyBBN

        git clone --recursive https://github.com/ckald/pyBBN.git

 1. Install dependencies using favorite OS package manager and [pip](https://pypi.python.org/pypi/pip)  \
    (usually preinstalled on Linux and Mac OS X along with the standard CPython 2.7+ interpreter)

     * Ubuntu and aptitute: `sudo aptitude install libblas-dev liblapack-dev python-dev gfortran`
     * Mac OS X and [Homebrew](https://brew.sh/): `brew install scipy`

```
    sudo pip install -r requirements.txt
```

 1. Compile modified KAWANO nucleosynthesis code and run unit tests

        cd KAWANO/
        make all

 1. Compile code extensions and run unit tests

        make all


## Usage

To run some script in the context of the code, one should temporarily alter the `PYTHONPATH`\
environment variable:

    PYTHONPATH=. python some_script.py

For example, to run the cosmic neutrino background temperature test:

    PYTHONPATH=. python tests/cosmic_neutrino_temperature

"""
