# -*- coding: utf-8 -*-

"""
## Installation

 1. Download the project source or clone the repository https://github.com/ckald/pyBBN

        git clone https://github.com/ckald/pyBBN.git



 1. Install dependencies using OS package manager and [pip](https://pypi.python.org/pypi/pip)  \
    (usually preinstalled on Linux and Mac OS X along with the standard CPython 2.7+ interpreter)

     * Ubuntu: `sudo aptitude install libblas-dev liblapack-dev python-dev gfortran`
     * Mac OS X: `brew install scipy`

```
    sudo pip install -r requirements.txt
```

## Usage

To run some script in the context of the code, one should temporarily alter the `PYTHONPATH`\
environment variable:

    PYTHONPATH=. python some_script.py

For example, to run the cosmic neutrino background temperature test:

    PYTHONPATH=. python tests/cosmic_neutrino_temperature/__init__.py

"""
