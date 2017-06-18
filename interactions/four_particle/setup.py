""" Cython extension setup.py file.

    To build an extension, run

        cd <extension_folder>
        python3 setup.py build_ext --inplace
"""

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sourcefiles = ['integral.pyx']

extensions = [Extension("integral", sourcefiles,
                        include_dirs=[numpy.get_include()],
                        language="c++",
                        libraries=["stdc++"]
                        )]

setup(
    cmdclass={'build_ext': build_ext},
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize(extensions)
)
