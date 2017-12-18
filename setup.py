
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import sys

# to run
# py setup.py build_ext --inplace
# py setup.py install

# setup bed parser
ext_modules = [Extension("mscentipede", sources=["mscentipede.pyx"])]
ext_modules = cythonize(ext_modules)

setup(
    name = 'msCentipede',
    author = 'Anil Raj',
    version = '1.0',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)
