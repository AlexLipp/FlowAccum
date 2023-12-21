# To install run: python setup.py build_ext --inplace
from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
import numpy

extensions = [
    Extension(
        "cfuncs",
        ["cfuncs.pyx"],
        language="c++",  # Use C++ compiler
    )
]

setup(ext_modules=cythonize(extensions), include_dirs=[numpy.get_include()])
