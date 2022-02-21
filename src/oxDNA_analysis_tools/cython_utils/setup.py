from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    include_dirs=[np.get_include()],
    ext_modules = cythonize("parse_conf.pyx", annotate=True, compiler_directives={'language_level' : "3"})
)