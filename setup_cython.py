from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("fions_cy.pyx")
)

# python setup.py build_ext --inplace (in CLI to build Cython file)
