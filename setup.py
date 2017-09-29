from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Levin sim',
  ext_modules = cythonize('levin.pyx'),
)
