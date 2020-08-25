from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules=[ Extension("SDFT",
              ["SDFT.pyx"],  # .pyx file name
              libraries=["m"],  # include m library
              extra_compile_args = ["-ffast-math"])]

setup(
  name = "SDFT",  # create .c and .so file named SDFT
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules,
  include_dirs=[numpy.get_include()])  # cython can include numpy
