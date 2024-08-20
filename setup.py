from setuptools import setup, Extension
import os

SRC_PATH = os.path.join(os.path.dirname(__file__))
CHOOCH_SOURCE_FILES = [
    "printbanner.c", "minmax.c", "mucal.c", "fdprime.c",
    "smooth.c", "fits.c", "normalize.c", "checks.c", "integrate.c",
    "selwavel.c", "toplot.c", "savwin.c"]

CHOOCH_SOURCE_FILES = [os.path.join(SRC_PATH, x) for x in CHOOCH_SOURCE_FILES]
PYCHOOCH_SOURCE_FILES = ["PyChooch.c"] + CHOOCH_SOURCE_FILES

setup(
  name="PyChooch",
  version="5.0.9",
  author=("Chooch: Gywndalf Evans, Diamond Light Source;"
          "Python binding: Matias Guijarro, ESRF"),
  ext_modules=[Extension(
    'PyChooch', PYCHOOCH_SOURCE_FILES, libraries=['gsl', 'gslcblas', 'plplot'])],
)
