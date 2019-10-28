# setup.py
# --------

from setuptools   import setup, find_packages, Extension

import numpy
import version

try:
   import setup_cleaner
   CleanCommand = setup_cleaner.CleanCommand
except ImportError:
   CleanCommand = None

cmd_classes = {}

if CleanCommand is not None:
   cmd_classes['clean'] = CleanCommand

# Define cython extension(s)

#src_cohen  = ["tftb/tfr/cohen.c"]
#src_linear = ['tftb/tfr/linear.c']
#extensions = [Extension("tftb/tfr/linear", src_linear,
#                        include_dirs = [numpy.get_include()]),
#              Extension("tftb/tfr/cohen", src_cohen,
#                                      include_dirs = [numpy.get_include()]),
#              ]

extensions = None

setup(name = 'fourierin4py',
      version = version.get_version(),
      description = "Computing Fourier integrals in Python",
      long_description = """A port of R's fourierin package to python.""",
      classifiers = [], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords = '',
      author = 'Christian Marquardt',
      author_email = 'christian@marquardt.sc',
      url = 'www.marquardt.sc',
      license = 'GPL',
      packages = find_packages(exclude = ['examples', 'tests']),
      include_package_data = True,
      zip_safe = False,
      setup_requires = ['pytest'],
      cmdclass = cmd_classes,
      ext_modules = extensions,
      install_requires = [
          # -*- Extra requirements: -*-
      ],
      entry_points = {
          # normal parameters, ie. console_scripts[]
          'distutils.commands': [
              ' clean = setup_cleaner:CleanCommand']
      }
)