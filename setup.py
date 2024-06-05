from setuptools import setup
from unWISExLens_lklh import __author__, __version__
import os

file_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(file_dir)

setup(name="unWISExLens_lklh",
      version=__version__,
      author=__author__,
      license='BSD 2-Clause',
      description='Likelihood software for unWISE cross-correlation analysis with CMB lensing',
      zip_safe=False,  # set to false if you want to easily access bundled package data files
      packages=['unWISExLens_lklh'],
      package_data={'unWISExLens_lklh': ['*.yaml', '*.bibtex', 'data/*', 'data/**/*']},
      # tests_require=['camb>=1.0.5']
      )