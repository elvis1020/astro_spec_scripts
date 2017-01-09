#!/usr/bin/python

from setuptools import setup, find_packages
from glob import glob

setup(
    name = 'astro_spec_scripts',
    packages = find_packages(),
    version = '0.17.01.09',
    license = 'GNU GPLv3',
    platforms = 'any',
    description = 'Tools for Astronomy; Spectroscopy; ESO data reduction; Spectral plotting; Spectral synthesis; Atomic linelists',
    author = 'Elvis Cantelli',
    author_email = 'elvis.shock@gmail.com',
    url = 'https://github.com/elvis1020/astro_spec_scripts', # use the URL to the github repo
    keywords= ['Tools for Astronomy', 'Spectroscopy', 'ESO data reduction', 'Spectral plotting', 'Spectral synthesis', 'Atomic linelists'],
    install_requires = ['numpy', 'scipy', 'astropy', 'matplotlib', 'PyAstronomy', 'pyraf', 'lmfit'], 
    scripts = glob('plotting/*.py')
)
