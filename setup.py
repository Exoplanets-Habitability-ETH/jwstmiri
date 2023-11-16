#!/usr/bin/env python
from setuptools import setup, find_packages


setup(
    name='jwstmiri',
    version='0.0.1',
    packages=['data', 'plots', 'tests', 'pipeline'],
    url='https://github.com/Exoplanets-Habitability-ETH/jwstmiri/',
    license='MIT',
    author='Polychronis Patapis',
    author_email='patapisp@ethz.ch',
    description='Toolkit for analysis of jwst exoplanet data',
    install_requires=["jwst >= 1.12.5"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)

