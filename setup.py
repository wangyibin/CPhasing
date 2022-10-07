#!/usr/bin/env python
# -*- coding:utf-8 -*-


import os.path as op
from setuptools import setup, find_packages
from setup_helper import SetupHelper

classifiers = [
    "Development Status :: 2 - Pre-Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
]

NAME = "CPhasing"
h = SetupHelper(initfile="cphasing/__init__.py", readmefile="README.md")
setup_dir = op.abspath(op.dirname(__file__))
requirements = [x.strip() for x in open(op.join(setup_dir, "requirements.txt"))]



setup(
    name=NAME,
    author=h.author,
    author_email=h.email,
    version=h.version,
    license=h.license,
    classifiers=classifiers,
    packages=find_packages(),
    zip_safe=False,
    url='http://github.com/wangyibin/CPhasing',
    description="Genome scaffold by Hi-C",
    long_description=h.long_description,
    long_description_content_type='test/markdown',
    install_requires=requirements
)