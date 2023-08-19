#!/usr/bin/env python
# -*- coding:utf-8 -*-

import io 
import re
import os
import os.path as op
from setuptools import setup, find_packages

classifiers = [
    "Development Status :: 2 - Pre-Alpha",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
]

NAME = "cphasing"
setup_dir = op.abspath(op.dirname(__file__))
requirements = [x.strip() for x in open(op.join(setup_dir, "requirements.txt"))]

def _read(*args, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *args)
    encoding = kwargs.pop("encoding", "utf-8")
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text

def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("cphasing", "__init__.py"),
        re.MULTILINE,
    ).group(1)

    return version

def get_author():
    author = re.search(
        r'^__author__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("cphasing", "__init__.py"),
        re.MULTILINE,
    ).group(1)

    return author 

def get_email():
    email = re.search(
        r'^__email__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("cphasing", "__init__.py"),
        re.MULTILINE,
    ).group(1)

    return email 

def get_license():
    license = re.search(
        r'^__license__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("cphasing", "__init__.py"),
        re.MULTILINE,
    ).group(1)

    return license


setup(
    name=NAME,
    author=get_author(),
    author_email=get_email(),
    version=get_version(),
    license=get_license(),
    classifiers=classifiers,
    packages=find_packages(),
    scripts=['bin/cphasing', 'bin/cphasing-cli'],
    zip_safe=False,
    url='http://github.com/wangyibin/CPhasing',
    description="Phasing and scaffolding polyploid genomes based on Pore-C or Hi-C data",
    install_requires=requirements
)