#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Custom exception of cphasing
"""


import argparse
import logging
import os
import os.path as op
import sys



class IsNotPairs(Exception):
    """
    Exception for not pairs.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
    