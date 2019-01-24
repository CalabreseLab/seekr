# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:16:50 2016

@author: jessime

A couple of functions for autodetecting if code is being run in the notebook.
The appropriate tqdm progressbar will be returned.
"""

from tqdm import tnrange, trange, tqdm, tqdm_notebook


import sys


def _is_kernel():
    if 'IPython' not in sys.modules:
        # IPython hasn't been imported, definitely not
        return False
    from IPython import get_ipython
    # check for `kernel` attribute on the IPython instance
    return getattr(get_ipython(), 'kernel', None) is not None


def my_tqdm():
    return tqdm_notebook if _is_kernel() else tqdm


def my_trange():
    return tnrange if _is_kernel() else trange
