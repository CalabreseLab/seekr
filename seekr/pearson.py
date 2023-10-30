#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description
-----------
Generate a matrix of Pearson similarities from two kmer count files.

Examples
--------
Generate a small, plain text .csv file:

```
import pandas as pd
from seekr.pearson import pearson

dist = pearson(counts1=example_2mers,
               counts2=example_2mers)
pd.DataFrame(dist).to_csv('example_vs_example_2mers.csv')

Notes
-----

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

import numpy as np

def pearson(counts1, counts2, row_standardize=True, outfile=None):
    """Calculates a column standardized Pearson correlation matrix"""
    if row_standardize:
        counts1 = (counts1.T - np.mean(counts1, axis=1)).T
        counts1 = (counts1.T / np.std(counts1, axis=1)).T
        counts2 = (counts2.T - np.mean(counts2, axis=1)).T
        counts2 = (counts2.T / np.std(counts2, axis=1)).T

    # Take the inner product and save
    dist = np.inner(counts1, counts2) / counts1.shape[1]
    if outfile:
        np.save(outfile, dist)
    return dist
