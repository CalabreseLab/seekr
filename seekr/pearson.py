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
import pandas as pd

from pickle import dump

try:
    from scipy.stats import t
except ImportError:
    pass
from tqdm import tqdm, tqdm_notebook

from .kmer_counts import BasicCounter

def pearson(counts1, counts2, row_standardize=True, outfile=None):
    """Calculates a column standardized Pearson correlation matrix"""
    if row_standardize:
        counts1 = (counts1.T - np.mean(counts1, axis=1)).T
        counts1 = (counts1.T / np.std(counts1, axis=1)).T
        counts2 = (counts2.T - np.mean(counts2, axis=1)).T
        counts2 = (counts2.T / np.std(counts2, axis=1)).T

    #Take the inner product and save
    dist = np.inner(counts1, counts2)/counts1.shape[1]
    if outfile:
        np.save(outfile, dist)
    return dist

def pvalue(dist, N):
    """Calculate the p-values (two tailed) of a Pearson correlation matrix

    Parameters
    ----------
    dist : ndarray
        Pearson correlation matrix
    N : int
        Sample size. (The length of each array used to find the R value)

    Returns
    -------
    pvals : ndarray
        p-value for each R value

    Reference
    ---------
    http://stats.stackexchange.com/questions/120199/calculate-p-value-for-the-correlation-coefficient
    http://stackoverflow.com/questions/17559897/python-p-value-from-t-statistic
    """
    pvals = dist / np.sqrt((1-dist**2) / (N -2))
    pvals = t.sf(np.abs(pvals), N - 2)*2
    return pvals

class StreamDist(object):
    """Pearson's R from fasta files too large for memory."""

    def __init__(self, infasta=None, query=None, outfile=None,
                 mean=None, std_dev=None, k=6, n=0,
                 norm=True, binary=True, nb=False):
        self.infasta = infasta
        self.query = query
        self.outfile = outfile
        self.mean = mean
        self.std_dev = std_dev
        self.k = k
        self.norm = norm
        self.binary = binary
        self.nb = nb

        self.counter = BasicCounter(k=self.k, silent=True)
        self.query_counts = None
        self.dist = None
        self.names = None
        self.n = n
        self.i = 0

        #Annoying, but necessary for init of self.dist
        if self.mean is not None or self.std_dev is not None:
            assert self.n > 0, 'Please provide number of sequences in fasta.'

    def get_names(self):
        """Stream over the file again and pull out header lines."""
        names = []
        with open(self.infasta) as infile:
            bar = self._progress()
            for line in bar(infile):
                if line[0] == '>':
                    names.append(line.strip())
        self.names = names

    def _progress(self):
        if self.nb:
            return tqdm_notebook
        else:
            return tqdm

    def _single_count(self, line):
        seq = line.strip().upper()
        row = np.zeros(4**self.k, dtype=np.float32)
        self.counter.seqs = seq
        row = self.counter.occurrences(row, seq)
        return row

    def _stream_seqs(self, func):
        """Perform a function on the counts of each fasta sequence"""
        with open(self.infasta) as infasta:
            bar = self._progress()
            if self.n == 0:
                total = None
            else:
                total = self.n * 2
            for line in bar(infasta, total=total):
                if line[0] != '>':
                    counts = self._single_count(line)
                    func(counts)

    def make_query_counts(self):
        with open(self.query) as query:
            seq = query.readlines()[1].strip().upper()
        self.query_counts = self._single_count(seq)
        if self.norm:
            self.query_counts = self.col_norm(self.query_counts)
        self.query_counts = self.row_norm(self.query_counts)

    def online_moments(self, counts):
        """Update the mean and std. vectors"""
        self.n += 1
        delta = counts - self.mean
        self.mean = self.mean + delta/self.n
        self.std_dev = self.std_dev + delta*(counts - self.mean)

    def calc_norm_vectors(self):
        self.mean = np.zeros(4**self.k, dtype=np.float32)
        self.std_dev = np.zeros(4**self.k, dtype=np.float32)
        self._stream_seqs(self.online_moments)
        self.std_dev = self.std_dev /self.n
        self.std_dev = np.sqrt(self.std_dev)

    def norm_vectors(self):
        """Either load or count norm vectors if they are not passed as arrays."""
        if self.mean is None and self.std_dev is None:
            self.calc_norm_vectors()
        elif isinstance(self.mean, str) and isinstance(self.std_dev, str):
            self.mean = np.load(self.mean).astype(dtype=np.float32)
            self.std_dev = np.load(self.std_dev).astype(dtype=np.float32)
        else:
            pass

    def col_norm(self, counts):
        """Column normalize"""
        counts = counts - self.mean
        counts = counts / self.std_dev
        return counts

    def row_norm(self, counts):
        """Row normalize"""
        counts = (counts.T - np.mean(counts)).T
        counts = (counts.T / np.std(counts)).T
        return counts

    def calc_dist(self, counts):
        if self.norm:
            counts = self.col_norm(counts)
        counts = self.row_norm(counts)
        #self.dist[self.i] = np.inner(self.query_counts, counts)/len(counts)
        self.dist[self.i] = np.sum(self.query_counts * counts)/len(counts)
        self.i += 1

    def save(self):
        """Saves the counts appropriately based on current settings"""
        if self.outfile is not None:
            if not self.binary:
                np.savetxt(self.outfile, self.dist, delimiter=',')
            else:
                try:
                    np.save(self.outfile, self.dist)
                except AttributeError:
                    outfile = open(self.outfile, 'w+b')
                    dump(self.dist.tolist(), outfile)
                    outfile.close()

    def make_dist(self):
        """Main"""
        if self.norm:
            self.norm_vectors()
        self.make_query_counts()
        self.dist = np.zeros(self.n, dtype=np.float32)
        self._stream_seqs(self.calc_dist)
        self.save()


if __name__ == '__main__':
    cmd_line_pearson()
