#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description
-----------
Generate a matrix of Pearson similarites from two kmer count files.

Examples
--------
The default setting accept two numpy files and output a third numpy file.

$ python pearson.py /path/to/kc_out.npy /path/to/kc_out.npy -o /path/to/out.npy

The only other options besides the `-o` flag control binary versus .csv input and output. If you have a non binary input file (i.e. a .csv file) and also want a non binary output file, you can do:

$ python pearson.py /path/to/kc_out.npy /path/to/kc_out.npy -o /path/to/out.npy -nbi -nbo

Notes
-----
For more sophisticated options, you cannot use the command-line, but need python instead.

Issues
------
Any issues can be reported to https://github.com/CalabreseLab #TODO

---
"""

import sys
import argparse
import numpy as np
import pandas as pd

from pickle import dump

try:
    from scipy.stats import t
except ImportError:
    pass
from kmer_counts import BasicCounter
from tqdm import tqdm, tqdm_notebook

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


def cmd_line_pearson():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('counts1', help='full path of first count file produced by kmer_counts.py')
    parser.add_argument('counts2', help='full path of second count file produced by kmer_counts.py. This can be the same path as the first counts file.')
    parser.add_argument('-o', '--outfile', default='pearson.seekr', help='path of file to save similarities to')
    parser.add_argument('-nbi', '--nonbinary_input', action='store_true', help='select if the input will be a csv file')
    parser.add_argument('-nbo', '--nonbinary_output', action='store_true', help='select if output should be a csv file')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    names1 = None
    names2 = None
    if args.nonbinary_input:
        counts1 = pd.read_csv(args.counts1, index_col=0)
        counts2 = pd.read_csv(args.counts2, index_col=0)
        names1 = counts1.index.values
        names2 = counts2.index.values
    else:
        counts1 = np.load(args.counts1)
        counts2 = np.load(args.counts2)

    if args.nonbinary_output:
        dist = pearson(counts1, counts2)
        dist = pd.DataFrame(dist, names1, names2)
        dist.to_csv(args.outfile)
    else:
        pearson(counts1, counts2, outfile=args.outfile)


if __name__ == '__main__':
    cmd_line_pearson()
