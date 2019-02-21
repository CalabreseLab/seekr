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

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pickle import dump

try:
    from scipy.stats import t
except ImportError:
    pass

from seekr.my_tqdm import my_tqdm
from seekr.fasta_reader import Reader
from seekr.kmer_counts import BasicCounter
from seekr.utils import get_adj


def pearson(counts1, counts2, row_standardize=True, outfile=None):
    """Calculates a column standardized Pearson correlation matrix"""
    if row_standardize:
        counts1 = (counts1.T - np.mean(counts1, axis=1)).T
        counts1 = (counts1.T / np.std(counts1, axis=1)).T
        counts2 = (counts2.T - np.mean(counts2, axis=1)).T
        counts2 = (counts2.T / np.std(counts2, axis=1)).T

    # Take the inner product and save
    dist = np.inner(counts1, counts2)/counts1.shape[1]
    if outfile:
        np.save(outfile, dist)
    return dist


def visualize_distro(adj, out_path, sample=None):
    adj = get_adj(adj)
    if isinstance(adj, pd.DataFrame):
        adj = adj.values
    if sample is not None:
        if not 0 < sample <= 1:
            raise ValueError('Value of sample must satisfy: 0 < sample <= 1')
        size = int(len(adj) * sample)
        rows = np.random.choice(len(adj), size, replace=False)
        cols = np.random.choice(len(adj), size, replace=False)
        adj = adj[rows][:, cols]
    flat = adj.ravel()
    mean = flat.mean()
    std = flat.std()
    std1 = mean + std
    std2 = mean + (2*std)
    ax = sns.distplot(flat, label=f'Distro (n={flat.size})')
    y_max = ax.get_ylim()[1]
    plt.plot((mean, mean), (0, y_max/2), label='Mean')
    plt.plot((std1, std1), (0, y_max/2), label='Mean + 1 std. dev.')
    plt.plot((std2, std2), (0, y_max/2), label='Mean + 2 std. dev.')
    plt.xlabel('r-value')
    plt.legend()
    plt.savefig(out_path, bbox_inches='tight', dpi=600)
    return mean, std


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


class DomainPearson:
    """Calculate r-values and percentiles over windows of target sequences for a list of queries.

    Parameters
    ----------
    query_path: str (default=None)
        Path to fasta file containing transcripts of interest (e.g. Xist-2kb).
    target_path: str (default=None)
        Path to second fasta file which will be tiled to find domains similar to query transcripts.
    r_values_path: str (default=None)
        Path to csv file containing pairwise comparisons between queries and target tiles.
    percentiles_path: str (default=None)
        Path to csv file containing the percentile equivalent of each element in r_value_path's csv.
    mean: bool, np.array, str (default=True)
        Set the mean to 0 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated mean array.
        Can be produced by `seekr_norm_vectors`.
    std: bool, np.array, str (default=True)
        Set the std. dev. to 1 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated std array.
        Can be produced by `seekr_norm_vectors`.
    log2: bool (default=True)
        If False, do not apply a log2 transform to the count matrix.
    k: int (default=6)
        Size of kmer to be counted.
    window: int (default=1000)
        Size of tile/domain to be created from target transcripts for comparison against queries.
    slide: int (default=100)
        Number of basepairs to move along target transcript before creating another tile/domain.

    Attributes
    ----------
    row_labels: List[str]
        Names of query transcripts from fasta file.
    column_labels List[str]
        Names of tiles as created at "{target fasta header}_{tile number}".
    r_values: List[ndarray] | DataFrame
        n by m array of Pearson r-values.
        n is the number of query transcripts and m is the total number of target tiles/domains.
    percentiles: DataFrame
        Percentile equivalent of each element in r_values, across each row.

    Notes
    -----
    This code is inefficient.
    Email me (jessime.kirk@gmail.com) if you're running into memory/speed issues.
    There are 'clever' ways of doing this that are a lot messier and harder to maintain.
    """

    def __init__(self, query_path=None, target_path=None, r_values_path=None, percentiles_path=None,
                 mean=True, std=True, log2=True, k=6, window=1000, slide=100):
        self.query_path = query_path
        self.target_path = target_path
        self.r_values_path = r_values_path
        self.percentiles_path = percentiles_path
        self.mean = mean
        self.std = std
        self.log2 = log2
        self.k = k
        self.window = window
        self.slide = slide

        self.row_labels = []
        self.column_labels = []
        self.r_values = []
        self.percentiles = None

    def make_query_counts(self):
        """Create normalized counts for all query transcripts.

        Returns
        -------
        query_counts: ndarray
            Normalized kmer profiles for each query transcript.
        """
        counter = BasicCounter(self.query_path, k=self.k, mean=self.mean, std=self.std,
                               log2=self.log2)
        counter.get_counts()
        return counter.counts

    def get_target_tile_counts(self, target):
        """Create normalized counts for all tiles in a single target transcript.

        If a target sequence is smaller than window size, the whole sequence is used for comparison.

        Parameters
        ----------
        target: str
            Sequence from which to create tiles and kmer counts.

        Returns
        -------
        target_counts: ndarray
            Normalized kmer profiles for each tile in target transcript.
        """
        tiles = []
        upper = max(1, len(target) - self.window + 1)
        for i in range(0, upper, self.slide):
            tiles.append(target[i:i + self.window])
        counter = BasicCounter(k=self.k, mean=self.mean, std=self.std, log2=self.log2, silent=True)
        counter.seqs = tiles
        counter.get_counts()
        return counter.counts

    def update_column_labels(self, target_header, target_counts):
        """Add new column labels from all tiles in latest target transcript.

        Parameters
        ----------
        target_header: str
            Name of target from fasta file.
        target_counts: ndarray
            Normalized kmer profiles for each tile in target transcript.
        """
        self.column_labels.extend(target_header + f'_{i}' for i in range(len(target_counts)))

    def compare_query_target(self, query_counts, target_header, target):
        """Calculare r-values for """
        target_counts = self.get_target_tile_counts(target)
        self.update_column_labels(target_header, target_counts)
        r_vals = pearson(query_counts, target_counts)
        self.r_values.append(r_vals)

    def r_values2df(self):
        """Convert r_values from list of numpy arrays to a labeled DataFrame."""
        self.r_values = np.hstack(self.r_values)
        self.r_values = pd.DataFrame(data=self.r_values,
                                     index=self.row_labels,
                                     columns=self.column_labels)

    def calc_percentiles(self):
        """Calculate percentile scores of each element in r_values, across each row/query.

        Returns
        -------
        percentiles: DataFame
            Percentile equivalent of each element in r_values, across each row.
        """
        percentiles = self.r_values.rank(axis=1, pct=True) * 100
        return percentiles

    def save(self):
        """Save csv files of r_values and percentiles, using limited precision float values."""
        if self.r_values_path is not None:
            self.r_values.to_csv(self.r_values_path, float_format='%.4f')
        if self.percentiles_path is not None:
            self.percentiles.to_csv(self.percentiles_path, float_format='%.4f')

    def run(self):
        self.row_labels = Reader(self.query_path).get_headers()
        query_counts = self.make_query_counts()
        target_reader = Reader(self.target_path)
        target_headers_seqs = target_reader.get_data(tuples_only=True)
        for target_header, target in my_tqdm()(target_headers_seqs):
            self.compare_query_target(query_counts, target_header, target)
        self.r_values2df()
        self.percentiles = self.calc_percentiles()
        self.save()


class StreamDist:
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
