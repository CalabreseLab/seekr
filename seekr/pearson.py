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
from seekr.kmer_counts import BasicCounter, Log2
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
        cols = np.random.choice(adj.shape[1], size, replace=False)
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
    reference_path: str (default=None)
        Path to third fasta file containing sequences to be used for comparison when calculating
        percentile values of the r-values between the query and targets (e.g. mouse transcriptome).
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
    log2: Log2 (default=Log2.post)
        Log2 transformation can occur pre- or post-standardization, or not at all.
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
    r_values: DataFrame
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

    def __init__(self, query_path=None, target_path=None, reference_path=None, r_values_path=None,
                 percentiles_path=None, mean=True, std=True, log2=Log2.post,
                 k=6, window=1000, slide=100):
        self.query_path = query_path
        self.target_path = target_path
        if reference_path is None and percentiles_path is not None:
            msg = "To calculate percentiles, pass a path to a reference fasta file."
            raise ValueError(msg)
        self.reference_path = reference_path
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
        self.r_values = None
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
            end = i + self.window
            tiles.append(target[i: end])
        tiles[-1] += target[end:]
        counter = BasicCounter(k=self.k, mean=self.mean, std=self.std, log2=self.log2)
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
        return r_vals

    def r_values2df(self, r_values):
        """Convert r_values from list of numpy arrays to a labeled DataFrame."""
        self.r_values = np.hstack(r_values)
        self.r_values = pd.DataFrame(data=self.r_values,
                                     index=self.row_labels,
                                     columns=self.column_labels)

    def percentileofscore(self, a, score):
        """Same as scipy's percentileofscore where kind='rank'.

        See: https://github.com/scipy/scipy/blob/master/scipy/stats/stats.py#L1846

        Notes
        -----
        This is copied here to avoid another large dependency.
        If scipy eventually becomes necessary, remove this code.
        """
        left = np.count_nonzero(a < score)
        right = np.count_nonzero(a <= score)
        pct = (right + left + (1 if right > left else 0)) * 50.0/len(a)
        return pct

    def calc_percentiles(self, query_counts):
        """Find percentiles of elements in r_values, across each row/query, relative to a reference.

        Notes
        -----
        Repeatedly calling `percentileofscore` is probably really slow.
        There should be some way of speeding things up, similar to `calc_internal_percentiles`.
        I don't have enough time to look into it now.

        Returns
        -------
        percentiles_df: DataFame
            Percentile equivalent of each element in r_values, across each row.
        """
        counter = BasicCounter(infasta=self.reference_path, k=self.k, mean=self.mean, std=self.std,
                               log2=self.log2)
        counter.get_counts()
        query_vs_ref_rvals = pearson(query_counts, counter.counts)[0]
        percentiles_df = pd.DataFrame(columns=self.r_values.columns)
        for index, row in self.r_values.iterrows():
            percentiles = []
            for query_vs_tile_rval in row:
                percentiles.append(self.percentileofscore(query_vs_ref_rvals,
                                                          query_vs_tile_rval))
            percentiles_df.loc[index] = percentiles
        return percentiles_df

    def calc_internal_percentiles(self):
        """Calculate percentile scores of each element in r_values, across each row/query.

        Notes
        -----
        This method doesn't need a reference to calculate percentiles.
        It currently isn't called durning DomainPearson.run()

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
        r_values = []
        for target_header, target in my_tqdm()(target_headers_seqs):
            r_values.append(self.compare_query_target(query_counts, target_header, target))
        self.r_values2df(r_values)
        self.percentiles = self.calc_percentiles(query_counts)
        self.save()
