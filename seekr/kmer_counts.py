#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description
-----------
Generates a kmer count matrix of m rows by n columns,
where m is the number of transcripts in a fasta file and n is 4^kmer.

Examples
--------
Generate a small, plain text .csv file:

```
from seekr.kmer_counts import BasicCounter

counter = BasicCounter(infasta='example.fa',
                       outfile='example_2mers.csv',
                       k=2,
                       binary=False)
counts = counter.make_count_file()
```

Notes
-----


Issues
------
Any issues can be reported to https://github.com/CalabreseLab

---
"""

import numpy as np

from enum import Enum
from collections import defaultdict
from itertools import product
from pandas import DataFrame

from seekr.my_tqdm import my_tqdm
from seekr.fasta_reader import Reader


Log2 = Enum('Log2', ('pre', 'post', 'none'))


class BasicCounter:
    """Generates overlapping kmer counts for a fasta file

    Parameters
    ----------
    infasta: str (default=None)
        Full path to fasta file to be counted
    outfile: str (default=None)
        Full path to the counts file to be saved
    k: int (default=6)
        Size of kmer to be counted
    binary: bool (default=True)
        Saves as numpy array if True, else saves as csv
    mean: bool, np.array, str (default=True)
        Set the mean to 0 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated mean array.
    std: bool, np.array, str (default=True)
        Set the std. dev. to 1 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated std array.
    log2: Log2 (default=Log2.post)
        Log2 transformation can occur pre- or post-standardization, or not at all.
    leave: bool (default=True)
        Set to False if get_counts is used within another tqdm loop
    silent: bool (default=False)
        Set to True to turn off tqdm progress bar
    alphabet: str (default='AGTC')
        Valid letters to include in kmer.

    Attributes
    ----------
    counts: np.array
        Matrix of kmer counts. Dimensions are equal to # of transcripts by # of kmers.
    kmers: list
        Str elements of all kmers of size k
    map: dict
        Mapping of kmers to column values
    alpha_len: int
        Length of alphabet
    """

    def __init__(self, infasta=None, outfile=None, k=6,
                 binary=True, mean=True, std=True, log2=Log2.post,
                 leave=True, silent=False, label=False, alphabet='AGTC'):
        self.infasta = infasta
        self.seqs = None
        if infasta is not None:
            self.seqs = Reader(infasta).get_seqs()
        self.outfile = outfile
        self.k = k
        self.binary = binary
        self.mean = mean
        if isinstance(mean, str):
            self.mean = np.load(mean)
        self.std = std
        if isinstance(std, str):
            self.std = np.load(std)
        self.log2 = log2 
        self.leave = leave
        self.silent = silent
        self.label = label
        self.counts = None
        self.alpha_len = len(alphabet)
        self.kmers = [''.join(i) for i in product(alphabet, repeat=k)]
        self.map = {k:i for k,i in zip(self.kmers, range(self.alpha_len**k))}

        if self.seqs is not None:
            if len(self.seqs) == 1 and self.std is True:
                err = ('You cannot standardize a single sequence. '
                       'Please pass the path to an std. dev. array, '
                       'or use raw counts by setting std=False.')
                raise ValueError(err)

        if not isinstance(self.log2, Log2):
            raise TypeError(f'log2 must be one of {list(Log2)}')

    def occurrences(self, row, seq):
        """Counts kmers on a per kilobase scale"""
        counts = defaultdict(int)
        length = len(seq)
        increment = 1000/(length-self.k+1) #fix, original 1000/length, incorrect divisor, should be 1000/(l-k+1)
        for c in range(length-self.k+1):
            kmer = seq[c:c+self.k]
            counts[kmer] += increment
        for kmer, n in counts.items():
            if kmer in self.map:
                row[self.map[kmer]] = n
        return row

    def _progress(self):
        """Determine which iterator to loop over for counting"""
        if self.silent:
            return self.seqs

        if not self.leave:
            tqdm_seqs = my_tqdm()(self.seqs, desc='Kmers', leave=False)
        else:
            tqdm_seqs = my_tqdm()(self.seqs)

        return tqdm_seqs

    def center(self):
        """Mean center counts by column"""
        if self.mean is True:
            self.mean  = np.mean(self.counts, axis=0)
        self.counts -= self.mean

    def standardize(self):
        """Divide out the standard deviations from columns of the count matrix"""
        if self.std is True:
            self.std = np.std(self.counts, axis=0)
        self.counts /= self.std
        if np.isnan(self.counts).any():
            print(('\nWARNING: You have `np.nan` values in your counts '
                   'after standardization. This is likely due to '
                   'a kmer not appearing in any of your sequences. '
                   'Try: \n1) using a smaller kmer size, \n2) beginning '
                   'with a larger set of sequences, \n3) passing '
                   'precomputed normalization vectors from a larger '
                   'data set (e.g. GENCODE).'))

    def log2_norm(self):
        """Apply a log2 transform to the count matrix"""
        self.counts += 1
        self.counts = np.log2(self.counts)

    def get_counts(self):
        """Generates kmer counts for a fasta file"""
        self.counts = np.zeros([len(self.seqs), self.alpha_len**self.k], dtype=np.float32)
        seqs = self._progress()

        for i, seq in enumerate(seqs):
            self.counts[i] = self.occurrences(self.counts[i], seq)
        if self.log2 == Log2.pre:
            self.log2_norm()
        if self.mean is not False:
            self.center()
        if self.std is not False:
            self.standardize()
        if self.log2 == Log2.post:
            self.counts += np.abs(np.min(self.counts))
            self.log2_norm()

    def save(self, names=None):
        """Saves the counts appropriately based on current settings.

        There are four output methods for the counts:
        1. Binary. This saves just the counts as a binary numpy array.
        2. No labels. Saves in plain text, but without any labels.
        3. Default names. If no names are provided, fasta headers will be used as labels.
        4. Custom names. Provide a list of names if you want to label lncRNAs with your own names.

        Parameters
        ----------
        names : [str] (default=None)
            Unique names for rows of the Dataframe.
        """
        err_msg = ('You cannot label a binary file. '
                   'Set only one of "binary" or "label" as True. '
                   'If you used `-b` from the command line, '
                   'try also using `-rl`.')
        assert not (self.binary and self.label), err_msg
        assert self.outfile is not None, 'Please provide an outfile location.'
        if self.binary:
            np.save(self.outfile, self.counts)
        elif self.label:
            if names is None:
                names = Reader(self.infasta).get_headers()
            df = DataFrame(data=self.counts, index=names, columns=self.kmers)
            df.to_csv(self.outfile)
        else:
            np.savetxt(self.outfile, self.counts, delimiter=',', fmt='%1.6f')

    def make_count_file(self, names=None):
        """Wrapper function for the most common way to generate count files.

        Given a numpy file name, it will save a numpy file where counts have been:
        cast as a dense array, centered, and standardized.

        Parameters
        ----------
        names : [str] (default=None)
            lncRNA names to pass to self.save

        Returns
        -------
        counts: np.array
            Matrix of kmer counts. Dimensions are equal to # of transcripts by # of kmers.
        """
        self.get_counts()
        if self.outfile is not None:
            self.save(names)
        return self.counts