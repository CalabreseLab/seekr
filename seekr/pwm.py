# TODO (Dan) What should this file be called?
import pandas as pd
import numpy as np

from collections import defaultdict
from itertools import product
from pathlib import Path


class CountsWeighter:
    """Weight kmer counts by a collection of PWMs.

    Parameters
    ----------
    pwm_dir: str (default=None)
        Path to directory containing pwm files patterned *.txt.
    counts: str | ndarray | DataFrame (default=None)
        (Path to) kmer counts matrix.
    k: int
        Length of kmer.
    out_path: str (default=None)
        Path to .csv file for weighted counts

    Attributes
    ----------
    k_sub: int (default=4)
        Length of sub_kmers to use if motif length is less than k.
        Currently hard-coded to 4.
    kmers: list
        Str elements of all kmers of size k
    df: pd.DataFrame (default=None)
        # TODO (Dan) One line description of contents
    """
    def __init__(self, pwm_dir=None, counts=None, k=5, out_path=None):
        self.pwm_dir = pwm_dir
        if pwm_dir is not None:
            self.pwm_dir = Path(pwm_dir)
        self.counts = counts
        self.k = k
        self.out_path = out_path

        self.k_sub = 4
        self.kmers = [''.join(p) for p in product('AGTC', repeat=self.k)]
        self.df = None
        if counts is not None:  # get_counts depends on self.kmers
            self.counts = self.get_counts(counts)

    def get_counts(self, counts):
        """Load kmer counts matrix from .csv or .npy file, if necessary.

        Parameters
        ----------
        counts: str | ndarray | DataFrame
            (Path to) kmer counts matrix.

        Returns
        -------
        counts: DataFrame
            Kmer counts matrix describing transcript kmer profiles.
        """
        counts_types = (str, pd.DataFrame, np.ndarray)
        err_msg = f'adj must be one of {counts_types}, not {type(counts)}.'
        assert type(counts) in counts_types, err_msg
        if isinstance(counts, str):
            try:
                counts = pd.read_csv(counts, index_col=0)
            except UnicodeDecodeError:
                counts = np.load(counts)
        if isinstance(counts, np.ndarray):
            counts = pd.DataFrame(counts, columns=self.kmers)
        return counts

    def gen_pwm_dicts(self):
        """Search directory for non-empty pwm files, and load them.

        Yields
        ------
        pwm_path: str
            Path to non-empty pwm file.
        pwm: Dict[str: Dict[int: float]]
            Position weight matrix as a nested dict (for fast lookups).
            Outer keys are nucleotides. Inner keys are positions. Inner values are weights.
        """
        for pwm_path in self.pwm_dir.glob('*.txt'):
            try:
                pwm = pd.read_csv(pwm_path, sep='\t')
            except pd.errors.EmptyDataError:
                print(f'The motif file {pwm_path} is empty. Skipping.')
                continue
            pwm.drop('Pos', axis=1, inplace=True)
            pwm = pwm.rename(columns={'U': 'T'}).to_dict()
            yield pwm_path, pwm

    def set_kmer2weight(self, kmer2weight, pwm, key_kmer, iter_kmer, n_kmers):
        for frame in range(n_kmers):
            weight = 1
            for pos, nucleotide in enumerate(iter_kmer):
                weight *= pwm[nucleotide][pos + frame]
            kmer2weight[key_kmer] += weight

    def build_weights_dict(self, pwm):
        """Create a mapping between all kmers their weights in a given PWM.

        Parameters
        ----------
        pwm: Dict[str: Dict[int: float]]
            Position weight matrix.

        Returns
        -------
        kmer2weight: Dict
            Keys are all kmers, values are the kmer's weight in a given pwm.
        """
        kmer2weight = defaultdict(int)
        motif_len = len(pwm['A'])
        if motif_len < self.k:
            n_kmers = motif_len - 4 + 1
            for kmer in self.kmers:
                for sub_kmer in (kmer[i:i+self.k_sub] for i in range(self.k-self.k_sub+1)):
                    self.set_kmer2weight(kmer2weight, pwm, kmer, sub_kmer, n_kmers)
        else:
            for kmer in self.kmers:
                n_kmers = motif_len - self.k + 1
                self.set_kmer2weight(kmer2weight, pwm, kmer, kmer, n_kmers)
        return kmer2weight

    def weight_counts(self, kmer2weight):
        """Weight each kmer in counts matix by corresponding kmer weight found in a PWM.

        Parameters
        ----------
        kmer2weight: Dict
            Keys are all kmers, values are the kmer's weight in a given pwm.

        Returns
        -------
        scores_sums: np.array
            Each element corresponds to sum of weighted kmer profile for given transcript.
            Length of scores_sums is equal to length of counts.
        """
        sorted_weights = np.array([kmer2weight[k] for k in self.counts.columns])
        weighted_z_scores = self.counts.values.copy() * sorted_weights
        scores_sums = weighted_z_scores.sum(axis=1)
        return scores_sums

    def save(self):
        """Save df to csv file."""
        if self.out_path is not None:
            self.df.to_csv(self.out_path, float_format='%.4f')

    def run(self):
        """TODO (Dan) Update based on description of self.df"""
        score_dict = {}
        for pwm_path, pwm in self.gen_pwm_dicts():
            kmer2weight = self.build_weights_dict(pwm)
            score_dict[pwm_path.name] = self.weight_counts(kmer2weight)
        self.df = pd.DataFrame.from_dict(score_dict, orient='index', columns=self.counts.index)
        self.save()
