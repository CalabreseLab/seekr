import numpy as np
import pandas as pd
import pkg_resources

from seekr import console_scripts


class TestConsoleScripts:

    def test_run_kmer_counts(self, tmpdir):
        infasta = 'tests/data/example.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        outfile = str(tmpdir.join('2mers.npy'))
        console_scripts._run_kmer_counts(fasta=infasta,
                                         outfile=outfile,
                                         kmer='2',
                                         nonbinary=True,
                                         centered=True,
                                         standardized=True,
                                         label=False,
                                         mean_vector=None,
                                         std_vector=None)
        kmers = np.load(outfile)
        expected = 'tests/data/example_2mers.npy'
        expected = pkg_resources.resource_filename('seekr', expected)
        expected = np.load(expected)
        assert np.allclose(kmers, expected)

    def test_run_kmer_counts_raw_csv(self, tmpdir):
        infasta = 'tests/data/example.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        outfile = str(tmpdir.join('3mers.csv'))
        console_scripts._run_kmer_counts(fasta=infasta,
                                         outfile=outfile,
                                         kmer='3',
                                         nonbinary=False,
                                         centered=False,
                                         standardized=False,
                                         label=False,
                                         mean_vector=None,
                                         std_vector=None)
        kmers = pd.read_csv(outfile, header=None)
        expected = 'tests/data/example_3mers_raw.csv'
        expected = pkg_resources.resource_filename('seekr', expected)
        expected = pd.read_csv(expected, header=None)
        assert np.allclose(kmers.values, expected.values)

    def test_run_kmer_counts_vectors(self, tmpdir):
        infasta = 'tests/data/example.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        mean_vector = 'tests/data/example_mean.npy'
        mean_vector = pkg_resources.resource_filename('seekr', mean_vector)
        std_vector = 'tests/data/example_std.npy'
        std_vector = pkg_resources.resource_filename('seekr', std_vector)
        outfile = str(tmpdir.join('2mers_vectors.npy'))
        console_scripts._run_kmer_counts(fasta=infasta,
                                         outfile=outfile,
                                         kmer='2',
                                         nonbinary=True,
                                         centered=False,
                                         standardized=False,
                                         label=False,
                                         mean_vector=mean_vector,
                                         std_vector=std_vector)
        kmers = np.load(outfile)
        expected = 'tests/data/example_2mers.npy'
        expected = pkg_resources.resource_filename('seekr', expected)
        expected = np.load(expected)
        assert np.allclose(kmers, expected)

    def test_run_norm_vectors(self, tmpdir):
        infasta = 'tests/data/example.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        mean = str(tmpdir.join('mean.npy'))
        std = str(tmpdir.join('std.npy'))
        console_scripts._run_norm_vectors(fasta=infasta,
                                          mean_vector=mean,
                                          std_vector=std,
                                          kmer=2)
        mean = np.load(mean)
        std = np.load(std)
        expected_mean = 'tests/data/example_mean.npy'
        expected_mean = pkg_resources.resource_filename('seekr', expected_mean)
        expected_mean = np.load(expected_mean)
        expected_std = 'tests/data/example_std.npy'
        expected_std = pkg_resources.resource_filename('seekr', expected_std)
        expected_std = np.load(expected_std)
        assert np.allclose(mean, expected_mean)
        assert np.allclose(std, expected_std)
