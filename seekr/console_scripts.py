import sys
import argparse
import numpy as np
import pandas as pd

from .kmer_counts import BasicCounter
from .pearson import pearson


KMER_COUNTS_DOC = """
Description
-----------
Generates a kmer count matrix of m rows by n columns,
where m is the number of transcripts in a fasta file and n is 4^kmer.

Examples
--------
The default settings produce a binary, normalized numpy file:
    $ kmer_counts /path/to/rnas.fa -o /path/to/out.npy

To get a human readable csv file, set the nonbinary flag:
    $ kmer_counts /path/to/rnas.fa -o /path/to/out.csv -nb

If you want to add default labels, also set the label flag:
    $ kmer_counts /path/to/rnas.fa -o /path/to/out.csv -nb -lb

You can change also change the size of the kmer you're using, and prevent normalization:
    $ kmer_counts /path/to/rnas.fa -o /path/to/out.npy -k 4 -uc -us

Notes
-----
For more sophisticated options, you cannot use the command line, but need python instead.

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

PEARSON_DOC = """
Description
-----------
Generate a matrix of Pearson similarities from two kmer count files.

Examples
--------
The default settings accept two numpy files and output a third numpy file.

$ pearson /path/to/kc_out.npy /path/to/kc_out.npy -o /path/to/out.npy

The only other options besides the `-o` flag control binary versus .csv input and output. 
If you have a non binary input file (i.e. a .csv file) and also want a non binary output file, you can do:

$ pearson /path/to/kc_out.npy /path/to/kc_out.npy -o /path/to/out.npy -nbi -nbo

Notes
-----
For more sophisticated options, you cannot use the command line, but need python instead.

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

NORM_VECTORS_DOC = """
Description
-----------
Generate two .npy files from a .fa file to use as normalization vectors for other .fa files.

Examples
--------
The default setting accept a single fasta file.

$ norm_vectors gencode.fa

If you want to specify paths for the output files, or choose a different kmer size:

$ norm_vectors gencode.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

def _parse_args_or_exit(parser):
    """"""
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def _run_kmer_counts(fasta, outfile, kmer, nonbinary, centered, standardized,
                     label, mean_vector, std_vector):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    mean = mean_vector or centered
    std = std_vector or standardized
    counter = BasicCounter(fasta, outfile, int(kmer),
                           nonbinary, mean, std, label=label)
    counter.make_count_file()


def console_kmer_counts():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=KMER_COUNTS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='full path of fasta file')
    parser.add_argument('-o', '--outfile', default='counts.seekr',
                        help='name of file to save counts to')
    parser.add_argument('-k', '--kmer', default=6,
                        help='length of kmers you want to count')
    parser.add_argument('-nb', '--nonbinary', action='store_false',
                        help='select if output should be a csv file')
    parser.add_argument('-uc', '--uncentered', action='store_false',
                        help='select if output should not have the mean subtracted')
    parser.add_argument('-us', '--unstandardized', action='store_false',
                        help='select if output should not be divided by the standard deviation')
    parser.add_argument('-lb', '--label', action='store_true',
                        help='select to save with fasta header labels.')
    parser.add_argument('-mv', '--mean_vector', default=None,
                        help='optional path to mean vector numpy file')
    parser.add_argument('-sv', '--std_vector', default=None,
                        help='optional path to std vector numpy file')
    args = _parse_args_or_exit(parser)
    _run_kmer_counts(args.fasta, args.outfile, args.kmer, args.nonbinary, args.uncentered,
                     args.unstandardized, args.label, args.mean_vector, args.std_vector)


def _run_pearson(counts1, counts2, outfile, nonbinary_input, nonbinary_output):
    names1 = None
    names2 = None
    if nonbinary_input:
        counts1 = pd.read_csv(counts1, index_col=0)
        counts2 = pd.read_csv(counts2, index_col=0)
        names1 = counts1.index.values
        names2 = counts2.index.values
    else:
        counts1 = np.load(counts1)
        counts2 = np.load(counts2)

    if nonbinary_output:
        dist = pearson(counts1, counts2)
        dist = pd.DataFrame(dist, names1, names2)
        dist.to_csv(outfile)
    else:
        pearson(counts1, counts2, outfile=outfile)


def console_pearson():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=PEARSON_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('counts1',
                        help='full path of a count file produced by kmer_counts.py')
    parser.add_argument('counts2',
                        help=('full path of a second count file produced by kmer_counts.py. '
                              'This can be the same path as the first counts file.'))
    parser.add_argument('-o', '--outfile', default='pearson.seekr',
                        help='path of file to save similarities to')
    parser.add_argument('-nbi', '--nonbinary_input', action='store_true',
                        help='select if the input will be a csv file')
    parser.add_argument('-nbo', '--nonbinary_output', action='store_true',
                        help='select if output should be a csv file')
    args = _parse_args_or_exit(parser)
    _run_pearson(args.counts1, args.counts2, args.outfile,
                 args.nonbinary_input, args.nonbinary_output)


def _run_norm_vectors(fasta, mean_vector, std_vector, kmer):
    counter = BasicCounter(fasta, k=int(kmer))
    counter.get_counts()
    np.save(mean_vector, counter.mean)
    np.save(std_vector, counter.std)


def console_norm_vectors():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=NORM_VECTORS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='full path of fasta file')
    parser.add_argument('-mv', '--mean_vector', default='mean.npy',
                        help='path to output mean vector')
    parser.add_argument('-sv', '--std_vector', default='std.npy',
                        help='path to output standard deviation vector')
    parser.add_argument('-k', '--kmer', default=6,
                        help='length of kmers you want to count')
    args = _parse_args_or_exit(parser)
    _run_norm_vectors(args.fasta, args.mean_vector, args.std_vector, args.kmer)
