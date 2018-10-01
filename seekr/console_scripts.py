import sys
import argparse
import numpy as np
import pandas as pd

from .kmer_counts import BasicCounter
from .pearson import pearson

def console_kmer_counts():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='full path of fasta file')
    parser.add_argument('-o', '--outfile', default='counts.seekr', help='name of file to save counts to')
    parser.add_argument('-k', '--kmer', default=6, help='length of kmers you want to count')
    parser.add_argument('-nb', '--nonbinary', action='store_false', help='select if output should be a csv file')
    parser.add_argument('-uc', '--uncentered', action='store_false', help='select if output should not have the mean subtracted')
    parser.add_argument('-us', '--unstandardized', action='store_false', help='select if output should not be divided by the standard deviation')
    parser.add_argument('-lb', '--label', action='store_true', help='select to save with fasta header labels.')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    counter = BasicCounter(args.fasta, args.outfile, int(args.kmer),
                           args.nonbinary, args.uncentered,
                           args.unstandardized, label=args.label)
    counter.make_count_file()


def console_pearson():
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
