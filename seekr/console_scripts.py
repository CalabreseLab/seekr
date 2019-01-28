import sys
import argparse
import numpy as np
import pandas as pd

from seekr import fasta
from seekr import graph
from seekr.kmer_counts import BasicCounter
from seekr.pearson import pearson


DOWNLOAD_GENCODE_DOC = """
Description
-----------
Download fasta files from https://www.gencodegenes.org/

The one parameter that must be passed is 'biotype'.
Its value must be one of:
* 'all' : Nucleotide sequences of all transcripts on the reference chromosomes
* 'pc' : Nucleotide sequences of coding transcripts on the reference chromosomes
* 'lncRNA' : Nucleotide sequences of long non-coding RNA transcripts on the reference chromosomes

Examples
--------
To download all human transcripts of the latest release into a fasta file:
    $ seekr_download_gencode all
    
To do the same for mouse:
    $ seekr_download_gencode all -s mouse
    
To get lncRNAs from the M5 release of mouse:
    $ seekr_download_gencode lncRNA -s mouse -r M5
    
If you want to leave the fasta file gzipped:
    $ seekr_download_gencode all -z

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

CANONICAL_GENCODE_DOC = """
Description
-----------
Filter GENCODE fasta file for only transcripts ending in 01.

This is based on the common names provided by GENCODE.
No strict guarantees are made about the relationship between genes and transcripts.

Examples
--------
To filter transcripts ending in 01, an input and output fasta file are required:
    $ seekr_canonical rnas.fa rnas01.fa

If you want to specifically find transcripts with the ending 001:
    $ seekr_canonical rnas.fa rnas01.fa -z 2

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

KMER_COUNTS_DOC = """
Description
-----------
Generates a kmer count matrix of m rows by n columns,
where m is the number of transcripts in a fasta file and n is 4^kmer.

Examples
--------
The default settings take a .fa file and produce a csv file:
    $ seekr_kmer_counts rnas.fa -o out.csv

It's usually helpful to label the csv file with transcript names:
    $ seekr_kmer_counts rnas.fa -o out.csv -lb

To get a compact and efficient .npy file, set the binary flag:
    $ seekr_kmer_counts rnas.fa -o out.npy -b

You can change also change the size of the kmer you're using, and prevent normalization:
    $ seekr_kmer_counts rnas.fa -o out.csv -k 4 -uc -us -nl

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
The default settings accept two csv files and output a third csv file.
    $ seekr_pearson kc_out.csv kc_out.csv -o out.csv

The only other options besides the `-o` flag control binary versus plain text input and output. 
If you have a binary input file (i.e. a .npy file) and also want a binary output file, you can do:
    $ seekr_pearson kc_out.npy kc_out.npy -o out.npy -bi -bo

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
    $ seekr_norm_vectors gencode.fa

If you want to specify paths for the output files, or choose a different kmer size:
    $ seekr_norm_vectors gencode.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

GEN_RAND_RNAS_DOC = """
Description
-----------
Given a .fa file, create a second .fa file with random RNAs based on original RNAs.

Users control how similar the synthetic RNAs are to the originals.
Gross scale of kmer content can be controlled by setting kmer conservation size.
Fine scale control of similarity can be set by number of SNP mutations.

"""


GRAPH_DOC = """
Description
-----------
Find communities of transcripts from an adjacency matrix.

Examples
--------
The default setting accept a single csv file.
This csv file should be the product of seekr_pearson, or some other adjacency matrix.
A gml file contain the graph and communities will be produced.
    $ seekr_graph adj.csv -g graph.gml

For a cleaner csv file of just community information:
    $ seekr_graph adj.csv -g graph.gml -c communities.csv

To threshold edges at a value besides 0, and change the resolution parameter (gamma):
    $ seekr_graph adj.csv -g graph.gml -t .1 -r 1.5
    
To change the cap of the number of communities found, and set the seed:
    $ seekr_graph adj.csv -g graph.gml -n 10 -s 0
    
Numpy files are also valid input:
    $ seekr_graph adj.npy -g graph.gml

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


def _run_download_gencode(biotype, species, release, out_path, unzip):
    # Note: This function is separated for testing purposes.
    downloader = fasta.Downloader()
    downloader.get_gencode(biotype, species, release, out_path, unzip)


def console_download_gencode():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=KMER_COUNTS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('biotype',
                        help=("Name of Genocde set to download. "
                              "Must be one of ('all', 'pc', 'lncRNA')."))
    parser.add_argument('-s', '--species', default='human',
                        help=" Name of species. Must be one of: ('human' or 'mouse').")
    parser.add_argument('-r', '--release', default=None,
                        help=("Name of specific release to download (e.g. 'M5'). "
                              "If None, download latest release."))
    parser.add_argument('-o', '--out_path', default=None,
                        help="Path to location for fasta file. Default will save by release name.")
    parser.add_argument('-z', '--zip', action='store_false',
                        help="Set if you do not want to gunzip fasta file after downloading.")
    args = _parse_args_or_exit(parser)
    _run_download_gencode(args.biotype, args.species, args.release, args.out_path, args.zip)


def _run_canonical_gencode(in_fasta, out_fasta, zeros):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    maker = fasta.Maker(in_fasta, out_fasta)
    maker.filter1(zeros)


def console_canonical_gencode():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=KMER_COUNTS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_fasta', help='Full path of fasta file.')
    parser.add_argument('out_fasta', help='Full path of filtered fasta file.')
    parser.add_argument('-z', '--zeros', help='Number of zeroes needed to be considered canonical.')
    args = _parse_args_or_exit(parser)
    _run_canonical_gencode(args.in_fasta, args.out_fasta, args.zeros)


def _run_kmer_counts(fasta, outfile, kmer, binary, centered, standardized,
                     log2, label, mean_vector, std_vector):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    mean = mean_vector or centered
    std = std_vector or standardized
    counter = BasicCounter(fasta, outfile, int(kmer), binary,
                           mean, std, log2, label=label)
    counter.make_count_file()


def console_kmer_counts():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=KMER_COUNTS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='Full path of fasta file.')
    parser.add_argument('-o', '--outfile', default='counts.seekr',
                        help='Name of file to save counts to.')
    parser.add_argument('-k', '--kmer', default=6,
                        help='Length of kmers you want to count.')
    parser.add_argument('-b', '--binary', action='store_true',
                        help='Select if output should be a .npy file.')
    parser.add_argument('-uc', '--uncentered', action='store_false',
                        help='Select if output should not have the mean subtracted.')
    parser.add_argument('-us', '--unstandardized', action='store_false',
                        help='Select if output should not be divided by the standard deviation.')
    parser.add_argument('-nl', '--no_log2', action='store_false',
                        help='Select if output should not be log2 transformed.')
    parser.add_argument('-lb', '--label', action='store_true',
                        help='Select to save with fasta header labels.')
    parser.add_argument('-mv', '--mean_vector', default=None,
                        help='Optional path to mean vector numpy file.')
    parser.add_argument('-sv', '--std_vector', default=None,
                        help='Optional path to std vector numpy file.')
    args = _parse_args_or_exit(parser)
    _run_kmer_counts(args.fasta, args.outfile, args.kmer, args.binary, args.uncentered,
                     args.unstandardized, args.log2, args.label, args.mean_vector, args.std_vector)


def _run_pearson(counts1, counts2, outfile, binary_input, binary_output):
    # Note: This function is separated for testing purposes.
    names1 = None
    names2 = None
    if binary_input:
        counts1 = np.load(counts1)
        counts2 = np.load(counts2)
    else:
        counts1 = pd.read_csv(counts1, index_col=0)
        counts2 = pd.read_csv(counts2, index_col=0)
        names1 = counts1.index.values
        names2 = counts2.index.values

    if binary_output:
        pearson(counts1, counts2, outfile=outfile)
    else:
        dist = pearson(counts1, counts2)
        dist = pd.DataFrame(dist, names1, names2)
        dist.to_csv(outfile)


def console_pearson():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=PEARSON_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('counts1',
                        help='Full path of a count file produced by kmer_counts.py.')
    parser.add_argument('counts2',
                        help=('Full path of a second count file produced by kmer_counts.py. '
                              'This can be the same path as the first counts file.'))
    parser.add_argument('-o', '--outfile', default='pearson.seekr',
                        help='Path of file to save similarities to.')
    parser.add_argument('-bi', '--binary_input', action='store_true',
                        help='Select if the input will be a .npy file.')
    parser.add_argument('-bo', '--binary_output', action='store_true',
                        help='Select if output should be a .npy file.')
    args = _parse_args_or_exit(parser)
    _run_pearson(args.counts1, args.counts2, args.outfile,
                 args.binary_input, args.binary_output)


def _run_norm_vectors(fasta, mean_vector, std_vector, kmer):
    counter = BasicCounter(fasta, k=int(kmer))
    counter.get_counts()
    np.save(mean_vector, counter.mean)
    np.save(std_vector, counter.std)


def console_norm_vectors():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=NORM_VECTORS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='path to .fa file')
    parser.add_argument('-mv', '--mean_vector', default='mean.npy',
                        help='path to output mean vector')
    parser.add_argument('-sv', '--std_vector', default='std.npy',
                        help='path to output standard deviation vector')
    parser.add_argument('-k', '--kmer', default=6,
                        help='length of kmers you want to count')
    args = _parse_args_or_exit(parser)
    _run_norm_vectors(args.fasta, args.mean_vector, args.std_vector, args.kmer)


def _run_graph(adj, gml_path, csv_path, louvain, limit, resolution, n_comms, seed):
    # Note: This function is separated for testing purposes.
    leiden = not louvain
    if seed is not None:
        seed = int(seed)
    maker = graph.Maker(adj, gml_path, csv_path, leiden, limit, resolution, n_comms, seed)
    maker.make_gml_csv_files()


def console_graph():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=GRAPH_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('adj',
                        help='Path to either .csv or .npy file, representing adjacency matrix')
    parser.add_argument('-g', '--gml_path', default=None,
                        help='Path to output graph file in .gml format')
    parser.add_argument('-c', '--csv_path', default=None,
                        help='Path to output community file in .csv format')
    parser.add_argument('-l', '--louvain',  action='store_true',
                        help='If set, use Louvain for community detection instead of Leiden.')
    parser.add_argument('-t', '--threshold', default=0,
                        help=('Value for thresholding adjacency matrix. '
                              'Below this limit, all edges are 0.'))
    parser.add_argument('-r', '--resolution', default=1,
                        help=' Resolution parameter for community detection algorithm')
    parser.add_argument('-n', '--n_comms', default=5,
                        help='Number of communities to find. This does not count a null community.')
    parser.add_argument('-s', '--seed', default=None)
    args = _parse_args_or_exit(parser)
    _run_graph(args.adj, args.gml_path, args.csv_path, args.louvain, args.threshold,
               args.resolution, args.n_comms, args.seed)


def _run_gen_rand_rnas(fasta, mean_vector, std_vector, kmer):
    # Note: This function is separated for testing purposes.
    pass


def console_gen_rand_rnas():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=GEN_RAND_RNAS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='path to .fa file')
    parser.add_argument('out_fasta', help='path to new .fa file')
    parser.add_argument('-k', '--kmer', default=1,
                        help='length of kmers you want to conserve')
    parser.add_argument('-m', '--mutations', default=0,
                        help='number of SNP mutations to make in RNA')
    parser.add_argument('-g', '--group', action='store_false',
                        help='set to concatenate RNAs before shuffling and mutating.')
    args = _parse_args_or_exit(parser)
    _run_norm_vectors(args.fasta, args.out_fasta, args.mutations, args.group)


def console_seekr_help():
    intro = ('Welcome to SEEKR! \n'
             'Below is a description of all SEEKR commands.\n'
             'For additional help see the README at: \n'
             'https://github.com/CalabreseLab/seekr.\n\n')
    print(intro)
    cmds = ['seekr_download',
            'seekr_norm_vectors',
            'seekr_kmer_counts',
            'seerk_pearson',
            'seekr_graph']
    docs = [DOWNLOAD_GENCODE_DOC, NORM_VECTORS_DOC, KMER_COUNTS_DOC, PEARSON_DOC, GRAPH_DOC]
    for c, d in zip(cmds, docs):
        print(f"{'='*20}\n{c}\n{'='*20}\n{d}")
    conclusion = ('To see a full description of flags and defaults, '
                  'run any of the commands listed above, without any parameters '
                  '(e.g. "$ seekr_graph").')
