import sys
import argparse
import numpy as np
import pandas as pd

from seekr import fasta
from seekr import graph	
from seekr import pearson
from seekr.kmer_counts import BasicCounter, Log2

# TODO (Dan) fix names
from seekr.pwm import CountsWeighter

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
    $ seekr_canonical_gencode rnas.fa rnas01.fa

If you want to specifically find transcripts with the ending 001:
    $ seekr_canonical_gencode rnas.fa rnas01.fa -z 2

To enforce one isoform per ENSG id (specifically, the smallest 01 isoform):
    $ seekr_canonical_gencode rnas.fa rnas01_1per.fa -u

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
The default settings take a .fa file and produce a labeld csv file:
    $ seekr_kmer_counts rnas.fa -o out.csv

To get a compact and efficient .npy file, set the binary flag:
    $ seekr_kmer_counts rnas.fa -o out.npy -b

You can change also change the size of the kmer you're using, and prevent normalization:
    $ seekr_kmer_counts rnas.fa -o out.csv -k 4 -uc -us -nl

If you ever do not want labels on a csv file:
    $ seekr_kmer_counts rnas.fa -o out.csv -rl

Notes
-----
For more sophisticated options, you cannot use the command line, but need python instead.

To pass --log 1 argument for pre-zscore log-transform of k-mer counts, seekr_norm_vectors MUST be
    run with the -cl flag. This log transforms the reference counts for appropriate mean and std calcs

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


VISUALIZE_DISTRO_DOC = """
Description
-----------
Generate an image showing the distribution of all Pearson r-values.
This can be useful for determining a threshold for the adjacency matrix.

Examples
--------
You must pass an adjacency matrix and an output path.
    $ seekr_visualize_distro adj.csv adj.pdf

For large arrays, it's likely sufficient to visualize a portion of the adjacency matrix.
You can pass a float between 0 and 1 to the `-s` flag:
    $ seekr_visualize_distro adj.csv adj.pdf -s .1

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

If pre-zscore log transform is desired, you must pass the -cl flag to log transform
    the reference k-mer counts 

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

Examples
--------
The two required positional arguments are an input and output path to fasta files.
This will shuffle the nucleotides for each RNA:
    $ seekr_gen_rand_rnas rnas.fa rnas_rand.fa

To conserve kmer content for k > 1, choose a different kmer size:
    $ seekr_gen_rand_rnas rnas.fa rnas_rand.fa -k 2

Shuffling kmers is random. To reproduce output between runs, set a seed:
    $ seekr_gen_rand_rnas rnas.fa rnas_rand.fa -s 0

It may be useful to conserve the kmer content of the whole fasta file.
Setting the `--group` flag loses conservation of individual sequences,
in preference for producing RNAs with a kmer frequency equal to background frequency.
Using `--group` still produces the same number of output RNAs as input.
Note: this may segfault on large fasta files.
    $ seekr_gen_rand_rnas rnas.fa rnas_rand.fa -k 2 -g

In some cases, it is useful to have more fine-grained control over final kmer content.
Ex: when conserving large kmers, it may be impossible to shuffle shorter seqs.
Ex: if you want to produce a sequence with an exact Pearson's r-value to the original.
A number of random SNP mutations can be made in addition to shuffling.
Use the --mutation flag to set the approximate number of SNP mutations.
Note: Because the new nucleotide may be the same as the old, -m is approximate.
    $ seekr_gen_rand_rnas rnas.fa rnas_rand.fa -m 100

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""


GRAPH_DOC = """
Description
-----------
Find communities of transcripts from an adjacency matrix.

Examples
--------
Default setting accept a csv file, and a threshold value.
The csv file should be the product of seekr_pearson, or some other adjacency matrix.
The threshold is the value below which edges are removed from the graph.
seekr_pearson_distro can be run to suggest a value for the threshold.
A gml file contain the graph and communities will be produced.
    $ seekr_graph adj.csv .13 -g graph.gml

For a cleaner csv file of just community information:
    $ seekr_graph adj.csv .5 -g graph.gml -c communities.csv

To change the resolution parameter (gamma) for louvain/leidenalg:
    $ seekr_graph adj.csv .1 -g graph.gml -r 1.5

To change the cap of the number of communities found, and set the seed:
    $ seekr_graph adj.csv .1 -g graph.gml -n 10 -s 0

Numpy files are also valid input:
    $ seekr_graph adj.npy .1 -g graph.gml

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""


PWM_DOC = """
Description
-----------
Weight kmer profiles by protein binding PWMs to infer protein binding likelihood.

Examples
--------
A standard run of this tool needs three things:
1. A directory of PWM files
2. A counts file (produced from seekr_kmer_counts)
3. An output path
    $ seekr_pwm path/to/pwms/ kc_out.csv -o pwm_weight_sums.csv

Numpy files can also be passed as input, but .csv files are the only output:
    $ seekr_pwm path/to/pwms/ kc_out.npy -o pwm_weight_sums.csv

The kmer size can also be passed.
It should match the counts file.
Unlike most other seekr tools k=5 is the default for this tool.
    $ seekr_pwm path/to/pwms/ kc_out.npy -k 6 -o pwm_weight_sums.csv

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""


DOMAIN_PEARSON_DOC = """
Description
-----------
# Find domains of similarity between query transcripts and tiles of target transcripts.

Examples
--------
This tool requires several pieces of data:
1. A fasta file containing query sequences.
2. A second fasta file containing target sequences which will be tiled.
3. A mean vector for normalization (e.g. from `seekr_norm_vectors`).
4. A std vector for standardization (e.g. from `seekr_norm_vectors`).

For brevity in the documentation below,
we will assume that these required data have been stored in a variable:
    $ REQUIRED="queries.fa targets.fa mean.npy std.npy"

To see the r-values, pass a location for storing them in a csv file.
    $ seekr_domain_pearson $REQUIRED -r r_values.csv

Intepretation of r-value elements can be aided by viewing them as percentiles.
If you want percentiles, you must also pass a reference fasta path:
    $ seekr_domain_pearson $REQUIRED -r r_values.csv -p percentiles.csv -rp reference.fa

Parameters you might pass to `seekr_kmer_counts` can also be passed.
If you change --kmer, ensure that your mean.npy and std.npy files match:
    $ seekr_domain_pearson $REQUIRED -r r_values.csv -nl -k 5

You can also change the size of the domain,
and how far you slide along the target sequence before creating another domain:
    $ seekr_domain_pearson $REQUIRED -r r_values.csv -w 1200 -s 150

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


def _run_canonical_gencode(in_fasta, out_fasta, zeros, unique_per_gene):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    maker = fasta.Maker(in_fasta, out_fasta)
    maker.filter1(zeros, unique_per_gene)


def console_canonical_gencode():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=CANONICAL_GENCODE_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_fasta', help='Full path of fasta file.')
    parser.add_argument('out_fasta', help='Full path of filtered fasta file.')
    parser.add_argument('-z', '--zeros', default=1,
                        help='Number of zeroes needed to be considered canonical.')
    parser.add_argument('-u', '--unique_per_gene', action='store_true',
                        help='Set to enforce a limit of one isoform per ENSG id.')
    args = _parse_args_or_exit(parser)
    _run_canonical_gencode(args.in_fasta, args.out_fasta, args.zeros, args.unique_per_gene)


def _run_kmer_counts(fasta, outfile, kmer, binary, centered, standardized,
                     log2, remove_labels, mean_vector, std_vector, alphabet):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    mean = mean_vector or centered
    std = std_vector or standardized
    label = not remove_labels
    counter = BasicCounter(fasta, outfile, kmer, binary,
                           mean, std, Log2[log2], label=label, alphabet=alphabet)
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
                        help='Set if output should be a .npy file.')
    parser.add_argument('-uc', '--uncentered', action='store_false',
                        help='Set if output should not have the mean subtracted.')
    parser.add_argument('-us', '--unstandardized', action='store_false',
                        help='Set if output should not be divided by the standard deviation.')
    parser.add_argument('-l', '--log2', default=Log2.post.name,
                        choices=[l2.name for l2 in Log2],
                        help='Decided if and when to log transform counts')
    parser.add_argument('-rl', '--remove_labels', action='store_true',
                        help='Set to save without index and column labels.')
    parser.add_argument('-mv', '--mean_vector', default=None,
                        help='Optional path to mean vector numpy file.')
    parser.add_argument('-sv', '--std_vector', default=None,
                        help='Optional path to std vector numpy file.')
    parser.add_argument('-a', '--alphabet', default='AGTC',
                        help='Valid letters to include in kmer.')
    args = _parse_args_or_exit(parser)
    _run_kmer_counts(args.fasta, args.outfile, int(args.kmer), args.binary, args.uncentered,
                     args.unstandardized, args.log2, args.remove_labels, args.mean_vector,
                     args.std_vector, args.alphabet)


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
        pearson.pearson(counts1, counts2, outfile=outfile)
    else:
        dist = pearson.pearson(counts1, counts2)
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
                        help='Set if the input will be a .npy file.')
    parser.add_argument('-bo', '--binary_output', action='store_true',
                        help='Set if output should be a .npy file.')
    args = _parse_args_or_exit(parser)
    _run_pearson(args.counts1, args.counts2, args.outfile,
                 args.binary_input, args.binary_output)


def _run_visualize_distro(adj, out_path, sample):
    if sample is not None:
        sample = float(sample)
    mean, std = pearson.visualize_distro(adj, out_path, sample)
    print('Mean: ', mean)
    print('Std. Dev.: ', std)


def console_visualize_distro():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=VISUALIZE_DISTRO_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('adj',
                        help='Path to either .csv or .npy file, representing adjacency matrix')
    parser.add_argument('out_path', help='Full path of a output image.')
    parser.add_argument('-s', '--sample', default=None,
                        help='Float representing random portion of adj to visualize.')
    args = _parse_args_or_exit(parser)
    _run_visualize_distro(args.adj, args.out_path, args.sample)


def _run_norm_vectors(fasta, mean_vector, std_vector, log2, kmer):
    counter = BasicCounter(fasta, k=int(kmer), log2=Log2[log2])
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
    parser.add_argument('-l', '--log2', default=Log2.post.name,
                        choices=[l2.name for l2 in Log2],
                        help='Decided if and when to log transform counts')
    parser.add_argument('-k', '--kmer', default=6,
                        help='length of kmers you want to count')
    args = _parse_args_or_exit(parser)
    _run_norm_vectors(args.fasta, args.mean_vector, args.std_vector, args.log2, int(args.kmer))


def _run_graph(adj, threshold, gml_path, csv_path, louvain, resolution, n_comms, seed):
    # Note: This function is separated for testing purposes.
    leiden = not louvain
    if seed is not None:
        seed = int(seed)
    maker = graph.Maker(adj, gml_path, csv_path, leiden, threshold, resolution, n_comms, seed)
    maker.make_gml_csv_files()


def console_graph():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=GRAPH_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('adj',
                        help='Path to either .csv or .npy file, representing adjacency matrix')
    parser.add_argument('threshold', help=('Value for thresholding adjacency matrix. '
                                           'Below this limit, all edges are 0.'))
    parser.add_argument('-g', '--gml_path', default=None,
                        help='Path to output graph file in .gml format')
    parser.add_argument('-c', '--csv_path', default=None,
                        help='Path to output community file in .csv format')
    parser.add_argument('-l', '--louvain',  action='store_true',
                        help='If set, use Louvain for community detection instead of Leiden.')
    parser.add_argument('-r', '--resolution', default=1,
                        help=' Resolution parameter for community detection algorithm')
    parser.add_argument('-n', '--n_comms', default=5,
                        help='Number of communities to find. This does not count a null community.')
    parser.add_argument('-s', '--seed', default=None,
                        help='An integer to create reproducible results between runs.')
    args = _parse_args_or_exit(parser)
    _run_graph(args.adj, float(args.threshold), args.gml_path, args.csv_path, args.louvain,
               float(args.resolution), int(args.n_comms), args.seed)


def _run_gen_rand_rnas(in_fasta, out_fasta, kmer, mutations, seed, group):
    # Note: This function is separated for testing purposes.
    # TODO do something with group?
    kmer = int(kmer)
    mutations = int(mutations)
    if seed is not None:
        seed = int(seed)
    individual = not group
    rand_maker = fasta.RandomMaker(in_fasta, out_fasta, kmer, mutations, seed, individual)
    rand_maker.synthesize_random()


def console_gen_rand_rnas():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=GEN_RAND_RNAS_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_fasta', help='path to .fa file')
    parser.add_argument('out_fasta', help='path to new .fa file')
    parser.add_argument('-k', '--kmer', default=1,
                        help='Length of kmers you want to conserve')
    parser.add_argument('-m', '--mutations', default=0,
                        help='Number of SNP mutations to make in RNA')
    parser.add_argument('-s', '--seed', default=None,
                        help='An integer to create reproducible results between runs.')
    parser.add_argument('-g', '--group', action='store_true',
                        help='Set to concatenate RNAs before shuffling and mutating.')
    args = _parse_args_or_exit(parser)
    _run_gen_rand_rnas(args.in_fasta, args.out_fasta, int(args.kmer), int(args.mutations),
                       args.seed, args.group)


def _run_pwms(pwm_dir, counts, kmer, out_path):
    # TODO (Dan) name this function
    # Note: This function is separated for testing purposes.
    counts_weighter = CountsWeighter(pwm_dir, counts, kmer, out_path)
    counts_weighter.run()


def console_pwm():
    # TODO (Dan) name this function
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=PWM_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('pwm_dir', help='Path to directory containing PWM files.')
    parser.add_argument('counts', help='Path to kmer_counts file.')
    parser.add_argument('-k', '--kmer', default=5,
                        help='Length of kmer.')
    parser.add_argument('-o', '--out_path',
                        help='Path to new csv file containing weighted count sums.')
    args = _parse_args_or_exit(parser)
    # TODO (Dan) update name
    _run_pwms(args.pwm_dir, args.counts, int(args.kmer), args.out_path)


def _run_domain_pearson(query_path, target_path, reference_path, mean, std, r_values,
                        percentiles, kmer, log2, window, slide):
    # Note: This function is separated for testing
    domain_pearson = pearson.DomainPearson(query_path, target_path, reference_path, r_values,
                                           percentiles, mean, std, log2, kmer, window, slide)
    domain_pearson.run()


def console_domain_pearson():
    assert sys.version_info[0] == 3, 'Python version must be 3.x'
    parser = argparse.ArgumentParser(usage=DOMAIN_PEARSON_DOC,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query_path',
                        help='Path to fa file containing transcripts of interest (e.g. Xist-2kb).')
    parser.add_argument('target_path', help=('Path to second fa file which will be tiled to find '
                                             'domains similar to query transcripts.'))
    parser.add_argument('mean', help='Path to npy file containing mean array for normalization.')
    parser.add_argument('std', help='Path to npy file containing std array for standardization.')
    parser.add_argument('-rp', '--reference_path', default=None,
                        help=('Path to third fasta file containing sequences to be used for '
                              'comparison when calculating percentile values of the r-values '
                              'between the query and targets (e.g. mouse transcriptome).'))
    parser.add_argument('-r', '--r_values', help='Path to new csv file for storing r-values.')
    parser.add_argument('-p', '--percentiles', help='Path to new csv file for storing percentiles.')
    parser.add_argument('-k', '--kmer', default=6,
                        help='Length of kmers you want to count.')
    parser.add_argument('-l', '--log2',
                        help='Choose a value of 1,2, or 3 for different log transformation options')
    parser.add_argument('-w', '--window', default=1000,
                        help=('Size of tile/domain to be created from target transcripts for '
                              'comparison against queries.'))
    parser.add_argument('-s', '--slide', default=100,
                        help=('Number of basepairs to move along target transcript before creating '
                              'another tile/domain.'))
    args = _parse_args_or_exit(parser)
    _run_domain_pearson(args.query_path, args.target_path, args.reference_path, args.mean, args.std,
                        args.r_values, args.percentiles, int(args.kmer), args.log2,
                        int(args.window), args.slide)


def console_seekr_help():
    intro = ('Welcome to SEEKR! \n'
             'Below is a description of all SEEKR commands.\n'
             'For additional help see the README at: \n'
             'https://github.com/CalabreseLab/seekr.\n\n')
    print(intro)
    cmds2doc = {'seekr_download_gencode': DOWNLOAD_GENCODE_DOC,
                'seekr_canonical_gencode': CANONICAL_GENCODE_DOC,
                'seekr_norm_vectors': NORM_VECTORS_DOC,
                'seekr_kmer_counts': KMER_COUNTS_DOC,
                'seekr_pearson': PEARSON_DOC,
                'seekr_visualize_distro': VISUALIZE_DISTRO_DOC,
                'seekr_graph': GRAPH_DOC,
                'seekr_gen_rand_rnas': GEN_RAND_RNAS_DOC,
                'seekr_pmw': PWM_DOC,
                'seekr_domain_pearson': DOMAIN_PEARSON_DOC}
    for c, d in cmds2doc.items():
        print(f"{'='*25}\n{c}\n{'='*25}\n{d}")
    conclusion = ('To see a full description of flags and defaults, '
                  'run any of the commands listed above, without any parameters '
                  '(e.g. "$ seekr_graph").')
    print(conclusion)
