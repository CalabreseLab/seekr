import sys
import argparse
import numpy as np
import pandas as pd

from seekr import fasta
from seekr import graph
from seekr import pearson
from seekr.kmer_counts import BasicCounter
from seekr import find_dist
from seekr import find_pval
from seekr import kmer_heatmap
from seekr import kmer_dendrogram
from seekr import kmer_leiden
from seekr import kmer_count_barplot
from seekr import kmer_msd_barplot
from seekr import kmer_textplot

from seekr.__version__ import __version__

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
    $ seekr_kmer_counts rnas.fa -o out.csv -k 4 -uc -us -l none

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

If pre-zscore log transform is desired, you must pass the `--log2 Log2.pre` flag to log transform
    the reference k-mer counts

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""

GEN_RAND_RNAS_DOC = """
Deprecated
----------
This function is deprecated by 2023. 
As the ushuffle python package is causing problem during installation, the function is removed from seekr package.
Please refer to fasta-shuffle-letters of the MEME Suite for generating random sequence while preserving kmer content.
https://meme-suite.org/meme/doc/fasta-shuffle-letters.html
fasta-shuffle-letters uses the same ushuffle algorithm in C platform. 

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
Find domains of similarity between query transcripts and tiles of target transcripts.

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
    $ seekr_domain_pearson $REQUIRED -r r_values.csv -l none -k 5

You can also change the size of the domain,
and how far you slide along the target sequence before creating another domain:
    $ seekr_domain_pearson $REQUIRED -r r_values.csv -w 1200 -s 150

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""


FIND_DIST_DOC = """
Description: 
Find the best fitted distribution to the input background sequences 

Details:
data to be fitted -- all possible pairwise kmer pearson correlation for the background sequences
can also choose to return the actual data by itself
the results (best fitted distribution or actual data) will be used to calculate p-values

Example:
calculate all possible pairwise pearson correlation scores of the background sequences -- all unique lncRNA of mouse vM25 from GENCODE
with kmer size of 4, log2 transform of 'Log2.post'
then subset 100000 of such data to fit models 'common10' with KS test as goodness of fit quantification
show progress bar for the model fitting step
plot the model fits against the actual data (modelfit.pdf) and save the fitted models to test_fitres.csv under current directory
    
    $ seekr_find_dist default -k 4 -l Log2.post -mdl common10 -sbt -sbs 100000 -fm -statm ks -pb -pf modelfit.pdf -o test

to not save the model fit plot change to --plotfit (that is not providing any value after --plotfit)
to save under other location change to -o /Users/username/Desktop/fitres.csv

to use user defined models list -mdl 'norm,expon,pareto'

minimal code with all settings to default
fitting the models:
    $ seekr_find_dist default -fm 
not fitting the models, directly output the background numpy array:
    $ seekr_find_dist default

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

FIND_PVAL_DOC = """
Description: 
calculte p values of the seekr.pearson correlation values for input sequence 1 vs input sequence 2 
p value is based on the output of find_dist, which is either a list of distributions or a npy array

Details:
this function connects the output of find_dist and the input sequnces of interests (2 input sequences)
given the background sequencs, find_dist calculate all possible pairwise seekr.pearson values
and then outputs either a list of fitted distributions or the npy array of the actual data 
find_pval firstly calculates the seekr.pearson of the two input sequences, which produces a correlation matrix 
find_pval then calculate the p value for each r value in the pearson correlation matrix based on the output of find_dist
the output of find_pval is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences

Example:
calculate p values for the seekr.pearson correlation between input sequences 1 and input sequences 2
with kmer size of 4, log2 transform of 'Log2.post', and corresponding normalization vectors saved as mean_4mer.npy and std_4mer.npy
with the bakcground distribution saved as a list of fitted distributions fitres.csv, among which the best fitted distribution is the first one
show progress bar for the p value calculation step and save the output to test_pval.csv under current directory
    
    $ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv -ft distribution -bf 1 -o test -pb

to save under other location change to -o /Users/username/Desktop/pval.csv

minimal code with all settings to default

    $ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_HEATMAP_DOC = """
Description:
Heatmap with both rows and columns dendrograms for input dataframe
visualize results of either pearson correlation r values of seekr.pearson or p-values of find_pval

Details:
Customizeable heatmap for easier visualization of the results of seekr.pearson or find_pval
takes in a dataframe with row and column names
performs hierarchial clustering on both rows and columns, and then reorder the cols and rows based on the dendrograms
can also use kmer_dendrgram to plot only the dendrograms (partial or full) to get a better idea of the clustering
for interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot and kmer_textplot for further analysis

Example:
plot heatmap for the input pval.csv file with a change of color at 0.05 (#ffffff, white)
p value 0 will be colored #1b7837 (green) and p value 1 will be colored #c51b7d (purple)
the distance metric is correlation and the linkage method is complete
the plot will be saved as test_kmer_heatmap.pdf under /Users/username/Desktop/
    
    $ seekr_kmer_heatmap pval.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o /Users/username/Desktop/test -hf pdf -hd 300

to save under the current directory change to -o test

minimal code with all settings to default

    $ seekr_kmer_heatmap pval.csv 0 1 

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_DENDROGRAM_DOC = """
Description:
Dendrpgram for the hierarchical clustering of the rows or columns for the input dataframe
visualize clustering results of either pearson correlation r values of seekr.pearson or p-values of find_pval

Details:
Customizeable dendrograms for easier visualization of the hierarchical clustering results of seekr.pearson or find_pval
takes in a dataframe with row and column names
performs hierarchial clustering on either rows or columns, is a better way to visualize the whole or partial clusters in kmer_heatmap
for interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot and kmer_textplot for further analysis

Example:
plot dendrogram for the rows of the input pval.csv file
with the distance metric of correlation and the linkage method of complete
the plot will be saved as test_kmer_dendrogram_row.pdf under current directory
    
    $ seekr_kmer_dendrogram pval.csv -dd row -distm correlation -linkm complete -ph 8 -hratio 0.5 -lfs 16 -o test -pf pdf -d 300

to save under other location change to -o /Users/username/Desktop/test

minimal code with all settings to default

    $ seekr_kmer_dendrogram pval.csv

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_LEIDEN_DOC = """
Description:
plot Leiden community network for input fasta seqeunces
it works better with small number of nodes (less than 50)
for better visualization and customization, please set savecsv to True and use Gephi with the saved nodes and edges files

Details:
take as input a fasta file with multiple sequences
calculate sequences distance matrix as seekr pearson correlation of the fasta file to itself
use Leiden community to call cluster
edge is weighted with the pearson correlation results
darker and thicker edges have higher weights
layout the nodes with networkx spring layout
the network is undirected

Example:
perform seekr pearson correlation on ldseq.fa with kmer size of 4 and corresponding normalization vectors: mean_4mer.npy and std_4mer.npy
correlation r values less than 0.1 will be set to 0
then perform Leiden community detection on the correlation matrix with algorithm RBERVertexPartition and resolution of 1.0
seed is not set so that the community assignment results could vary each time the code is run
the plot edge will be gradient colored based on the pearson correlation values
the plot will be saved as test_gradient_leiden_network.pdf under current directory
the corresponding nodes and edges files will be saved as test_nodes_leiden.csv and test_edges_leiden.csv under current directory
    
    $ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -a RBERVertexPartition -r 1.0 -pco 0.1 -sd -ec gradient -et 0.1 -lfs 12 -o test -s

to save under other location change to -o /Users/username/Desktop/test
to supress plotting, change to --outputname (that is not providing any value after --outputname or -o)
to set seed to get reproducible community assignments change -sd to -sd 1

minimal code with all settings to default

    $ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -s

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_COUNT_BARPLOT_DOC = """
Description:
barplot of the transformed or raw z-score for kmer words of the input sequences
limit input to 10 sequences -- hard to choose more than 10 distinct colors

Details:
plot the z-score that is the output of BasicCounter, where user can define whether and how to do log2 transform
order kmer words by the summed difference from the mean among all sequences, can be descending or ascending
user define to plot the top x kmer words (could be messy if there are too many words)
x axis label is the kmer words

Example:
Perform BasicCounter to count the kmers of sequences from test.fa file, with kmer size of 4 and corresponding normalization vectors: mean_4mer.npy and std_4mer.npy
then perform log2 post transform of the z score and arrange the kmer words by the summed difference from the mean among all sequences in ascending order, 
meaning the kmer words with the smallest summed difference from the mean among all sequences will be plotted first
and then only plot the top 10 kmer words into test_kmer_count_barplot.pdf under current directory
    
    $ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4 -l Log2.post -sm ascending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o test -pf pdf -d 300

to save under other location change to -o /Users/username/Desktop/test

minimal code with all settings to default

    $ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4 

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_MSD_BARPLOT_DOC = """
Description:
barplot of the mean of the transformed or raw z-score for kmer words across input sequences
error bar is the standard deviation of the transformed or raw z-score for kmer words across input sequences

Details:
the z-score that is the output of BasicCounter, where user can define whether and how to do log2 transform
get the mean and standard deviation (sd) of the z-score for each kmer word across input sequences
can choose to order kmer words by mean or sd, can either be descending or ascending
user define to plot the top x kmer words (could be messy if there are too many words)
x axis label is the kmer words

Example:
Perform BasicCounter to count the kmers of sequences from test.fa file, with kmer size of 4 and corresponding normalization vectors: mean_4mer.npy and std_4mer.npy
then perform log2 post transform of the z score and calculate the mean and std of the z score for each kmer word
then arrange the kmer words by the mean of the transformed z score in descending order, meaning the kmer words with the largest mean of the z score will be plotted first
and then only plot the top 10 kmer words into test_kmer_msd_barplot.pdf under current directory
    
    $ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o test -pf pdf -d 300

to save under other location change to -o /Users/username/Desktop/test

minimal code with all settings to default

    $ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4 

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_TEXTPLOT_DOC = """
Description:
highlight the input kmer words in the input sequences

Details:
plot 2 input sequences and several interested words (max 10)
limit the number of interested words to 10 as it is hard to distinguish more than 10 colors
align the sequence at the beginning and label sequence positions on the bottom
highlight in colors and bold the words of interest
for multiple words, if the words overlap, the overlapped characters will be highlighted in the color of the first word in the list
therefore arrange the words in the list based on the priority of the words, with the most important word at the beginning of the list

Example:
Plot and align the two input sequences with numbered positions. Highlight the words of interest in corresponding colors at all matched locations along the two sequences
the words of interest are 'ATTA' in '#d62728', 'AAAA' in '#e377c2', and 'ACTC' in '#ff7f0e'
wrap length is 60 characters and the textplot is saved as test_kmer_textplot.pdf under current directory
    
    $ seekr_kmer_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC' -cv '#d62728,#e377c2,#ff7f0e' -wl 60 -cs 0.1 -ls 0.2 -sfs 42 -nfs 40 -o test -pf pdf -d 300

to save under other location change to -o /Users/username/Desktop/test

minimal code with all settings to default

$ seekr_kmer_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC'

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

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
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_COUNTS_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("biotype", help=("Name of Genocde set to download. " "Must be one of ('all', 'pc', 'lncRNA')."))
    parser.add_argument(
        "-s", "--species", default="human", help=" Name of species. Must be one of: ('human' or 'mouse')."
    )
    parser.add_argument(
        "-r",
        "--release",
        default=None,
        help=("Name of specific release to download (e.g. 'M5'). " "If None, download latest release."),
    )
    parser.add_argument(
        "-o", "--out_path", default=None, help="Path to location for fasta file. Default will save by release name."
    )
    parser.add_argument(
        "-z", "--zip", action="store_false", help="Set if you do not want to gunzip fasta file after downloading."
    )
    args = _parse_args_or_exit(parser)
    _run_download_gencode(args.biotype, args.species, args.release, args.out_path, args.zip)


def _run_canonical_gencode(in_fasta, out_fasta, zeros, unique_per_gene):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    maker = fasta.Maker(in_fasta, out_fasta)
    maker.filter1(zeros, unique_per_gene)


def console_canonical_gencode():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(
        usage=CANONICAL_GENCODE_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("in_fasta", help="Full path of fasta file.")
    parser.add_argument("out_fasta", help="Full path of filtered fasta file.")
    parser.add_argument("-z", "--zeros", default=1, help="Number of zeroes needed to be considered canonical.")
    parser.add_argument(
        "-u", "--unique_per_gene", action="store_true", help="Set to enforce a limit of one isoform per ENSG id."
    )
    args = _parse_args_or_exit(parser)
    _run_canonical_gencode(args.in_fasta, args.out_fasta, args.zeros, args.unique_per_gene)


def _run_kmer_counts(
    fasta, outfile, kmer, binary, centered, standardized, log2, remove_labels, mean_vector, std_vector, alphabet
):
    # Note: This function is separated from console_kmer_counts for testing purposes.
    mean = mean_vector or centered
    std = std_vector or standardized
    label = not remove_labels
    counter = BasicCounter(fasta, outfile, kmer, binary, mean, std, log2, label=label, alphabet=alphabet)
    counter.make_count_file()


def console_kmer_counts():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_COUNTS_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", help="Full path of fasta file.")
    parser.add_argument("-o", "--outfile", default="counts.seekr", help="Name of file to save counts to.")
    parser.add_argument("-k", "--kmer", default=6, help="Length of kmers you want to count.")
    parser.add_argument("-b", "--binary", action="store_true", help="Set if output should be a .npy file.")
    parser.add_argument(
        "-uc", "--uncentered", action="store_false", help="Set if output should not have the mean subtracted."
    )
    parser.add_argument(
        "-us",
        "--unstandardized",
        action="store_false",
        help="Set if output should not be divided by the standard deviation.",
    )
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="Decided if and when to log transform counts",
    )
    parser.add_argument(
        "-rl", "--remove_labels", action="store_true", help="Set to save without index and column labels."
    )
    parser.add_argument("-mv", "--mean_vector", default=None, help="Optional path to mean vector numpy file.")
    parser.add_argument("-sv", "--std_vector", default=None, help="Optional path to std vector numpy file.")
    parser.add_argument("-a", "--alphabet", default="AGTC", help="Valid letters to include in kmer.")
    args = _parse_args_or_exit(parser)
    _run_kmer_counts(
        args.fasta,
        args.outfile,
        int(args.kmer),
        args.binary,
        args.uncentered,
        args.unstandardized,
        args.log2,
        args.remove_labels,
        args.mean_vector,
        args.std_vector,
        args.alphabet,
    )


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
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=PEARSON_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("counts1", help="Full path of a count file produced by kmer_counts.py.")
    parser.add_argument(
        "counts2",
        help=(
            "Full path of a second count file produced by kmer_counts.py. "
            "This can be the same path as the first counts file."
        ),
    )
    parser.add_argument("-o", "--outfile", default="pearson.seekr", help="Path of file to save similarities to.")
    parser.add_argument("-bi", "--binary_input", action="store_true", help="Set if the input will be a .npy file.")
    parser.add_argument("-bo", "--binary_output", action="store_true", help="Set if output should be a .npy file.")
    args = _parse_args_or_exit(parser)
    _run_pearson(args.counts1, args.counts2, args.outfile, args.binary_input, args.binary_output)


def _run_visualize_distro(adj, out_path, sample):
    if sample is not None:
        sample = float(sample)
    mean, std = pearson.visualize_distro(adj, out_path, sample)
    print("Mean: ", mean)
    print("Std. Dev.: ", std)


def console_visualize_distro():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=VISUALIZE_DISTRO_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("adj", help="Path to either .csv or .npy file, representing adjacency matrix")
    parser.add_argument("out_path", help="Full path of a output image.")
    parser.add_argument("-s", "--sample", default=None, help="Float representing random portion of adj to visualize.")
    args = _parse_args_or_exit(parser)
    _run_visualize_distro(args.adj, args.out_path, args.sample)


def _run_norm_vectors(fasta, mean_vector, std_vector, log2, kmer):
    counter = BasicCounter(fasta, k=int(kmer), log2=log2)
    counter.get_counts()
    np.save(mean_vector, counter.mean)
    np.save(std_vector, counter.std)


def console_norm_vectors():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=NORM_VECTORS_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", help="path to .fa file")
    parser.add_argument("-mv", "--mean_vector", default="mean.npy", help="path to output mean vector")
    parser.add_argument("-sv", "--std_vector", default="std.npy", help="path to output standard deviation vector")
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="Decided if and when to log transform counts",
    )
    parser.add_argument("-k", "--kmer", default=6, help="length of kmers you want to count")
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
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=GRAPH_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("adj", help="Path to either .csv or .npy file, representing adjacency matrix")
    parser.add_argument(
        "threshold", help=("Value for thresholding adjacency matrix. " "Below this limit, all edges are 0.")
    )
    parser.add_argument("-g", "--gml_path", default=None, help="Path to output graph file in .gml format")
    parser.add_argument("-c", "--csv_path", default=None, help="Path to output community file in .csv format")
    parser.add_argument(
        "-l", "--louvain", action="store_true", help="If set, use Louvain for community detection instead of Leiden."
    )
    parser.add_argument("-r", "--resolution", default=1, help=" Resolution parameter for community detection algorithm")
    parser.add_argument(
        "-n", "--n_comms", default=5, help="Number of communities to find. This does not count a null community."
    )
    parser.add_argument("-s", "--seed", default=None, help="An integer to create reproducible results between runs.")
    args = _parse_args_or_exit(parser)
    _run_graph(
        args.adj,
        float(args.threshold),
        args.gml_path,
        args.csv_path,
        args.louvain,
        float(args.resolution),
        int(args.n_comms),
        args.seed,
    )


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
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=GEN_RAND_RNAS_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_fasta", help="path to .fa file")
    parser.add_argument("out_fasta", help="path to new .fa file")
    parser.add_argument("-k", "--kmer", default=1, help="Length of kmers you want to conserve")
    parser.add_argument("-m", "--mutations", default=0, help="Number of SNP mutations to make in RNA")
    parser.add_argument("-s", "--seed", default=None, help="An integer to create reproducible results between runs.")
    parser.add_argument(
        "-g", "--group", action="store_true", help="Set to concatenate RNAs before shuffling and mutating."
    )
    args = _parse_args_or_exit(parser)
    _run_gen_rand_rnas(args.in_fasta, args.out_fasta, int(args.kmer), int(args.mutations), args.seed, args.group)


def _run_pwms(pwm_dir, counts, kmer, out_path):
    # TODO (Dan) name this function
    # Note: This function is separated for testing purposes.
    counts_weighter = CountsWeighter(pwm_dir, counts, kmer, out_path)
    counts_weighter.run()


def console_pwm():
    # TODO (Dan) name this function
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=PWM_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pwm_dir", help="Path to directory containing PWM files.")
    parser.add_argument("counts", help="Path to kmer_counts file.")
    parser.add_argument("-k", "--kmer", default=5, help="Length of kmer.")
    parser.add_argument("-o", "--out_path", help="Path to new csv file containing weighted count sums.")
    args = _parse_args_or_exit(parser)
    # TODO (Dan) update name
    _run_pwms(args.pwm_dir, args.counts, int(args.kmer), args.out_path)


def _run_domain_pearson(
    query_path, target_path, reference_path, mean, std, r_values, percentiles, kmer, log2, window, slide
):
    # Note: This function is separated for testing
    domain_pearson = pearson.DomainPearson(
        query_path, target_path, reference_path, r_values, percentiles, mean, std, log2, kmer, window, slide
    )
    domain_pearson.run()


def console_domain_pearson():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=DOMAIN_PEARSON_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("query_path", help="Path to fa file containing transcripts of interest (e.g. Xist-2kb).")
    parser.add_argument(
        "target_path",
        help=("Path to second fa file which will be tiled to find " "domains similar to query transcripts."),
    )
    parser.add_argument("mean", help="Path to npy file containing mean array for normalization.")
    parser.add_argument("std", help="Path to npy file containing std array for standardization.")
    parser.add_argument(
        "-rp",
        "--reference_path",
        default=None,
        help=(
            "Path to third fasta file containing sequences to be used for "
            "comparison when calculating percentile values of the r-values "
            "between the query and targets (e.g. mouse transcriptome)."
        ),
    )
    parser.add_argument("-r", "--r_values", help="Path to new csv file for storing r-values.")
    parser.add_argument("-p", "--percentiles", help="Path to new csv file for storing percentiles.")
    parser.add_argument("-k", "--kmer", default=6, help="Length of kmers you want to count.")
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="Decided if and when to log transform counts",
    )
    parser.add_argument(
        "-w",
        "--window",
        default=1000,
        help=("Size of tile/domain to be created from target transcripts for " "comparison against queries."),
    )
    parser.add_argument(
        "-s",
        "--slide",
        default=100,
        help=("Number of basepairs to move along target transcript before creating " "another tile/domain."),
    )
    args = _parse_args_or_exit(parser)
    _run_domain_pearson(
        args.query_path,
        args.target_path,
        args.reference_path,
        args.mean,
        args.std,
        args.r_values,
        args.percentiles,
        int(args.kmer),
        args.log2,
        int(args.window),
        args.slide,
    )


def _run_console_seekr_help(version):
    if version:
        print(__version__)
        sys.exit()

    intro = (
        f"Welcome to SEEKR! ({__version__})\n"
        "Below is a description of all SEEKR commands.\n"
        "For additional help see the README at: \n"
        "https://github.com/CalabreseLab/seekr.\n\n"
    )
    print(intro)
    cmds2doc = {
        "seekr_download_gencode": DOWNLOAD_GENCODE_DOC,
        "seekr_canonical_gencode": CANONICAL_GENCODE_DOC,
        "seekr_norm_vectors": NORM_VECTORS_DOC,
        "seekr_kmer_counts": KMER_COUNTS_DOC,
        "seekr_pearson": PEARSON_DOC,
        "seekr_visualize_distro": VISUALIZE_DISTRO_DOC,
        "seekr_graph": GRAPH_DOC,
        "seekr_gen_rand_rnas": GEN_RAND_RNAS_DOC,
        "seekr_pmw": PWM_DOC,
        "seekr_domain_pearson": DOMAIN_PEARSON_DOC,
        "find_dist": FIND_DIST_DOC,
        "find_pval": FIND_PVAL_DOC,
        "kmer_heatmap": KMER_HEATMAP_DOC,
        "kmer_dendrogram": KMER_DENDROGRAM_DOC,
        "kmer_leiden": KMER_LEIDEN_DOC,
        "kmer_count_barplot": KMER_COUNT_BARPLOT_DOC,
        "kmer_msd_barplot": KMER_MSD_BARPLOT_DOC,
        "kmer_textplot": KMER_TEXTPLOT_DOC,
    }
    for c, d in cmds2doc.items():
        print(f"{'='*25}\n{c}\n{'='*25}\n{d}")
    conclusion = (
        "To see a full description of flags and defaults, "
        "run any of the commands listed above, without any parameters "
        '(e.g. "$ seekr_graph").'
    )
    print(conclusion)


def console_seekr_help():
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--version", action="store_true", help="Print current version and exit.")
    args = parser.parse_args()
    _run_console_seekr_help(args.version)




def console_find_dist():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=FIND_DIST_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", 
                        help="full path of the fasta file contains background sequences or 'default' which uses mouse vM25 Long non-coding RNA unique transcript sequences.")
    parser.add_argument("-k", "--kmer", default=4, help="length of kmers you want to count.")
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="decided if and when to log transform counts",
        )
    parser.add_argument("-mdl", "--models", default='common10', help="groups of candidate models for fitting, options are: 'all','common10' (default), or a string of distributions separated by comma for example 'norm,expon,pareto'.")
    parser.add_argument(
        "-sbt", "--subsetting", action="store_true", help="whether or not to use a subset of the data for fitting or output."
        )
    parser.add_argument("-sbs", "--subset_size", default=100000, help="the size of the subset to use for fitting or output.")
    parser.add_argument(
        "-fm",
        "--fit_model",
        action="store_true",
        help="whether or not to fit the data to the distributions sepcified in --models.",
        )
    parser.add_argument(
        "-statm", 
        "--statsmethod", 
        default='ks', 
        choices=['ks', 'mse', 'aic', 'bic'],
        help="the stats method used to quantify the goodness of fit of the distributions vs data."
        )
    parser.add_argument(
        "-pb", "--progress_bar", action="store_true", help="whether or not to show the progress bar when fitting the data to the distributions."
        )
    parser.add_argument("-pf", "--plotfit", default=None, help="whether or not to plot the fitted distributions (red dash line) vs the actual data (blue solid line) for all distributions in --models, fitted_model_plots.pdf will be saved under current directory.")
    parser.add_argument("-o", "--outputname", default="test", help="path and name to save fitted results as a csv file.")
    args = _parse_args_or_exit(parser)

    if (args.models != 'common10' and args.models != 'all'):
        modelslist = args.models.split(',')
    else:
        modelslist = args.models

    find_dist.find_dist(
        args.fasta,
        int(args.kmer),
        args.log2,
        modelslist,
        args.subsetting,
        int(args.subset_size),
        args.fit_model,
        args.statsmethod,
        args.progress_bar,
        args.plotfit,
        args.outputname,
        )


def console_find_pval():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=FIND_PVAL_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("seq1file", help="the full path to the fasta file of input sequence 1.")
    parser.add_argument(
        "seq2file",
        help=(
            "the full path to the fasta file of input sequence 2."
            "This can be the same path as the seq1file."
        ),
        )
    parser.add_argument("mean_path", help="full path to the normalization mean vector numpy file.")
    parser.add_argument("std_path",  help="full path to the normalization std vector numpy file.")
    parser.add_argument("kmer", help="length of kmers you want to count, which should be consistent with the normalization mean and std vector numpy files.")
    parser.add_argument("fitres_file",  help="the full path to the output csv file of find_dist, which is either a list of distributions or a npy array.")
    parser.add_argument(
        "-ft", 
        "--fitres_type", 
        default='distribution', 
        choices=['distribution', 'npy'], 
        help="specify the type of fitres, which is either a list of distributions or a npy array."
        )
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="decided if and when to log transform counts",
        )
    parser.add_argument("-bf", "--bestfit", default=1, help="the index of user determined best fit distribution in the list of distributions returned by find_dist (only used when fitres is a list of distributions), using normal index starting from 1 (not 0).")
    parser.add_argument("-o", "--outputname", default='test', help="the full path to store the p value csv file.")
    parser.add_argument(
        "-pb", "--progress_bar", action="store_true", help="whether or not to show the progress bar during p value calculation."
        )
    args = _parse_args_or_exit(parser)

    if (args.fitres_type == 'distribution'):
        # Read the CSV file into a DataFrame
        fitres = pd.read_csv(args.fitres_file)
        # Convert the DataFrame back to a list of tuples
        fitres = [tuple(row) for row in fitres.values]
        # the params (3rd column) is a tuple but when reloaded this way is imported as a string
        # convert the params back to a tuple
        fitres = [(row[0], row[1], tuple(map(float, row[2][1:-1].split(',')))) for row in fitres]
    else:
        # Load the CSV data into a NumPy array
        fitres = np.loadtxt(args.fitres_file, delimiter=',')

    find_pval.find_pval(
        args.seq1file,
        args.seq2file,
        args.mean_path,
        args.std_path,
        int(args.kmer),
        fitres,
        args.log2,
        int(args.bestfit),
        args.outputname,
        args.progress_bar,
        )


def console_kmer_heatmap():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_HEATMAP_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("df_file", 
                        help="full path of the input csv file, which could be a csv of pearson correlation matrix or p value with row and column names.")
    parser.add_argument("datamin", help="the minimum possible value of the data. For r values, this is -1. For p values, this is 0.")
    parser.add_argument("datamax", help="the maximum possible value of the data. For r values, this is 1. For p values, this is 1.")
    parser.add_argument("-th", "--thresh_value", default=0.05, help="If using a 3 color palette, this corresponds to the middle color.")
    parser.add_argument(
        "-cr",
        "--color_range_str",
        default='#1b7837,#ffffff,#c51b7d',
        help="A comma-separated string of colors in hex format for the heatmap, has to include 2 or 3 hex colors."
        )
    parser.add_argument(
        "-distm", 
        "--distmetric", 
        default='correlation', 
        help="the distance metric to use for the dendrogram, common options are: 'euclidean', 'cityblock', 'correlation', 'cosine', 'jaccard', 'hamming'."
        )
    parser.add_argument(
        "-linkm", 
        "--linkmethod", 
        default="complete", 
        help="the linkage method to use for the dendrogram, common options are: 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'."
        )
    parser.add_argument("-wratio", "--hmapw_ratio", default=0.3, help="a ratio factor to control the heatmap width.")
    parser.add_argument("-hratio", "--hmaph_ratio", default=0.3, help="a ratio factor to control the heatmap height.")
    parser.add_argument("-xts", "--x_tick_size", default=16, help="the font size for the column labels.")
    parser.add_argument("-yts", "--y_tick_size", default=16, help="the font size for the column labels.")
    parser.add_argument("-cfs", "--cbar_font_size", default=16, help="the font size for the colorbar ticks.")
    parser.add_argument("-o", "--outputname", default="test", help="the path and name to save the output heatmap, will automatically be combined with the trailing part _kmer_heatmap.")
    parser.add_argument("-hf", "--hformat", default='pdf', help="the format of the output heatmap, other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'.")
    parser.add_argument("-hd", "--hdpi", default=300, help="the dpi of the output heatmap.")
    args = _parse_args_or_exit(parser)

    df=pd.read_csv(args.df_file, index_col=0)
    # Split the input hex color string into a list using ',' as the delimiter
    color_range = args.color_range_str.split(',')
    kmer_heatmap.kmer_heatmap(
        df,
        int(args.datamin),
        int(args.datamax),
        float(args.thresh_value),
        color_range,
        args.distmetric,
        args.linkmethod,
        float(args.hmapw_ratio),
        float(args.hmaph_ratio),
        int(args.x_tick_size),
        int(args.y_tick_size),
        int(args.cbar_font_size),
        args.outputname,
        args.hformat,
        int(args.hdpi),
        )
    

def console_kmer_dendrogram():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_DENDROGRAM_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("df_file", 
                        help="full path of the input csv file, which could be a csv of pearson correlation matrix or p value with row and column names.")
    parser.add_argument(
        "-dd", 
        "--dendro_direct", 
        default='row', 
        choices=['row', 'column'],
        help="the direction to perform hierarchical clustering.")
    parser.add_argument(
        "-distm", 
        "--distmetric", 
        default='correlation', 
        help="the distance metric to use for the dendrogram, common options are: 'euclidean', 'cityblock', 'correlation', 'cosine', 'jaccard', 'hamming'."
        )
    parser.add_argument(
        "-linkm", 
        "--linkmethod", 
        default="complete", 
        help="the linkage method to use for the dendrogram, common options are: 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'."
        )
    parser.add_argument("-ph", "--plot_ht", default=8, help="height of the dendrogram plot.")
    parser.add_argument("-hratio", "--wd_ratio", default=0.5, help="a ratio factor to control the dendrgram width.")
    parser.add_argument("-lfs", "--leaf_font_size", default=16, help="the font size of the leaves or labels.")
    parser.add_argument("-o", "--outputname", default="test", help="the path and name to save the output dendrogram, will automatically be combined with the trailing part _kmer_dendrogram_row/column depending on --dendro_direct.")
    parser.add_argument("-pf", "--pformat", default='pdf', help="the format of the output dendrogram, other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'.")
    parser.add_argument("-d", "--pdpi", default=300, help="the dpi of the output dendrogram.")
    args = _parse_args_or_exit(parser)

    df=pd.read_csv(args.df_file, index_col=0)
    kmer_dendrogram.kmer_dendrogram(
        df,
        args.dendro_direct,
        args.distmetric,
        args.linkmethod,
        int(args.plot_ht),
        float(args.wd_ratio),
        int(args.leaf_font_size),
        args.outputname,
        args.pformat,
        int(args.pdpi),
        )
    


def console_kmer_leiden():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_LEIDEN_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", 
                        help="full path of the input fasta sequence file, please ensure all sequences have unique headers, this will be used for node label.")
    parser.add_argument("mean_path", help="full path to the normalization mean vector numpy file.")
    parser.add_argument("std_path",  help="full path to the normalization std vector numpy file.")
    parser.add_argument("kmer", help="length of kmers you want to count, which should be consistent with the normalization mean and std vector numpy files.")
    parser.add_argument(
        "-a", 
        "--algo", 
        default='RBERVertexPartition', 
        choices=['ModularityVertexPartition', 'RBConfigurationVertexPartition','RBERVertexPartition', 'CPMVertexPartition', 'SurpriseVertexPartition','SignificanceVertexPartition'],
        help="the Leiden algorithm used for calling communities or clusters.")
    parser.add_argument("-r", "--rs", default=1.0, help="the resolution parameter used in Leiden algorithm, larger value leads to more clusters.")
    parser.add_argument("-pco", "--pearsoncutoff", default=0.0, help="set to 0 pearson correlation r values that are less than the pearsoncutoff.")
    parser.add_argument("-sd", "--seed", default=None, help="set to an integer to get reproducible community assignments.")
    parser.add_argument(
        "-ec", 
        "--edgecolormethod", 
        default='gradient', 
        choices=['gradient', 'threshold'],
        help="the method used to set edge color, 'gradient' -- gradient scale of grey for egdes based on edge weights (pearson correlation); 'threshold' -- edge color thresholded on edge weights, smaller than threshold is grey, otherwise black."
        )
    parser.add_argument("-et", "--edgethreshold", default=0.1, help="the threshold used in --edgecolormethod 'threshold'.")
    parser.add_argument("-lfs", "--labelfontsize", default=12, help="the font size of the node label.")
    parser.add_argument(
        "-o", 
        "--outputname", 
        default="test", 
        help=(
            "the name and path of your output file."
            "leave blank to suppress the plotting"
            "this will broadcast to all saved images (with trailing part _gradient_leiden_network.pdf or _threshold_leiden_network.pdf) and csv files (with trailing part _edges_leiden.csv or _nodes_leiden.csv)."
            )
            )
    parser.add_argument(
        "-s", 
        "--savecsv", 
        action="store_true", 
        help=(
            "whether to save the nodes and edges csv files." 
            "The saved edges and nodes csv files are readily importable to Gephi for a better network visualization."
            )
            )
    args = _parse_args_or_exit(parser)

    kmer_leiden.kmer_leiden(
        args.fasta,
        args.mean_path,
        args.std_path,
        int(args.kmer),
        args.algo,
        float(args.rs),
        float(args.pearsoncutoff),
        int(args.seed)
        args.edgecolormethod,
        float(args.edgethreshold),
        int(args.labelfontsize),
        args.outputname,
        args.savecsv,
        )



def console_kmer_count_barplot():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_COUNT_BARPLOT_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", 
                        help="full path of the input fasta sequence file with unique headers, limit to the first 10 sequences or less of the input file.")
    parser.add_argument("mean_path", help="full path to the normalization mean vector numpy file.")
    parser.add_argument("std_path",  help="full path to the normalization std vector numpy file.")
    parser.add_argument("kmer", help="length of kmers you want to count, which should be consistent with the normalization mean and std vector numpy files.")
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="decided if and when to log transform counts",
        )
    parser.add_argument(
        "-sm", 
        "--sortmethod", 
        default='ascending', 
        choices=['ascending', 'descending'],
        help="how to sort the kmer words based on summed difference from the mean among all sequences."
        )
    parser.add_argument("-tn", "--topkmernumber", default=10, help="the number of sorted kmer words to plot, if input number is more than total number of kmer words, plot all words.")
    parser.add_argument("-xls", "--xlabelsize", default=20, help="x axis label font size.")
    parser.add_argument("-yls", "--ylabelsize", default=20, help="y axis label font size.")
    parser.add_argument("-xts", "--xticksize", default=20, help="x tick label font size.")
    parser.add_argument("-yts", "--yticksize", default=20, help="y tick label font size.")
    parser.add_argument("-ls", "--lengendsize", default=12, help="legend font size.")
    parser.add_argument(
        "-o", 
        "--outputname", 
        default="test", 
        help=(
            "the name and path of your output barplot."
            "a trailing part '_kmer_count_barplot' is then automatically added to the name."
            )
            )
    parser.add_argument("-pf", "--pformat", default='pdf', help="the format of the output barplot, other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'.")
    parser.add_argument("-d", "--pdpi", default=300, help="the dpi of the output barplot.")
    args = _parse_args_or_exit(parser)

    kmer_count_barplot.kmer_count_barplot(
        args.fasta,
        args.mean_path,
        args.std_path,
        int(args.kmer),
        args.log2,
        args.sortmethod,
        int(args.topkmernumber),
        int(args.xlabelsize),
        int(args.ylabelsize),
        int(args.xticksize),
        int(args.yticksize),
        int(args.lengendsize),
        args.outputname,
        args.pformat,
        int(args.pdpi),
        )



def console_kmer_msd_barplot():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_MSD_BARPLOT_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", 
                        help="full path of the input fasta sequence file with unique headers.")
    parser.add_argument("mean_path", help="full path to the normalization mean vector numpy file.")
    parser.add_argument("std_path",  help="full path to the normalization std vector numpy file.")
    parser.add_argument("kmer", help="length of kmers you want to count, which should be consistent with the normalization mean and std vector numpy files.")
    parser.add_argument(
        "-l",
        "--log2",
        default='Log2.post',
        choices=['Log2.post', 'Log2.pre', 'Log2.none'],
        help="decided if and when to log transform counts",
        )
    parser.add_argument(
        "-ss", 
        "--sortstat", 
        default='mean', 
        choices=['mean', 'sd'],
        help="whether to sort the kmer words based on mean or sd."
        )
    parser.add_argument(
        "-sm", 
        "--sortmethod", 
        default='descending', 
        choices=['ascending', 'descending'],
        help="how to sort the kmer words based on selected --sortstat."
        )
    parser.add_argument("-tn", "--topkmernumber", default=10, help="the number of sorted kmer words to plot, if input number is more than total number of kmer words, plot all words.")
    parser.add_argument("-xls", "--xlabelsize", default=20, help="x axis label font size.")
    parser.add_argument("-yls", "--ylabelsize", default=20, help="y axis label font size.")
    parser.add_argument("-xts", "--xticksize", default=20, help="x tick label font size.")
    parser.add_argument("-yts", "--yticksize", default=20, help="y tick label font size.")
    parser.add_argument(
        "-o", 
        "--outputname", 
        default="test", 
        help=(
            "the name and path of your output barplot."
            "a trailing part '_kmer_msd_barplot' is then automatically added to the name."
            )
            )
    parser.add_argument("-pf", "--pformat", default='pdf', help="the format of the output barplot, other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'.")
    parser.add_argument("-d", "--pdpi", default=300, help="the dpi of the output barplot.")
    args = _parse_args_or_exit(parser)

    kmer_msd_barplot.kmer_msd_barplot(
        args.fasta,
        args.mean_path,
        args.std_path,
        int(args.kmer),
        args.log2,
        args.sortstat,
        args.sortmethod,
        int(args.topkmernumber),
        int(args.xlabelsize),
        int(args.ylabelsize),
        int(args.xticksize),
        int(args.yticksize),
        args.outputname,
        args.pformat,
        int(args.pdpi),
        )



def console_kmer_textplot():
    assert sys.version_info[0] == 3, "Python version must be 3.x"
    parser = argparse.ArgumentParser(usage=KMER_TEXTPLOT_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("seq1file", 
                        help="full path of the first fasta file, only the first sequence will be plotted.")
    parser.add_argument("seq2file", 
                        help="full path of the second fasta file, only the first sequence will be plotted.")
    parser.add_argument("words_str", help=(
        "words of interest, max 10 words in the format of a comma separated string, example: 'ATTA,AAAA,ACTC'." 
        "if inputed more than 10 , only the first 10 will be plotted."
        )
        )
    parser.add_argument(
        "-cv",
        "--color_vec_str",
        default='default',
        help=(
            "a string of comma separated hex color for the words of interest, for example: '#d62728,#e377c2,#ff7f0e'. "
            "if 'default', will use the 'tab10' pallette with rearrangements into a quasi-rainbow order."
            "please make sure the number of hex color is the same as the number of the words of interest."
        )
        )
    parser.add_argument("-wl", "--wraplen", default=60, help="how many words to wrap in each line.")
    parser.add_argument("-cs", "--char_spacing", default=0.1, help="space between characters in the plot.")
    parser.add_argument("-ls", "--line_spacing", default=0.2, help="line space between seq1, seq2 and number.")
    parser.add_argument("-sfs", "--seqfontsize", default=42, help="sequence character font size.")
    parser.add_argument("-nfs", "--numfontsize", default=40, help="sequence position number font size.")
    parser.add_argument(
        "-o", 
        "--outputname", 
        default="test", 
        help=(
            "the name and path of your output textplot."
            "a trailing part '_kmer_textplot' is then automatically added to the name."
            )
            )
    parser.add_argument("-pf", "--plotformat", default='pdf', help="the format of the output textplot, other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'.")
    parser.add_argument("-d", "--plotdpi", default=300, help="the dpi of the output textplot.")
    args = _parse_args_or_exit(parser)

    # Split the input hex color string into a list using ',' as the delimiter
    words = args.words_str.split(',')
    
    if (args.color_vec_str == 'default'):
        color_vec='default'
    else:
        color_vec = args.color_vec_str.split(',')

    kmer_textplot.kmer_textplot(
        args.seq1file,
        args.seq2file,
        words,
        color_vec,
        int(args.wraplen),
        float(args.char_spacing),
        float(args.line_spacing),
        int(args.seqfontsize),
        int(args.numfontsize),
        args.outputname,
        args.plotformat,
        int(args.plotdpi),
        )
    

