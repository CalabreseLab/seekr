import sys
import argparse
import numpy as np
import pandas as pd

from seekr import fasta
from seekr import pearson
from seekr.kmer_counts import BasicCounter
from seekr import filter_gencode
from seekr import find_dist
from seekr import find_pval
from seekr import adj_pval
from seekr import kmer_heatmap
from seekr import kmer_dendrogram
from seekr import kmer_leiden
from seekr import kmer_count_barplot
from seekr import kmer_msd_barplot
from seekr import kmer_textplot

from seekr.__version__ import __version__


DOWNLOAD_GENCODE_DOC = """
Description
-----------
Download fasta and gtf files from https://www.gencodegenes.org/

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

If you want to also download the gtf file:
    $ seekr_download_gencode all -g

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

FILTER_GENCODE_DOC = """
Description:
filter gencode fasta file by 1: length; 2: transcript feature type and Ensemble_canonical tag 3: transcript feature type and isoform number

Details:
please use gencode downloaded fasta and gtf files as input, other formats are not tested
please make sure the fasta and gtf files are from the same release and same species
if only length filtering is needed, please set -can F and omit -gtf, omit -iso
if canonical filtering is turned on, only sequences that has the feature type 'transcript' and has a tag 'Ensembl_canonical' will be kept
if isoform filtering is turned on, only sequences that has the feature type 'transcript' and has a transcript_name ('Gm37180-201') that contains the isoform number ('201') will be kept
if more than 50 transcript_ids in gtf file cannot be matched to the fasta headers, the program will abort
transcript_id in gtf file is located inside the 9th field, and the transcript_id is the value of the key 'transcript_id'
transcript_id in the fasta headers is the first element of the header split by '|'
transcript_name in gtf file is also located inside the 9th field
these are standard formats from gencode


Examples:
filter by length (>=500) and Ensemble_canonical tag and isoform number '201'

    $ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -gtf gencode.vM33.long_noncoding_RNAs.gtf -len 500 -can -iso 201 -o test

filter by length only (>=500)

    $ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -len 500 -o test

filter by Ensemble_canonical tag only

    $ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -gtf gencode.vM33.long_noncoding_RNAs.gtf -can -o test

filter by isoform number only

    $ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -gtf gencode.vM33.long_noncoding_RNAs.gtf -iso 201 -o test

to save under other location change to -o /Users/username/Desktop/test

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

to not save the model fit plot, leave out -pf and its argument
to save under other location change to -o /Users/username/Desktop/test
to not save the model fit result, leave out -o and its argument

to use user defined models list -mdl 'norm,expon,pareto'

minimal code with all settings to default
fitting the models and save fitted results:
    $ seekr_find_dist default -fm -o test
not fitting the models, directly output the background numpy array:
    $ seekr_find_dist default -o test

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

to save under other location change to -o /Users/username/Desktop/test
to not save the result, leave out -o and its argument

minimal code with all settings to default

    $ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv -o test

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

ADJ_PVAL_DOC = """
Description: 
adjust p-values for multiple comparisons

Details:
this function performs multiple comparison correction for p values calculated by find_pval
multiple comparison correction is performed by statsmodels.stats.multitest.multipletests function
the input is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences
the output is a dataframe of adjusted p values, with the same dimension as the input dataframe
if the input dataframe is a symmetric matrix: the row and column (or input1 and input2 in find_pval) are the same
the number of multiple comparisons is automatically corrected to be half the matrix size, excluding the diagonal
that is only the upper triangle of the matrix (excluding diagonal) is used for multiple comparison correction
in the output dataframe, the upper triangle of the matrix (excluding diagonal) is filled with the adjusted p values
while the lower triangle and the diagonal is filled with NaN
if the input dataframe is not a symmetric matrix the total matrix is used for multiple comparison correction

Example:
use bonferroni correction with family-wise error rate at 0.05 to adjust the p values in test_pval.csv
which is the output of the find_pval function and save the output as test_bonferroni_alpha0.05_adjusted_pval.csv
    
    $ seekr_adj_pval test_pval.csv bonferroni -a 0.05 -o test

to save under other location change to -o /Users/username/Desktop/test
to not save the result, leave out -o and its argument

minimal code with all settings to default

    $ seekr_adj_pval test_pval.csv bonferroni

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

    $ seekr_kmer_heatmap pval.csv 0 1 -o test

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
    
    $ seekr_kmer_dendrogram pval.csv -dd row -distm correlation -linkm complete -ph 8 -wratio 0.5 -lfs 16 -o test -pf pdf -d 300

to save under other location change to -o /Users/username/Desktop/test

minimal code with all settings to default

    $ seekr_kmer_dendrogram pval.csv -o test

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/seekr/
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

"""

KMER_LEIDEN_DOC = """
Description:
plot Leiden community network for input fasta seqeunces
it works better with small number of nodes (less than 50)
for better visualization and customization, please set csvfile to a customized name (string) and use Gephi with the saved nodes and edges files

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
seed is set so that the community assignment results are reproducible each time the code is run
the plot edge will be gradient colored based on the pearson correlation values
the plot will be saved as test_gradient_leiden_network.pdf under current directory
the corresponding nodes and edges files will be saved as testdata_nodes_leiden.csv and testdata_edges_leiden.csv under current directory
    
    $ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -a RBERVertexPartition -r 1.0 -pco 0.1 -sd -ec gradient -et 0.1 -lfs 12 -pn test -cf testdata

to save plot under other location change to -pn /Users/username/Desktop/test
to suppress plotting, leave out -pn and its argument
to suppress saving nodes and edges files, leave out -cf and its argument
to turn set seed off to test robustness of community assignments, remove -sd from the command

minimal code with all settings to default

    $ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -cf testdata

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

    $ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4 -o test

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

    $ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4 -o test

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


def _run_download_gencode(biotype, species, gtf, release, fasta_path, gtf_path, unzip):
    # Note: This function is separated for testing purposes.
    downloader = fasta.Downloader()
    downloader.get_gencode(biotype, species, gtf, release, fasta_path, gtf_path, unzip)


def console_download_gencode():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=DOWNLOAD_GENCODE_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("biotype", help=("Name of Genocde set to download. " "Must be one of ('all', 'pc', 'lncRNA')."))
    parser.add_argument(
        "-s", "--species", default="human", help=" Name of species. Must be one of: ('human' or 'mouse')."
    )
    parser.add_argument(
        "-g", "--gtf", action="store_true", help="Set if you want to download the accompanying gtf file (gencode.vXX.chr_patch_hapl_scaff.annotation.gtf.gz) as well."
    )
    parser.add_argument(
        "-r",
        "--release",
        default=None,
        help=("Name of specific release to download (e.g. 'M5'). " "If None, download latest release."),
    )
    parser.add_argument(
        "-fp", "--fasta_path", default=None, help="Path to location for fasta file. Default will save by release name."
    )
    parser.add_argument(
        "-gp", "--gtf_path", default=None, help="Path to location for gtf file if --gtf is set. Default will save by release name."
    )
    parser.add_argument(
        "-z", "--zip", action="store_false", help="Set if you do not want to gunzip fasta file after downloading."
    )
    args = _parse_args_or_exit(parser)
    _run_download_gencode(args.biotype, args.species, args.gtf, args.release, args.fasta_path, args.gtf_path, args.zip)



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
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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


def _run_norm_vectors(fasta, mean_vector, std_vector, log2, kmer):
    counter = BasicCounter(fasta, k=int(kmer), log2=log2)
    counter.get_counts()
    np.save(mean_vector, counter.mean)
    np.save(std_vector, counter.std)


def console_norm_vectors():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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



def console_filter_gencode():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=FILTER_GENCODE_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", 
                        help="full path of the fasta file to be filtered, best downloaded from gencode.")
    parser.add_argument(
        "-gtf", 
        "--gtf_path", 
        default=None, 
        help=(
            "the full path to the gtf file that matched the fasta file, best downloaded from gencode."
            "only needed when filtering by Ensemble_canonical tag (-can)."
        ),
            )
    parser.add_argument(
        "-len", 
        "--len_threshold", 
        default=0, 
        help="length threshold to filter fasta sequences, only sequences longer than or equal to this threshold will be kept."
        )
    parser.add_argument(
        "-can", 
        "--canonical", 
        action="store_true", 
        help="whether to filter by Ensemble_canonical tag."
        )
    parser.add_argument(
        "-iso", 
        "--isoform", 
        default="0", 
        help="isoform number (ex: '201') to filter fasta sequences, default '0' will turn off the filter."
        )
    parser.add_argument("-o", "--outputname", default="test", help="the path and name to save the filtered fasta file, will automatically be combined with the trailing part _filtered.fa.")
    args = _parse_args_or_exit(parser)

    filter_gencode.filter_gencode(
        args.fasta,
        args.gtf_path,
        int(args.len_threshold),
        args.canonical,
        args.isoform,
        args.outputname,
        )


def console_find_dist():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    parser.add_argument("-o", "--outputname", default=None, help="path and name to save fitted results as a csv file.")
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
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    parser.add_argument("-o", "--outputname", default=None, help="the full path to store the p value csv file.")
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

def console_adj_pval():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=ADJ_PVAL_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pval_path", help="the full path to the .csv file that contains the p values to be adjusted, should be the output of find_pval.")
    parser.add_argument("method", help="the method used for multiple comparison correction, refer to statsmodels.stats.multitest.multipletests for all options.")
    parser.add_argument("-a", "--alpha", default=0.05, help="the desired family-wise error rate: the significance level for the entire set of tests.")
    parser.add_argument("-o", "--outputname", default=None, help="the full path to store the adjusted p value csv file.")
    args = _parse_args_or_exit(parser)


    # Read the CSV file into a DataFrame with the headers and index names
    pvals = pd.read_csv(args.pval_path, header=0, index_col=0)

    adj_pval.adj_pval(
        pvals,
        args.method,
        float(args.alpha),
        args.outputname,
        )

def console_kmer_heatmap():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    parser.add_argument("-wratio", "--wd_ratio", default=0.5, help="a ratio factor to control the dendrgram width.")
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
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    parser.add_argument(
        "-sd", 
        "--setseed", 
        action="store_true", 
        help="set to True to set seed, which gives reproducible results for community assignment."
            )
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
        "-pn", 
        "--plotname", 
        default=None, 
        help=(
            "the name and path of your output file."
            "this will broadcast to all saved images (with trailing part _gradient_leiden_network.pdf or _threshold_leiden_network.pdf)."
            )
            )
    parser.add_argument(
        "-cf", 
        "--csvfile", 
        default=None, 
        help=(
            "the name and path of the nodes and edges csv files." 
            "a trailing part _edges_leiden.csv and _nodes_leiden.csv will be added to the files."
            "the saved edges and nodes csv files are readily importable to Gephi for a better network visualization."
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
        args.setseed,
        args.edgecolormethod,
        float(args.edgethreshold),
        int(args.labelfontsize),
        args.plotname,
        args.csvfile,
        )



def console_kmer_count_barplot():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    parser.add_argument("-ls", "--legendsize", default=12, help="legend font size.")
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
        int(args.legendsize),
        args.outputname,
        args.pformat,
        int(args.pdpi),
        )



def console_kmer_msd_barplot():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
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
        "seekr_filter_gencode": FILTER_GENCODE_DOC,
        "seekr_norm_vectors": NORM_VECTORS_DOC,
        "seekr_kmer_counts": KMER_COUNTS_DOC,
        "seekr_pearson": PEARSON_DOC,
        "seekr_find_dist": FIND_DIST_DOC,
        "seekr_find_pval": FIND_PVAL_DOC,
        "seekr_adj_pval": ADJ_PVAL_DOC,
        "seekr_kmer_heatmap": KMER_HEATMAP_DOC,
        "seekr_kmer_dendrogram": KMER_DENDROGRAM_DOC,
        "seekr_kmer_leiden": KMER_LEIDEN_DOC,
        "seekr_kmer_count_barplot": KMER_COUNT_BARPLOT_DOC,
        "seekr_kmer_msd_barplot": KMER_MSD_BARPLOT_DOC,
        "seekr_kmer_textplot": KMER_TEXTPLOT_DOC,
    }
    for c, d in cmds2doc.items():
        print(f"{'='*25}\n{c}\n{'='*25}\n{d}")
    conclusion = (
        "To see a full description of flags and defaults, "
        "run any of the commands listed above, without any parameters "
        '(e.g. "$ seekr_kmer_leiden").'
    )
    print(conclusion)


def console_seekr_help():
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--version", action="store_true", help="Print current version and exit.")
    args = parser.parse_args()
    _run_console_seekr_help(args.version)

