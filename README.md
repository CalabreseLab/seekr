# SEEKR

[![Build Status](https://travis-ci.com/CalabreseLab/seekr.svg?branch=master)](https://travis-ci.com/CalabreseLab/seekr)
[![Build Status](https://img.shields.io/pypi/v/seekr.svg)](https://pypi.python.org/pypi/seekr)

Find communities of nucleotide sequences based on *k*-mer frequencies.

## Installation

 To use this library, you have to have >Python3.6 on your computer.

 Once you have Python, run:

 ```
 $ pip install seekr
 ```

 which will make both the command line tool and the python module available.

### CentOS

Users have been successful in installing `seekr` from source on CentOS:

```
conda create --name seekr_source python=3.8
conda activate seekr_source
git clone https://github.com/CalabreseLab/seekr.git
python3 setup.py install
conda install python-igraph
conda install louvain
```

See [this issue](https://github.com/CalabreseLab/seekr/issues/10) for further discussion.

## Usage

You can either use SEEKR from the command line or as a python module.
The package is broken up into a set of tools, each of which perform a single task.
From the command line, all of the functions will begin with `seekr_`.
For example, you can use `seekr_kmer_counts` to generate a normalized *k*-mer count matrix of `m` rows by `n` columns,
where `m` is the number of transcripts in a fasta file and `n` is 4^*k*-mer.
Then  `seekr_pearson` can be used to calculate how well correlated all pairwise combinations of sequences are.

To see all tools and some examples, run:

```
$ seekr
```

### Quickstart

To get a .csv file of communities for every transcript in a small .fa file called  [`example.fa`](https://raw.githubusercontent.com/CalabreseLab/seekr/master/seekr/tests/data/example.fa),
(where RNAs have been normalized to a data set of canonical transcripts from [GENCODE](https://www.gencodegenes.org/),
we would run:

```
$ seekr_download_gencode lncRNA
$ seekr_canonical_gencode v29_lncRNA.fa v29-01.fa # Name may change with GENCODE updates.
$ seekr_norm_vectors v29-01.fa
$ seekr_kmer_counts example.fa -o 6mers.csv -mv mean.npy -sv std.npy
$ seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv
$ seekr_graph example_vs_self.csv 0.13 -g example_vs_self.gml -c communities.csv
$ cat example_vs_self.csv
```

This quickstart procedure produces a number of other files beside the file `communities.csv`.
See below to learn about the other files produced along the way.

**Notes:**

* Some advanced usages are not available from the command line and require that you import the module.
* We'll use [`example.fa`](https://raw.githubusercontent.com/CalabreseLab/seekr/master/seekr/tests/data/example.fa)
as a small sample set,
if you want to download that file and follow along.
* GENCODE is a high quality source for human and mouse lncRNA annotation.
Fasta files can be found [here](https://www.gencodegenes.org/releases/current.html).
  * In the examples below we'll generically refer to `gencode.fa`.
    Any sufficiently large fasta file can be used, as needed.

Here are some examples if you just want to get going.

### Command line examples

#### seekr_download

Browsing [GENCODE](https://www.gencodegenes.org/) is nice if you want to explore fasta file options.
But if you know what you want, you can just download it from the command line.
This tool is also helpful on remote clusters.

To download all human transcripts of the latest release into a fasta file, run:

```
    $ seekr_download_gencode all
```

GENCODE also stores mouse sequences. You can select mouse using the `--species` flag:

```
    $ seekr_download_gencode all -s mouse
```

For consistency across experiments, you may want to stick to a particular release of GENCODE. To get lncRNAs from the M5 release of mouse, use `--release`:

```
    $ seekr_download_gencode lncRNA -s mouse -r M5
```

Finally, if you do not want the script to automatically unzip the file, you can leave the fasta file gzipped with `--zip`:

```
    $ seekr_download_gencode all -z
```

#### seekr_canonical_gencode

GENCODE fasta files provide multiple transcripts per genomic loci.
To reduce *k*-mer redundancy due to these isoforms,
we can filter for transcripts ending in "01",
as indicated by the sequence headers:

```
$ seekr_canonical_gencode v22_lncRNAs.fa v22-01.fa
```

#### seekr_kmer_counts

Let's make a small `.csv` file of counts.
We'll set a couple flags:
* `--kmer 2` so we only have 16 *k*-mers
* `--outfile out_counts.csv`.
This file will contain the log2-transformed z-scores of *k*-mer counts per kb.

```
$ seekr_kmer_counts example.fa -o out_counts.csv -k 2
$ cat out_counts.csv
```

You can also see the output of this command
[here](https://github.com/CalabreseLab/seekr/blob/master/seekr/tests/data/example_2mers.csv).

Three options are available for log transformation, using the `--log2` flag.
Pass `--log2 pre` for log transformation of length normalized *k*-mer counts, with a +1 pseudo-count,
pass `--log2 post` for log transformation of z-scores following count standardization (this is the default),
and pass `--log2 none` for no log transformation.

If we want to avoid normalization, we can produce *k*-mer counts per kb by setting the `--log2 non`, `--uncentered` and `--unstandardized` flags:

```
$ seekr_kmer_counts example.fa -o out_counts.csv -k 2 --log2 none -uc -us
```

Similarly, if we want a more compact, efficient numpy file,
we can add the `--binary` and `--remove_label` flags:

```
$ seekr_kmer_counts example.fa -o out_counts.npy -k 2 --binary --remove_label
```

**Note:** This numpy file is binary, so you won't be able to view it directly.

What happens if we also remove the `--kmer 2` option?

```
$ seekr_kmer_counts example.fa -o out_counts.npy
~/seekr/seekr/kmer_counts.py:143: RuntimeWarning: invalid value encountered in true_divide
  self.counts /= self.std

WARNING: You have `np.nan` values in your counts after standardization.
This is likely due to a *k*-mer not appearing in any of your sequences. Try:
1) using a smaller *k*-mer size,
2) beginning with a larger set of sequences,
3) passing precomputed normalization vectors from a larger data set (e.g. GENCODE).

```

The code runs, but we get a warning.
That's because we're normalizing 4096 columns of *k*-mers.
Most of those *k*-mers never appear in any of our 5 lncRNAs.
This necessarily results in division by 0.
If we use a much larger set of [sequences](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.lncRNA_transcripts.fa.gz),
this same line works fine:

```
$ kmer_counts gencode.fa -o gencode_counts.npy
```

But what should you do if you're only interested in specific sequences?

#### seekr_norm_vectors

An effective way to find important *k*-mers in a small number of RNAs is to
count their *k*-mers, but normalize their counts to mean and
standard deviation vectors produced from a larger set of transcripts.
We can produce these vectors once, then use them on multiple smaller sets
of RNAs of interest. To produce the vectors, run:

**Note:** If `--log2 post` is passed in *seekr_kmer_counts*, then the `--log2 post` flag must be passed to *seekr_norm_vectors*. 
This is so that the log-transformed *k*-mer counts are standardized against reference *k*-mer counts that are also log transformed.

```
$ seekr_norm_vectors gencode.fa
```

If you run `ls`, you should see `mean.npy` and `std.npy` in your directory.

To specify the path of these output files,
use the `--mean_vector` and `--std_vector` flags:

```
$ seekr_norm_vectors gencode.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy
```

Now, we can use these vectors to analyze our RNAs of interest:

```
$ kmer_counts example.fa -o out_5mers_gencode_norm.csv -k 5 -mv mean_5mers.npy -sv std_5mers.npy
```

#### seekr_pearson

To find Pearson correlations between *k*-mer count profiles, run `seekr_pearson`.
Running the program and options are similar to `seekr_kmers_counts`.
Input files for `seekr_pearson` will always be the output files from
one or more runs of `kmer_counts`.
The default settings accept two csv files and output a third csv file.

```
$ seekr_pearson out_counts.csv out_counts.csv -o example_vs_self.csv
$ cat example_vs_self.csv
```

The only other options besides the `-o` flag control binary versus plain text input and output.
If you have a binary input file (i.e. a .npy file),
and also want a binary output file, you can use the `--binary_input` and `--binary_output` flags:

```
$ seekr_pearson out_counts.npy out_counts.npy -o example_vs_self.npy -bi -bo
```

If we want to compare counts between two files
(e.g. RNAs between mouse and human),
that is also possible:

```
$ seekr_pearson human_6mers.npy mouse_6mers.npy -o human_vs_mouse.npy
```

#### seekr_graph

We can treat the results of `seekr_pearson` as an [adjacency matrix](https://en.wikipedia.org/wiki/Adjacency_matrix), and use it to find communities of transcripts with the [Louvain algorithm](https://arxiv.org/abs/0803.0476) or the [Leiden algorithm](https://arxiv.org/abs/1810.08473).

The default setting accept a csv file and a threshold value.
This csv file should be the product of seekr_pearson, or some other adjacency matrix.
The threshold is the value below which edges are removed from the graph.
A [gml](https://gephi.org/users/supported-graph-formats/gml-format/) file contain the graph and communities will be produced.

```
$ seekr_graph example_vs_self.csv .13 -g example.gml
```

Numpy files are also valid input:

```
$ seekr_graph example_vs_self.npy .13 -g graph.gml
```

GML files are plain text, so you can view them if you want.
But they contain all of the information describing the graph.
Often you just want to know what transcripts belong to what community.
To get a csv file mapping transcripts to communities, run:

```
$ seekr_graph example_vs_self.csv .13 -g example.gml -c communities.csv
```

The value for thresholding the adjacency matrix is very experiment dependent.
Because a reasonable default is difficult to predict, it is a required parameter.
For example, if you are using `k=5`, it likely makes sense to increase the threshold (e.g. to .3).
`seekr_pearson_distro` can be run to suggest a value for the threshold.
Similarly, the community finding algorithm allows you course control of community size with its resolution parameter (gamma).
Testing ranges from .1 to 5 is reasonable, where values from 1-3 have been most useful in our experience.

```
    $ seekr_graph example_vs_self.csv .3 -g graph.gml -r 1.5
```

A third tunable parameters is a cap of the number of communities found.
And finally, since community finding is partially random,
you can make your results reproducible by setting a seed value:

```
    $ seekr_graph example_vs_self.csv .13 -g graph.gml -n 10 -s 0
```

#### seekr_gen_rand_rnas

It's often useful to understand what we might expect "at random".
One way to think about "random" with respect to *k*-mers and RNA sequences,
is to think about what would happen if the nucleotide or *k*-mer contents was conserved, but shuffled.
`seekr_gen_rand_rnas` provides a way to conserve but shuffle *k*-mers.

To conserve dinucleotide content, for example, run:

```
$ seekr_gen_rand_rnas example.fa example_rand.fa -k 2 -s 0
```

The `--kmer` flag sets the size of the *k*-mer,
and the `--seed` flag makes sure you can reproduce the resulting sequences.
We could now repeat our experiment with the new `example_rand.fa` file.

### Module example

For larger or more specific workloads, it may be better to use the `seekr` module.
In this example, we'll calculate similarities between two example fasta files,
(e.g., XIST and a set of RNAs we think could be similar to XIST)
using the normalization vectors from the human GENCODE set.
We'll use all *k*-mers from 3 to 6, and label transcripts with unique labels.

```python
import numpy as np
import pandas as pd
from seekr.kmer_counts import BasicCounter
from seekr.pearson import pearson
from seekr.fasta_reader import Reader


gencode = 'gencode.fa'
xist = 'xist.fa'
lncRNAs = 'other_lncs.fa'

# Make sure each lncRNA in other_lncs.fa has a unique name
headers = Reader(lncRNAs).get_headers()
names = [h.strip('>') + f'_{i}' for i, h in enumerate(headers)]

for k in range(3, 7):
    # Make normalization vectors
    gencode_counter = BasicCounter(gencode, k=k)
    gencode_counter.get_counts()
    mean_path = f'mean_{k}mers.npy'
    std_path = f'std_{k}mers.npy'
    np.save(mean_path, gencode_counter.mean)
    np.save(std_path, gencode_counter.std)

    # Count *k*-mers
    xist_counter = BasicCounter(xist,
                                outfile=f'{k}mers_xist.npy',
                                mean=mean_path,
                                std=std_path,
                                k=k)
    lncs_counter = BasicCounter(lncRNAs,
                                outfile=f'{k}mers_lncs.npy',
                                mean=mean_path,
                                std=std_path,
                                k=k)
    xist_counter.make_count_file(names=['XIST'])
    lncs_counter.make_count_file(names=names)

    # Find similarities
    sim = pearson(xist_counter.counts,
                  lncs_counter.counts,
                  outfile=f'xist_vs_lncs_{k}mers.npy')

    # Save labeled .csv file of similarities.
    sim_df = pd.DataFrame(sim, ['XIST'], names)
    sim_df.to_csv(f'xist_vs_lncs_{k}mers.csv')

```

Each loop will write six files to disk:

* `mean_{k}mers.npy`: Mean vector for GENCODE human lncRNAs.
Once this has been saved, the first portion of the code doesn't need to be run again.
* `std_{k}mers.npy`: Standard deviation vector for GENCODE human lncRNAs.
* `{k}mers_xist.npy`: Normalized *k*-mer profile for Xist.
* `{k}mers_lncs.npy`: Normalized *k*-mer profile for other lncRNAs of interest.
* `xist_vs_lncs_{k}mers.npy`: Pearson's r values for all pairwise comparisons between Xist and the other lncRNAs.
* `xist_vs_lncs_{k}mers.csv`: Labeled, plain text version of pairwise comparisons.

## Seekr2.0.0 functions update

In seekr2.0.0 we updated 8 functions spanning from stastics to visualizations. The statistics part focuses on fitting a background distribution and then using this background distribution to calculate p values for pearson's correlation r values. In this way, p value could be used straightforwardly to call pairs of significant similarities. The visualization part focuses on plotting the r or p value from various aspects as well as rendering sequence similarities based on seekr mode. 

### Statistics

#### find_dist
Find the best fitted model to the input background sequences. Data to be fitted are all possible pairwise kmer pearson correlation for the background sequences. User can define the list of models for fitting. The output is usually a fitted model, but could also be the actual data (pearson correlation r values as an numpy array). This is helpful if the background distribution is small, in which case fitting to any model would not give a descent result. The results (best fitted model or actual data) will be used to calculate p-values.

If fit_model is True, it returns a dataframe contains the distributions name, goodness of fit (like D stats in ks) and fitted params for all models.
If fit_model is False, it returns the actual data (with or without subsetting based on subsetting argument) as an numpy array, without fitting to any distributions.
If plotfit is given, a plot with the given name and path will be saved, showing the fitted distributions (red dash line) vs the actual data (blue solid line) for all distributions in models.
If outputname is given, a csv file (outputname_fitres.csv) will be saved for either a list of fitted models or an numpy array.

Python example:
```python
from seekr import find_dist

fitres = find_dist.find_dist(inputseq='default', k_mer=4, models='common10', 
                             subsetting=True, subset_size = 10000, 
                             fit_model=True, statsmethod='ks',progress_bar=True, 
                             plotfit='testmodefitplot.pdf', outputname='test')
```

Console example (minimal and full versions):
```
$ seekr_find_dist default

$ seekr_find_dist default -k 4 -l Log2.pre -mdl common10 -sbt -sbs 100000 -fm -statm ks -pb -pf modelfit.pdf -o test
```

#### find_pval
Calculte p values of the seekr.pearson correlation r values for input sequences 1 vs input sequences 2, sequences 1 and sequences 2 can be the same. The calculation is based on the output of find_dist, which is either a list of models or an numpy array. This function connects the output of find_dist and the input sequnces of interests. Given the background sequencs, find_dist calculates all possible pairwise seekr.pearson values and then outputs either a list of fitted models or the numpy array of the actual data. find_pval firstly calculates the seekr.pearson of the two input sequences, which produces a correlation matrix. p value is then calculated for each r value in the pearson correlation matrix based on the output of find_dist.

The output of find_pval is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences. If save_file is True, the dataframe is saved to a csv file named as outputname_pval.csv. 

Python example:
```python
from seekr import find_pval

pvals=find_pval.find_pval(seq1file='test1.fa', seq2file='test2.fa', 
                          mean_path='mean_4mers.npy', std_path='std_4mers.npy',
                          k_mer=4, fitres=fitres, log2='Log2.post', 
                          bestfit=1, outputname='test', progress_bar=True)
```

Console example (minimal and full versions):
```
$ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv

$ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv -ft npy -bf 1 -o test -pb
```

### Visualization

#### kmer_heatmap
Customizeable heatmap (outputname_kmer_heatmap.hformat) for easier visualization of the results of seekr.pearson (r values) or find_pval (p values). It takes in a dataframe with row and column names and plot the heatmap and dendrograms for both rows and columns. kmer_dendrgram is a good alternative to plot only the dendrograms (partial or full) to get a better idea of the clustering. For interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot and kmer_textplot for further analysis.

Python example:
```python
from seekr import kmer_heatmap

kmer_heatmap.kmer_heatmap(df, datamin=-1, datamax=1, thresh_value=0.00,
                          color_range=['#1b7837', '#ffffff', '#c51b7d'],
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.4, hmaph_ratio=0.3, x_tick_size=10, 
                          y_tick_size=10, cbar_font_size=26, 
                          outputname='test',hformat='pdf', hdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_heatmap pval.csv 0 1

$ seekr_kmer_heatmap pval.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o /Users/username/Desktop/test -hf pdf -hd 300
```

#### kmer_dendrogram
Customizeable dendrograms (outputname_kmer_dendrogram_row/column.p format) for easier visualization of the hierarchical clustering results of seekr.pearson (r values) or find_pval (p values). Performs hierarchial clustering on either rows or columns of a dataframe, which is a better way to visualize the whole or partial clusters in kmer_heatmap. For interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot and kmer_textplot for further analysis.

Python example:
```python
from seekr import kmer_dendrogram

kmer_dendrogram.kmer_dendrogram(df, dendro_direct='row', 
                                distmetric='correlation', linkmethod='complete', 
                                plot_ht=8, wd_ratio=0.5, leaf_font_size = 16, 
                                outputname='test', pformat='pdf', pdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_dendrogram pval.csv

$ seekr_kmer_dendrogram pval.csv -dd column -distm correlation -linkm complete -ph 8 -hratio 0.5 -lfs 16 -o test -pf pdf -d 300
```

#### kmer_leiden
Plot Leiden community network (outputname_gradient/threshold_leiden_network.pdf ) for input fasta seqeunces: calculte sequences distance matrix as seekr pearson correlation of the fasta file to itself, then use Leiden algorithms to call cluster. The resulting undirected network edges are weighted with the pearson correlation r values: darker and thicker edges have higher weights. It works better with small number of nodes (less than 50). For better visualization and customization, please set savecsv to True and import the saved nodes and edges files (outputname_edges_leiden.csv and outputname_nodes_leiden.csv) directly to Gephi.

Python example:
```python
from seekr import kmer_leiden

kmer_leiden.kmer_leiden(inputfile='testld.fa',mean='mean_4mers.npy',std='std_4mers.npy', 
                        k=4, algo='RBERVertexPartition', rs=1.0,
                        edgecolormethod='gradient', edgethreshold=0.1, labelfontsize=12,
                        outputname='test', savecsv=True)
```

Console example (minimal and full versions):
```
$ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -s

$ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -a RBERVertexPartition -r 1.0 -ec threshold -et 0.1 -m 12 -o test -s
```

#### kmer_count_barplot
Barplot (outputname_kmer_count_barplot.pformat) of the transformed or raw z-score for kmer words of the input sequences (limit to 10). Calculate the z-score with BasicCounter, where user can define whether and how to do log2 transform. Order kmer words by the summed difference from the mean among all sequences, in descending or ascending order. Plot the top x kmer words with bars differently colored for each input sequences.

Python example:
```python
from seekr import kmer_count_barplot

kmer_count_barplot.kmer_count_barplot(inputfile='test.fa', mean='mean_4mers.npy', 
                                       std = 'std_4mers.npy', log2 = 'Log2.post', 
                                       k=4, sortmethod='ascending', topkmernumber=10,
                                       xlabelsize=20, ylabelsize=20, 
                                       xticksize=20, yticksize=20, lengendsize=12,
                                       outputname='test', pformat='pdf', pdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4

$ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4 -l Log2.post -sm ascending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o test -pf pdf -d 300
```

#### kmer_msd_barplot
Barplot (outputname_kmer_msd_barplot.pformat) of the mean of the transformed or raw z-score for kmer words across input sequences. Error bar represents the standard deviation. Order kmer words by mean or sd, in descending or ascending order. Plot the top x kmer words with each bar corresponds to one kmer word.

Python example:
```python
from seekr import kmer_msd_barplot

kmer_msd_barplot.kmer_msd_barplot(inputfile='test.fa', mean='mean_4mers.npy',
                                  std = 'std_4mers.npy', log2 = 'Log2.post',
                                  k=4, sortstat='sd', sortmethod='ascending', 
                                  topkmernumber=10,xlabelsize=20, ylabelsize=20, 
                                  xticksize=20, yticksize=20, 
                                  outputname='test', pformat='pdf', pdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4

$ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o test -pf pdf -d 300
```

#### kmer_textplot
Highlight the interested kmer words in the input sequences (outputname_kmer_textplot.plotformat). Align the 2 input sequences at the beginning and label sequence positions on the bottom. Highlight in colors the words of interest. For multiple words, if the words overlap, the overlapped characters will be highlighted in the color of the first word in the list. Therefore arrange the words of interest based on the priority of the words, with the most important word at the beginning of the list.

Python example:
```python
from seekr import kmer_textplot

kmer_textplot.kmer_textplot(seq1file='test.fa', seq2file='test2.fa', 
                            words=['CTCT','GTAG','AAAA','GCGC'],
                            color_vec='default', wraplen=60, 
                            char_spacing=0.1, line_spacing=0.3,
                            seqfontsize=72, numfontsize=40, 
                            outputname='test', plotformat='pdf', plotdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC'

$ seekr_kmer_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC' -cv '#d62728,#e377c2,#ff7f0e' -wl 60 -cs 0.1 -ls 0.2 -sfs 42 -nfs 40 -o test -pf pdf -d 300
```


## Additional considerations

This section is a "not-so-quickstart", providing more complete views on selection of input data and parameter selection.

Some general advice for thinking about how to use SEEKR.
One challenge that we continually face in the lab is there are few ground truths in the lncRNA field and thus it is often unclear how to decide on the best parameters for sequence comparisons using SEEKR.
Below are some points that may be useful – these are also discussed in the conclusions of [PMID 31097619](https://pubmed.ncbi.nlm.nih.gov/31097619/)

### Selection of a set of sequences to use for the calculation of standardization vectors

In our experience, one of the most useful features of SEEKR is that it provides a metric of relative similarity.
Consider two lncRNAs, lncRNA-X, which is from mouse, and lncRNA-Y, which is from human.
One way to compare these two lncRNAs using SEEKR would be to calculate their *k*-mer profiles and compare these profiles in relation to all known mouse lncRNAs.
To do this, one would first use all mouse lncRNAs as an input for the "seekr_norm_vectors" script, to determine the mean and standard deviation of counts for all *k*-mers in all mouse lncRNAs.
Then, users would take those standardization vectors along with the sequences of lncRNA-X and lncRNA-Y and use them in the “kmer_counts” script to calculate *k*-mer profiles of lncRNA-X and lncRNA-Y in relation to all mouse lncRNAs.
Finally, users would employ the “seekr_pearson” script to determine how similar lncRNA-X and lncRNA-Y were to each other relative to all other mouse lncRNAs.
The point here is that the set of sequences used to calculate standardization vectors is a key variable – perhaps analogous to a reference gene in a quantitative PCR experiment or a loading control in a western blot.
Users should think through what comparison they are interested in performing and choose their set of standardization sequences accordingly.
In the example above, where all mouse lncRNAs were used for standardization, users might discover that “At the level of *k*-mers, lncRNA-X is more similar to lncRNA-Y than it is similar to 99% of other mouse lncRNAs”.
But changing the set of sequences for standardization can change the question being asked.
For example, users might be interested in comparing lncRNA-X and lncRNA-Y relative to all known human enhancer RNAs (eRNAs).
Here, using all human eRNAs as the standardization set, users might make an additional insightful discovery, perhaps: “At the level of *k*-mers, lncRNA-X is no more similar to lncRNA-Y than it is similar to the average human eRNA”.
Relatedly, when users are comparing two large groups of sequences (let us call these “set A” and “set B”), we again recommend thinking about what set of reference sequences would be best for standardization.
In most cases, users will probably want to use the same set of reference sequences to standardize *k*-mer counts in set A and in set B.  
Perhaps set A and set B should again be standardized relative to all mouse lncRNAs; or, both set A and set B be should be standardized relative to the sequences in set B.
But standardizing set A relative to itself, then standardizing set B relative to itself, then using “seekr_pearson” to compare the two sets of sequences might yield a non-sensical comparison.

### Selection of *k*-mer length

In our experience, the most robust biological trends have been relatively insensitive to the length of *k*-mer used in SEEKR.
Still, when deciding on a length of k to use for comparisons, we recommend using a *k*-mer length for which 4^k is similar to the length of the average feature or key feature that is being compared.
The reason for this is that as the length of k increases, so does the number of zero values for *k*-mer counts in a given sequence.
For example, there are 16384 possible 7-mers.
If users were interested in finding lncRNAs that are similar to lncRNA-X, which is 500 nucleotides long, a *k*-mer length of 7 would not be ideal, because the vector of 7-mer counts that corresponds to lncRNA-X would be dominated by zero values.
In this example, unless users had a specific rationale for searching 7-mers, a *k*-mer length of 4 (256 possible *k*-mers) or 5 (1024 possible *k*-mers) would provide the basis for a stronger comparison.

## Issues and Help

If you have questions about how you can use seekr in your own research, please send an email to jmcalabr@med.unc.edu

For full documentation, run:

```
$ seekr
```

Any suggestions, questions, or problems can be directed to our
[GitHub Issues page](https://github.com/CalabreseLab/seekr/issues).

Please also see [the pre-print](https://github.com/CalabreseLab/seekr/blob/logchanges/methods_mol_bio_seekr-v20.pdf) to a methods paper we wrote last year. This paper was originally scheduled to appear in Methods in Molecular Biology in 2020 but its publication date may be delayed.

## Citation

If you use this work, please cite:

```
Kirk, J. M., Kim, S. O., Inoue, K., Smola, M. J., Lee, D. M., Schertzer, M. D., … Calabrese, J. M. (2018). Functional classification of long non-coding RNAs by k -mer content. Nature Genetics, 50(10), 1474–1482. https://doi.org/10.1038/s41588-018-0207-8
```
