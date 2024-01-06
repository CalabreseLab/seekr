# SEEKR

[![Build Status](https://travis-ci.com/CalabreseLab/seekr.svg?branch=master)](https://travis-ci.com/CalabreseLab/seekr)
[![Build Status](https://img.shields.io/pypi/v/seekr.svg)](https://pypi.python.org/pypi/seekr)

Find communities of nucleotide sequences based on *k*-mer frequencies.

## Installation

 To use this library, you have to have >=Python3.9 on your computer.

 Once you have Python, run:

 ```
 $ pip install seekr
 ```

 which will make both the command line tool and the python module available.

### CentOS

Users have been successful in installing `seekr` from source on CentOS:

```
conda create --name seekr_source python=3.9
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
$ seekr_download_gencode lncRNA -g
$ seekr_filter_gencode v33_lncRNA.fa -gtf v33_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -o v33 # Name may change with GENCODE updates.
$ seekr_norm_vectors v33_filtered.fa
$ seekr_kmer_counts example.fa -o 6mers.csv -mv mean.npy -sv std.npy
$ seekr_pearson 6mers.csv 6mers.csv -o example_vs_self.csv
$ cat example_vs_self.csv
```

This quickstart procedure produces a number of other files beside the file `example_vs_self.csv`.
See below to learn about the other files produced along the way.

**Notes:**

* We'll use [`example.fa`](https://raw.githubusercontent.com/CalabreseLab/seekr/master/seekr/tests/data/example.fa)
as a small sample set,
if you want to download that file and follow along.
* GENCODE is a high quality source for human and mouse lncRNA annotation.
Fasta files can be found [here](https://ftp.ebi.ac.uk/pub/databases/gencode/).
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

If you do not want the script to automatically unzip the file, you can leave the fasta file gzipped with `--zip`:

```
    $ seekr_download_gencode all -z
```

Finally, if you want to download the gtf file from the same species and same release for further useage in sequence filtering with seekr_filter_gencode, you can do so with `--gtf`:
```
    $ seekr_download_gencode all -g
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
Pass `--log2 Log2.pre` for log transformation of length normalized *k*-mer counts, with a +1 pseudo-count,
pass `--log2 Log2.post` for log transformation of z-scores following count standardization (this is the default),
and pass `--log2 Log2.none` for no log transformation.

If we want to avoid normalization, we can produce *k*-mer counts per kb by setting the `--log2 Log2.none`, `--uncentered` and `--unstandardized` flags:

```
$ seekr_kmer_counts example.fa -o out_counts.csv -k 2 --log2 Log2.none -uc -us
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

**Note:** If `--log2 Log2.post` is passed in *seekr_kmer_counts*, then the `--log2 Log2.post` flag must be passed to *seekr_norm_vectors*. 
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

### Module example

For larger or more specific workloads, it may be better to use the `seekr` module in python.
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

In seekr2.0.0 we updated 9 functions spanning from sequence processing, stastics to visualizations. The statistics part focuses on fitting a background distribution and then using this background distribution to calculate p values for pearson's correlation r values. In this way, p value could be used straightforwardly to call pairs of significant similarities. The visualization part focuses on plotting the r or p value from various aspects as well as rendering sequence similarities based on seekr mode. 

### Sequence processing

#### filter_gencode
Filter gencode fasta file by 1: length; 2: transcript feature type and Ensemble_canonical tag 3: transcript feature type and isoform number 4: remove duplicated sequences. Please use gencode downloaded fasta and gtf files as input, other formats are not tested. Make sure the fasta and gtf files are from the same release and same species. If canonical filtering is turned on, only sequences that has the feature type 'transcript' and has a tag 'Ensembl_canonical' will be kept. For older version of gencode gtfs that does not have the tag 'Ensembl_canonical', please use isoform filter. If isoform filtering is turned on, only sequences that has the feature type 'transcript' and has a transcript_name ('Gm37180-201') that contains the isoform number (201) will be kept. Isoform number support regular expression, for example, `'[0-9]01'` will match 101, 201, 301 up to 901. Please **quote the argument** if using regular expression in the command line to avoid errors. If rm_dup is set to True, sequences that are exactly the same will be removed and only the first occurence will be kept. If more than 50 transcript_ids in gtf file cannot be matched to the fasta headers, there will be a warining. transcript_id in gtf file is located inside the 9th field. transcript_id in the fasta headers is the first element of the header split by `|`. transcript_name is also located inside the 9th field. These are standard formats from gencode. Please refer to [Gencode data format website](https://www.gencodegenes.org/pages/data_format.html) for further details. 


Python example:
```python
from seekr import filter_gencode

headers, seqs = filter_gencode.filter_gencode(fasta_path='gencode.vM33.lncRNA_transcripts.fa', 
                                              gtf_path='gencode.vM33.long_noncoding_RNAs.gtf',
                                              len_threshold=500, canonical=True, isoform='201',
                                              rm_dup=False, outputname='test')
```

Console examples:

filter by length only (>=500)
```
$ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -len 500 -o test
```

filter by length (>=500) and Ensemble_canonical tag and isoform number 201
```
$ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -gtf gencode.vM33.long_noncoding_RNAs.gtf -len 500 -can -iso 201 -o test
```

filter by Ensemble_canonical tag and remove duplicated sequences
```
$ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -gtf gencode.vM33.long_noncoding_RNAs.gtf -can -rd -o test
```

filter by isoform number with regular expression, remeber to **quote the full argument** when using regular expression
```
$ seekr_filter_gencode gencode.vM33.lncRNA_transcripts.fa -gtf gencode.vM33.long_noncoding_RNAs.gtf -iso '[0-9]01' -o test
```


### Statistics

#### find_dist
Find the best fitted model to the input background sequences. Data to be fitted are all possible pairwise kmer pearson correlation for the background sequences. User can define the list of models for fitting. The output is usually a fitted model, but could also be the actual data (pearson correlation r values as an numpy array). This is helpful if the background distribution is small, in which case fitting to any model would not give a descent result. The results (best fitted model or actual data) will be used to calculate p-values.

If fit_model is True, it returns a dataframe contains the distributions name, goodness of fit (like D stats in ks) and fitted params for all models.
If fit_model is False, it returns the actual data (with or without subsetting based on subsetting argument) as an numpy array, without fitting to any distributions.
If plotfit is given, a pdf plot with the given name and path will be saved, showing the fitted distributions (red dash line) vs the actual data (blue solid line) for all distributions in models.
If outputname is given, a csv file (outputname.csv) will be saved for either a list of fitted models or an numpy array.

Python example:
```python
from seekr import find_dist

fitres = find_dist.find_dist(inputseq='default', k_mer=4, models='common10', 
                             subsetting=True, subset_size = 10000, 
                             fit_model=True, statsmethod='ks',progress_bar=True, 
                             plotfit='modefitplot', outputname='test_fitres')
```

Console example (minimal and full versions):
```
$ seekr_find_dist default -o test_fitres

$ seekr_find_dist default -k 4 -l Log2.pre -mdl common10 -sbt -sbs 100000 -fm -statm ks -pb -pf modelfitplot -o test_fitres
```

#### find_pval
Calculate p values of the seekr.pearson correlation r values for input sequences 1 vs input sequences 2, sequences 1 and sequences 2 can be the same. The calculation is based on the output of find_dist, which is either a list of models or an numpy array. This function connects the output of find_dist and the input sequnces of interests. Given the background sequencs, find_dist calculates all possible pairwise seekr.pearson values and then outputs either a list of fitted models or the numpy array of the actual data. find_pval firstly calculates the seekr.pearson of the two input sequences, which produces a correlation matrix. p value is then calculated for each r value in the pearson correlation matrix based on the output of find_dist.

The output of find_pval is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences. If outputname is given, the dataframe is saved to a csv file named as outputname.csv. 

Python example:
```python
from seekr import find_pval

pvals=find_pval.find_pval(seq1file='test1.fa', seq2file='test2.fa', 
                          mean_path='mean_4mers.npy', std_path='std_4mers.npy',
                          k_mer=4, fitres=fitres, log2='Log2.post', 
                          bestfit=1, outputname='test_pval', progress_bar=True)
```

Console example (minimal and full versions):
```
$ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv -o test_pval

$ seekr_find_pval seqs1.fa seqs2.fa mean_4mer.npy std_4mer.npy 4 fitres.csv -ft npy -bf 1 -o test_pval -pb
```

#### adj_pval
Adjust p-values for multiple comparisons. This function performs multiple comparison correction for p values calculated by find_pval. Multiple comparison correction is performed by statsmodels.stats.multitest.multipletests function. The input is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences. The output is a dataframe of adjusted p values, with the same dimension as the input dataframe. If the input dataframe is a symmetric matrix: the row and column (or input1 and input2 in find_pval) are the same. The number of multiple comparisons is automatically corrected to be half the matrix size, excluding the diagonal. That is only the upper triangle of the matrix (excluding diagonal) is used for multiple comparison correction. In the output dataframe, the upper triangle of the matrix (excluding diagonal) is filled with the adjusted p values, while the lower triangle and the diagonal is filled with NaN. If the input dataframe is not a symmetric matrix the total matrix is used for multiple comparison correction

The output of adj_pval is a dataframe of adjusted p values, with the row names (input 1) and column names (input 2) same as the input dataframe. If outputname is given, the dataframe is saved to a csv file named as outputname.csv. 

Python example:
```python
from seekr import adj_pval

adjpvals=adj_pval.adj_pval(pvals, method='bonferroni', alpha=0.05, outputname='test_bonferroni_0.05_adj_pval')
```

Console example (minimal and full versions):
```
$ seekr_adj_pval test_pval.csv bonferroni

$ seekr_adj_pval test_pval.csv bonferroni -a 0.05 -o test_bonferroni_0.05_adj_pval
```


### Visualization

#### kmer_heatmap
Customizeable heatmap (outputname.hformat) for easier visualization of the results of seekr.pearson (r values) or find_pval (p values). It takes in a dataframe with row and column names and plot the heatmap with/without dendrograms for both rows and columns. kmer_dendrgram is a good alternative to plot only the dendrograms (partial or full) to get a better idea of the clustering. For interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot, kmer_comp_textplot and kmer_indi_textplot for further analysis.

Python example:
```python
from seekr import kmer_heatmap

kmer_heatmap.kmer_heatmap(df, datamin=-1, datamax=1, thresh_value=0.00,
                          color_range=['#1b7837', '#ffffff', '#c51b7d'], 
                          cluster=True,
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.4, hmaph_ratio=0.3, x_tick_size=10, 
                          y_tick_size=10, cbar_font_size=26, 
                          outputname='kmer_heatmap',hformat='pdf', hdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_heatmap pval.csv 0 1

$ seekr_kmer_heatmap pval.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -cl -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o /Users/username/Desktop/kmer_heatmap -hf pdf -hd 300
```

#### kmer_dendrogram
Customizeable dendrograms (outputname.pformat) for easier visualization of the hierarchical clustering results of seekr.pearson (r values) or find_pval (p values). Performs hierarchial clustering on either rows or columns of a dataframe, which is a better way to visualize the whole or partial clusters in kmer_heatmap. For interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot, kmer_comp_textplot and kmer_indi_textplot for further analysis.

Python example:
```python
from seekr import kmer_dendrogram

kmer_dendrogram.kmer_dendrogram(df, dendro_direct='row', 
                                distmetric='correlation', linkmethod='complete', 
                                plot_ht=8, wd_ratio=0.5, leaf_font_size = 16, 
                                outputname='kmer_dendrogram', pformat='pdf', pdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_dendrogram pval.csv

$ seekr_kmer_dendrogram pval.csv -dd column -distm correlation -linkm complete -ph 8 -wratio 0.5 -lfs 16 -o kmer_dendrogram -pf pdf -d 300
```

#### kmer_leiden
Plot Leiden community network (plotname.pdf ) for input fasta seqeunces: calculte sequences distance matrix as seekr pearson correlation of the fasta file to itself. After filtering the edges with pearsoncutoff, use Leiden algorithms to call cluster. Seed can be set to get reporducible community assignments. The resulting undirected network edges are weighted with the pearson correlation r values: darker and thicker edges have higher weights. It works better with small number of nodes (less than 50). For better visualization and customization, please set csvfile to a customized name (string) and import the saved nodes and edges files (csvfile_edges_leiden.csv and csvfile_nodes_leiden.csv) directly to Gephi.

Python example:
```python
from seekr import kmer_leiden

kmer_leiden.kmer_leiden(inputfile='testld.fa',mean='mean_4mers.npy',std='std_4mers.npy', 
                        k=4, algo='RBERVertexPartition', rs=1.0, pearsoncutoff=0, setseed=False,
                        edgecolormethod='gradient', edgethreshold=0.1, labelfontsize=12,
                        plotname='kmer_leiden', csvfile=None)
```

Console example (minimal and full versions):
```
$ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -cf testdata

$ seekr_kmer_leiden ldseq.fa mean_4mer.npy std_4mer.npy 4 -a RBERVertexPartition -r 1.0 -pco 0.1 -sd -ec threshold -et 0.2 -m 12 -pn kmer_leiden -cf testdata
```

#### kmer_count_barplot
Barplot (outputname.pformat) of the transformed or raw z-score for kmer words of the input sequences (limit to 10). Calculate the z-score with BasicCounter, where user can define whether and how to do log2 transform. Order kmer words by the summed difference from the mean among all sequences, in descending or ascending order. Plot the top x kmer words with bars differently colored for each input sequences.

Python example:
```python
from seekr import kmer_count_barplot

kmer_count_barplot.kmer_count_barplot(inputfile='test.fa', mean='mean_4mers.npy', 
                                       std = 'std_4mers.npy', log2 = 'Log2.post', 
                                       k=4, sortmethod='ascending', topkmernumber=10,
                                       xlabelsize=20, ylabelsize=20, 
                                       xticksize=20, yticksize=20, lengendsize=12,
                                       outputname='kmer_count_barplot', pformat='pdf', pdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4

$ seekr_kmer_count_barplot test.fa mean_4mer.npy std_4mer.npy 4 -l Log2.post -sm ascending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o kmer_count_barplot -pf pdf -d 300
```

#### kmer_msd_barplot
Barplot (outputname.pformat) of the mean of the transformed or raw z-score for kmer words across input sequences. Error bar represents the standard deviation. Order kmer words by mean or sd, in descending or ascending order. Plot the top x kmer words with each bar corresponds to one kmer word.

Python example:
```python
from seekr import kmer_msd_barplot

kmer_msd_barplot.kmer_msd_barplot(inputfile='test.fa', mean='mean_4mers.npy',
                                  std = 'std_4mers.npy', log2 = 'Log2.post',
                                  k=4, sortstat='sd', sortmethod='ascending', 
                                  topkmernumber=10,xlabelsize=20, ylabelsize=20, 
                                  xticksize=20, yticksize=20, 
                                  outputname='kmer_msd_barplot', pformat='pdf', pdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4

$ seekr_kmer_msd_barplot test.fa mean_4mer.npy std_4mer.npy 4 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o kmer_msd_barplot -pf pdf -d 300
```

#### kmer_comp_textplot
Highlight the interested kmer words and compare them in two input sequences (outputname.plotformat). Align the 2 input sequences at the beginning and label sequence positions on the bottom. Highlight in colors the words of interest. For multiple words, if the words overlap, the overlapped characters will be highlighted in the color of the first word in the list. Therefore arrange the words of interest based on the priority of the words, with the most important word at the beginning of the list.

Python example:
```python
from seekr import kmer_comp_textplot

kmer_comp_textplot.kmer_comp_textplot(seq1file='test1.fa', seq2file='test2.fa', 
                                      words=['CTCT','GTAG','AAAA','GCGC'],
                                      color_vec='default', wraplen=60, 
                                      char_spacing=0.1, line_spacing=0.2,
                                      seqfontsize=72, numfontsize=40, colorblockh=0.15, 
                                      outputname='comp_textplot', plotformat='pdf', plotdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_comp_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC'

$ seekr_kmer_comp_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC' -cv '#d62728,#e377c2,#ff7f0e' -wl 60 -cs 0.1 -ls 0.2 -sfs 72 -nfs 40 -cbh 0.15 -o comp_textplot -pf pdf -d 300
```

#### kmer_indi_textplot
Highlight the input kmer words across all input sequences and save for each sequence a separate plot. This is a similar function as kmer_comp_textplot. But instead of aligning and comparing two sequences, kmer_indi_textplot highlight the interested words in all input sequences and plot each sequence individually. For multiple words, if the words overlap, the overlapped characters will be highlighted in the color of the first word in the list. Therefore arrange the words of interest based on the priority of the words, with the most important word at the beginning of the list.

Python example:
```python
from seekr import kmer_indi_textplot

kmer_indi_textplot.kmer_indi_textplot(seqfile='test.fa', 
                                      words=['CTCT','GTAG','AAAA','GCGC'],
                                      color_vec='default', wraplen=60, 
                                      char_spacing=0.1, line_spacing=0.3,
                                      seqfontsize=72, numfontsize=40, colorblockh=0.3, 
                                      outputpath='', plotformat='pdf', plotdpi=300)
```

Console example (minimal and full versions):
```
$ seekr_kmer_indi_textplot seq.fa 'ATTA,AAAA,ACTC'

$ seekr_kmer_indi_textplot seq.fa 'ATTA,AAAA,ACTC' -cv '#d62728,#e377c2,#ff7f0e' -wl 60 -cs 0.1 -ls 0.3 -sfs 72 -nfs 40 -cbh 0.3 -op /Users/username/Desktop/indi_textplot/ -pf pdf -d 300
```

## Seekr Docker Image
We now provide a Docker image of seekr which hopefully avoids all the package dependencies and complications when install and run seekr through pip. 

### Docker Installation
Firstly, you should install Docker on your computer. Please refer to the [official website](https://www.docker.com/get-started/) for installation instructions.

After you have successfully installed Docker. Start/Run the application and make sure it is running properly – you should see the docker icon on your task bar with the status indicated.

### Pull Docker Image and Test Run
1. Start your command line tool: Terminal for MacOS and CMD for Windows. You can also use Powershell or Cygwin for Windows, but Cygwin might have interaction issues.

From the command line, pull the Docker Image:
```
docker pull calabreselab/seekr:latest
```
You can replace `latest` with a specific tag if needed.

2. Test Run the Docker Image
```
docker run -it --rm calabreselab/seekr:latest
```
The `-it` tag sets it to interactive mode. If you don't need to run the Docker container in interactive mode (i.e., you don't need a shell inside the container), you can simply omit the `-it` flag.
This will print the user manual out to the command line, which is basically the same as you run the command `seekr` directly from command line when you pip install seekr. 

### Run Docker Image from command line
You can run the seekr function from this Docker Image directly from command line with the following specified syntax.
```
docker run -v /path/to/your/files:/data calabreselab/seekr:latest seekr_kmer_comp_textplot /data/seq1.fa /data/seq2.fa 'ATTA,AAAA,ACTC' -o /data/comp_textplot
```
In this command:
* `-v /path/to/your/files:/data`: This mounts the directory `/path/to/your/files` from your host machine (where seq1.fa and seq2.fa are located) to the `/data` directory inside the Docker container. Replace `/path/to/your/files` with the actual path to your files.
* `seekr_kmer_comp_textplot /data/seq1.fa /data/seq2.fa 'ATTA,AAAA,ACTC' -o /data/test`: This is the command that gets executed inside the Docker container. Since we mounted our files to `/data` in the container, we reference them with `/data/seq1.fa` and `/data/seq2.fa`. 
* The `/data` folder is basically a mirror of the folder you specified in `/path/to/your/files`. So by specifying `-o /data/comp_textplot` (output with the path /data/ and plotname as in kmer_textplot.pdf) we can have the output files directly written in to the folder in `/path/to/your/files`.
* Please remember to **specify your output path to `/data/`** otherwise it will not be saved to your folder on local machine and it would be hard to locate it even inside the Docker Container Filesystem (in this case, when the Docker Container is removed, your data will be deleted as well). 

Examples of code mounts e:/test on Windows as the folder that contains the input and holds the output files:
```
docker run -v e:/test:/data calabreselab/seekr:latest seekr_kmer_comp_textplot /data/test1.fa /data/test2.fa 'ATTA,AAAA,ACTC' -o /data/comp_textplot
```
```
docker run -v e:/test:/data calabreselab/seekr:latest seekr_kmer_leiden /data/testld.fa /data/bkg_mean_4mers.npy /data/bkg_std_4mers.npy 4 -s -o /data/kmer_leiden
```
Basically you need to add: `docker run -v /path/to/your/files:/data calabreselab/seekr:latest` before the command line code for seekr (see above for examples of all functions).

### Run Docker Image with Jupyter Notebook
If you want to work with python codes instead of directly calling from the command line, you can choose to run seekr with Jupyter Notebook inside the Docker Container.

1.  Run Docker Container with Interactive Terminal:
```
docker run -it -p 8888:8888 -v /path/on/host:/path/in/container calabreselab/seekr:latest /bin/bash
```
This command will start the container and give you a bash terminal inside it. The `-p 8888:8888` flag maps the port *8888* from the container to your host machine so that you can access the Jupyter Notebook server.

`/path/on/host` is the path to a directory on your local machine where you want the notebooks and other data to be saved. `/path/in/container` is the directory inside the container where you want to access and save the files from Jupyter Notebook.

When you use Jupyter Notebook now and create or save notebooks, they will be stored at `/path/in/container` inside the Docker container. Due to the volume mount, these files will also be accessible and stored at `/path/on/host` on your host machine so that you can later access the code and files even when the container is removed. 

Example of code:
```
docker run -it -p 8888:8888 -v e:/test:/data calabreselab/seekr:latest /bin/bash
```

2. Manually start Jupyter Notebook. From the bash terminal inside the Docker container:
```
jupyter notebook --ip=0.0.0.0 --port=8888 --NotebookApp.token='' --NotebookApp.password='' --allow-root
```
* `--ip=0.0.0.0`: Allow connections from any IP address. This is important for accessing Jupyter from outside the container.
* `--port=8888`: Run Jupyter on this port. You can change this if needed, but remember to adjust the `-p` flag when starting the Docker container.
* `--allow-root`: Allow Jupyter to be run as the root user inside the Docker container.
* `--NotebookApp.token=''`: This disables the token-based security measure.
* `--NotebookApp.password=''`: This ensures that there's no password required to access the Jupyter server. 

Disabling the token and password for Jupyter Notebook reduces security. It's generally okay for local development, but you should avoid doing this in production or any publicly accessible server.

3. Access the Jupyter Notebook from your host machine's browser by entering this address:
<http://localhost:8888/>

4. Run Python functions as demonstrated above. It would be convenient if all input files can be copied over to the folder you have mounted: `/path/on/host`. Pay attention that when you specify your input and output file route, use the `/path/in/container` route, as that is a mirror to your local folder. 

Example of code:
```python
from seekr import kmer_comp_textplot

kmer_comp_textplot.kmer_comp_textplot(seq1file='/data/test1.fa', seq2file='/data/test2.fa', 
                                      words=['ATCG','GTAG','AAAA','GCGC'], 
                                      color_vec='default', wraplen=60,
                                      char_spacing=0.1, line_spacing=0.2,
                                      seqfontsize=72, numfontsize=40, colorblockh=0.15, 
                                      outputname='/data/kmer_textplot', plotformat='pdf',
                                      plotdpi=300)
```
Once you are done, you can click the **Shut Down** button under the **File** tab of Jupyter Notebook to shut down the instance or you can just click `Ctrl+C` twice from the command line to kill the process. 
Then you need to exit the Docker Container from the interative session by typing `exit` and hit enter.

### Cleanup (Optional)
* Clean Up Docker Container:
    + List all containers, including the stopped ones: `docker ps -a`
    + To remove a specific container: `docker rm CONTAINER_ID_OR_NAME`
    + To remove all stopped containers: `docker container prune`

* Clean Up Docker Image. If you want to remove the Docker image you used:
    + List all images: `docker images`
    + Remove a specific image: `docker rmi IMAGE_ID_OR_NAME`

You'd typically only remove an image if you're sure you won't be using it again soon, or if you want to fetch a fresh version of it from a repository like Docker Hub.

* Additional Cleanup. Docker also maintains a cache of intermediate images and volumes. Over time, these can accumulate. To free up space:
    + Remove unused data: `docker system prune`
    + To also remove unused volumes (be careful, as this might remove data you want to keep): `docker system prune --volumes`

Remember to always be cautious when cleaning up, especially with commands that remove data. Make sure you have backups of any essential data, and always double-check what you're deleting.


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
