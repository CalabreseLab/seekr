# SEEKR

[![Build Status](https://travis-ci.com/CalabreseLab/seekr.svg?branch=master)](https://travis-ci.com/CalabreseLab/seekr)
[![Build Status](https://img.shields.io/pypi/v/seekr.svg)](https://pypi.python.org/pypi/seekr)

Find communities of nucleotide sequences based on kmer frequencies.

A web portal is available at [seekr.org](http://seekr.org).

`The v1.0.1 update contains additional log transformation options and a tweak for length normalization`

## Updates
Includes a redesigned flag to indicate the method of k-mer standardization, and an additional option for *k*-mer standardization: --log2 [1,2,3] or -l [1,2,3]

1.	Default standardizationmethod is --log2 2. This is the same default standardization method used in SEEKR v1.0.0. For a given set of sequences, k-mers are counted, then length normalized (counts per kb of sequence), then z-scores for each k-mer are calculated, and then these z-scores are log2-tranformed. See PMID 31097619 for examples and an in-depth description of the rationale for using log2-transformed z-scores as a default.

2.	--log2 3 is the same as the no-log transformation option (-nl) of seekr v1.0.0. Here, after k-mers are counted and length-normalized, z-scores are calculated and used without any transformation. This was the approached used in our original SEEKR publication ([PMID 30224646](https://pubmed.ncbi.nlm.nih.gov/30224646/).

3. --log2 1 is the additional/new option for standardization. In this approach, k-mers are counted across a set of sequences and then length normalized (counts per kb of sequence). For each k-mer that has a zero-count value in each sequence, a pseudo-count of 1 is added; this allows the k-mer count values to be log2 transformed. Z-scores are calculated after log2-transformation of k-mer counts, and these z-scores can then be used directly for comparisons. In effect, this standardization method is not much different from the default option of --log2 2 (users can compare for themselves). It is, however, a slightly cleaner heuristic. In the time since our original two publications using seekr, we have noted that k-mer counts in the mouse and human transcriptomes tend to follow a log-normal distribution; thus, we currently favor this method of standardization. 

4. Includes small updates to “notes” and “Help” sections of README.md

5. Includes a small fix to the length normalization step of the core k-mer counting code. In v1.0.0, k-mer counts were normalized to the number of basepairs in each input sequence. In v1.0.1, k-mer counts are normalized to the number of k-mers in each input sequence ([length_of_sequence] -[k-mer_length]+1). This is a minor change for correctness that should not meaningfully affect results.


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
For example, you can use `seekr_kmer_counts` to generate a normalized kmer count matrix of `m` rows by `n` columns,
where `m` is the number of transcripts in a fasta file and `n` is 4^kmer.
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
 
Some general advice for thinking about how to use SEEKR. One challenge that we continually face in the lab is there are few ground truths in the lncRNA field and thus it is often unclear how to decide on the best parameters for sequence comparisons using SEEKR. Below are two points that may be useful – these are also discussed in the conclusions of [PMID 31097619](https://pubmed.ncbi.nlm.nih.gov/31097619/):

*	On the selection of a set of sequences to use for the calculation of standardization vectors: In our experience, one of the most useful features of SEEKR is that it provides a metric of relative similarity. “At the level of k-mers, sequence X is more similar to sequence Y than it is similar to 99% of other mouse lncRNAs”. Here, k-mer counts in sequence X and Y are being compared relative to the background k-mer frequencies found in all mouse lncRNAs. The calculated similarity of sequence X and Y derives from the standardization vectors used to calculate z-scores; in this example, those standardization vectors are the mean and standard deviation of the k-mer counts found all mouse lncRNAs. If users are trying to use SEEKR to determine how similar are two sets of sequences to each other, the key question to ask is, “how similar relative to what”? What set of sequences would be the most appropriate reference for the comparison? In many cases, using a large set of reference sequences to calculate standardization vectors is ideal (e.g. all mouse lncRNAs). However, by changing the set of reference sequences, users can change the question being asked. For example, if a scientist wanted to determine the similarity between sequences X and Y relative to sequence Z, then rather than all mouse lncRNAs, sequence Z should be used as a reference sequence. The reference sequence(s) can be controlled by using “seekr_norm_vectors” prior to “seekr_kmer_counts”.

*	On the selection of k-mer length: In our experience, the most robust biological trends have been relatively insensitive to the length of k-mer used in SEEKR. Still, when deciding on a length of k to use for comparisons, we recommend using a k-mer length for which 4^k is similar to the length of the average feature or key feature that is being compared. The reason for this is that as the length of k increases, so does the number of zero values for k-mer counts in a given sequence. For example, there are 16384 possible 7-mers. If users were interested in finding lncRNAs that are similar to lncRNA-X, which is 500 nucleotides long, a k-mer length of 7 would not be ideal, because the vector of 7-mer counts that corresponds to lncRNA-X would be dominated by zero values. In this example, unless users had a specific rationale for searching 7-mers, a k-mer length of 4 (256 possible k-mers) or 5 (1024 possible k-mers) would provide the basis for a stronger comparison.


### Help

Please also see a pre-print to a methods paper we wrote last year. This paper was originally scheduled to appear in Methods in Molecular Biology in 2020 but its publication date may be delayed. 

If you have questions about how you can use seekr in your own research, please send an email to jmcalabr@med.unc.edu

For full documentation, type...

```
$ seekr_kmer_counts
```

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
To reduce kmer redundancy due to these isoforms,
we can filter for transcripts ending in "01",
as indicated by the sequence headers:

```
$ seekr_canonical_gencode v22_lncRNAs.fa v22-01.fa
```

#### seekr_kmer_counts

Let's make a small `.csv` file of counts.
We'll set a couple flags:
* `--kmer 2` so we only have 16 kmers
* `--outfile out_counts.csv`.
This file will contain the log2-transformed z-scores of kmer counts per kb.

```
$ seekr_kmer_counts example.fa -o out_counts.csv -k 2
$ cat out_counts.csv
```

You can also see the output of this command
[here](https://github.com/CalabreseLab/seekr/seekr/tests/data/example_2mers.csv).

Three options are available for log transformation, using the --log2 flag. Pass `--log2 1` for log transformation of length normalized *k*-mer counts, with a +1 pseudo-count, pass `--log2 2` for log transformation of z-scores following count standardization (default), and pass `--log2 3` for no log transformation.

If we want to avoid normalization, we can produce kmer counts per kb by setting the `--log2 3`, `--uncentered` and `--unstandardized` flags:

```
$ seekr_kmer_counts example.fa -o out_counts.csv -k 2 --log2 3 -uc -us
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
This is likely due to a kmer not appearing in any of your sequences. Try:
1) using a smaller kmer size,
2) beginning with a larger set of sequences,
3) passing precomputed normalization vectors from a larger data set (e.g. GENCODE).

```

The code runs, but we get a warning.
That's because we're normalizing 4096 columns of kmers.
Most of those kmers never appear in any of our 5 lncRNAs.
This necessarily results in division by 0.
If we use a much larger set of [sequences](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.lncRNA_transcripts.fa.gz),
this same line works fine:

```
$ kmer_counts gencode.fa -o gencode_counts.npy
```

But what should you do if you're only interested in specific sequences?

#### seekr_norm_vectors

An effective way to find important kmers in a small number of RNAs is to
count their kmers, but normalize their counts to mean and
standard deviation vectors produced from a larger set of transcripts.
We can produce these vectors once, then use them on multiple smaller sets
of RNAs of interest. To produce the vectors, run:

**Note** If --log2 1 is passed in *seekr_kmer_counts*, then the -cl flag must be passed to *seekr_norm_vectors*. This is so that the log-transformed *k*-mer counts are standardized against reference *k*-mer counts that are also log transformed. 

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

To find Pearson correlations between kmer count profiles, run `seekr_pearson`.
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
One way to think about "random" with respect to kmers and RNA sequences,
is to think about what would happen if the nucleotide or kmer contents was conserved, but shuffled.
`seekr_gen_rand_rnas` provides a way to conserve but shuffle kmers.

To conserve dinucleotide content, for example, run:

```
$ seekr_gen_rand_rnas example.fa example_rand.fa -k 2 -s 0
```

The `--kmer` flag sets the size of the kmer,
and the `--seed` flag makes sure you can reproduce the resulting sequences.
We could now repeat our experiment with the new `example_rand.fa` file.

### Module example

For larger or more specific workloads, it may be better to use the `seekr` module.
In this example, we'll calculate similarities between two example fasta files,
(e.g., XIST and a set of RNAs we think could be similar to XIST)
using the normalization vectors from the human GENCODE set.
We'll use all kmers from 3 to 6, and label transcripts with unique labels.

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

    # Count kmers
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
* `{k}mers_xist.npy`: Normalized kmer profile for Xist.
* `{k}mers_lncs.npy`: Normalized kmer profile for other lncRNAs of interest.
* `xist_vs_lncs_{k}mers.npy`: Pearson's r values for all pairwise comparisons between Xist and the other lncRNAs.
* `xist_vs_lncs_{k}mers.csv`: Labeled, plain text version of pairwise comparisons.

## Issues

Any suggestions, questions, or problems can be directed to our
[GitHub Issues page](https://github.com/CalabreseLab/seekr/issues).

## Citation

If you use this work, please cite:

```
Kirk, J. M., Kim, S. O., Inoue, K., Smola, M. J., Lee, D. M., Schertzer, M. D., … Calabrese, J. M. (2018). Functional classification of long non-coding RNAs by k -mer content. Nature Genetics, 50(10), 1474–1482. https://doi.org/10.1038/s41588-018-0207-8
```
