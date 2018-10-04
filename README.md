# SEEKR

A library for counting small kmer frequencies in nucleotide sequences.

A web portal is available at [seekr.org](http://seekr.org).

## Installation

 To use this library, you have to have Python3.x on your computer. 
 If you don't have it installed, the easiest place to get it is from the 
 [Anaconda distribution](https://www.continuum.io/downloads).

 Once you have Python, run:

 ```commandline
 $ pip install seekr
 ```

 which will make both the command line tool and the python module available.

## Usage

You can either use SEEKR from the command line or as a python module. 
In either case, you will use `kmer_counts` to generate a kmer count matrix of `m` rows by `n` columns,
where `m` is the number of transcripts in a fasta file and `n` is 4^kmer. 
Then  `pearson` can be used to calculate how well correlated all pairwise combinations of sequences are.

**Notes:** 

* Some advanced usages are not available from the command line and require that you import the module.
* We'll use [`example.fa`](https://github.com/CalabreseLab/seekr/seekr/tests/data/example.fa) 
as a small sample set,
if you want to open that file and follow along.
* GENCODE is a high quality source for human and mouse lncRNA annotation.
Fasta files can be found [here](https://www.gencodegenes.org/releases/current.html).
  * In the examples below we'll generically refer to `gencode.fa`.
    Any sufficiently large fasta file can be used, as needed.

Here are some quick-start examples if you just want to get going.

### Command line examples

#### kmer_counts

Let's make a small `.csv` file we can view.
We'll set several flags:
* `--nonbinary` so the output is plain text
* `--kmer 2` so we only have 16 kmers
* `--label` so there are column and row labels

```commandline
$ kmer_counts example.fa -o out_counts.csv -k 2 -nb -lb
$ cat out_counts.csv

```

You can also see the output of this command 
[here](https://github.com/CalabreseLab/seekr/seekr/tests/data/example_2mers.csv).


If we want a more compact, efficient numpy file,
we can drop the `--nonbinary` and `--label` flags:

```commandline
$ kmer_counts example.fa -o out_counts.npy -k 2
```

**Note:** This numpy file is binary, so you won't be able to view it directly.

What happens if we also remove the `--kmer 2` option?

```commandline
$ kmer_counts example.fa -o out_counts.npy
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

```commandline
$ kmer_counts gencode.fa -o gencode_counts.npy
```

But what should you do if you're only interested in specific sequences?

#### norm_vectors

An effective way to find important kmers in a small number of RNAs is to
count their kmers, but normalize their counts to mean and 
standard deviation vectors produced from a larger set of transcripts.
We can produce these vectors once, then use them on multiple smaller sets
of RNAs of interest. To produce the vectors, run:

```commandline
$ norm_vectors gencode.fa 
```

If you run `ls`, you should see `mean.npy` and `std.npy` in your directory.

To specify the path of these output files,
use the `--mean_vector` and `--std_vector` flags:

```commandline
$ norm_vectors gencode.fa -k 7 -mv mean_7mers.npy -sv std_7mers.npy
```

Now, we can use these vectors to analyze our RNAs of interest:

```commandline
$ kmer_counts example.fa -o out_7mers_gencode_norm.npy -k 7 -mv mean_7mers.npy -sv std_7mers.npy
```

#### pearson

To find Pearson correlations between kmer count profiles, run `pearson`. 
Running the program and options are similar to `kmers_counts`. 
Input files for `pearson` will always be the output files from 
one or more runs of `kmer_counts`. 
The default setting accept two numpy files and output a third numpy file.

```commandline
$ pearson out_counts.npy out_counts.npy -o example_vs_self.npy
```

The only other options besides the `-o` flag control binary versus .csv input and output. 
If you have a non binary input file (i.e. a .csv file), 
and also want a non binary output file, you can do:

```commandline
$ pearson out_counts.csv out_counts.csv -o example_vs_self.csv -nbi -nbo
$ cat example_vs_example.csv

```

If we want to compare counts between two files 
(e.g. RNAs between mouse and human), 
that is also possible:

```commandline
$ pearson human_6mers.npy mouse_6mers.npy -o human_vs_mouse.npy
```

#### Summary

If we want to get a .csv file that has all pairwise comparisons of `example.fa`,
where RNAs have been normalized to `gencode.fa` using 6mers, we would run:

```commandline
$ norm_vectors gencode.fa
$ kmer_counts example.fa -o 6mers.npy -mv mean.npy -sv std.npy
$ pearson 6mers.npy 6mers.npy -o example_vs_self.csv -nbo
$ cat example_vs_self.csv
```
### Module example

For larger or more specific workloads, it may be better to use the `seekr` module.
In this example, we'll calculate similarities between two example fasta files,
(e.g., XIST and a set of RNAs we think could be similar to XIST)
using the normalization vectors from the human GENCODE set.
We'll use all kmers from 3 to 7, and label transcripts with unique labels.

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

for k in range(3, 8):
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

### Help

For full documentation of the parameters and flags, you can run `kmer_counts`  or  `pearson` without any arguments.

```
$ kmer_counts
```

## Issues

Any suggestions, questions, or problems can be directed to our 
[GitHub Issues page](https://github.com/CalabreseLab/seekr/issues).

## Citation

If you use this work, please cite:

```
Kirk, J. M., Kim, S. O., Inoue, K., Smola, M. J., Lee, D. M., Schertzer, M. D., … Calabrese, J. M. (2018). Functional classification of long non-coding RNAs by k -mer content. Nature Genetics, 50(10), 1474–1482. https://doi.org/10.1038/s41588-018-0207-8
```
