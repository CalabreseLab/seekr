# SEEKR

A library for counting small kmer frequencies in nucleotide sequences.

## Installation

 * To use this library, you have to have Python3.x on your computer. If you don't have it installed, the easiest place to get it is from the [Anaconda distribution](https://www.continuum.io/downloads). Downloading Anaconda will also provide you with most of the dependencies you need to use SEEKR.
 * Either download this repository as a .zip file, or use git to `clone` it.
 * Install any missing dependencies by running: `pip install -r requirements.txt`.

## Usage

You can either use SEEKR from the command-line or as a python module. In either case, you will use `kmer_counts.py` to generate a kmer count matrix of m rows by n columns,
where m is the number of transcripts in a fasta file and n is 4^kmer. Then  `pearson.py` can be used to calculate how well correlated all pairwise combinations of sequences are. ** Note: ** Some advanced usages are currently not available from the command-line and require that you import the the module.

Here are some quick-start examples if you just want to get going:

### Examples

#### kmer_counts

The default settings produce a binary, normalized numpy file:

```   
$ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.npy
```

To get a human readable csv file, set the nonbinary flag:

```
$ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.csv -nb
```

If you want to add default labels, also set the label flag:

```
$ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.csv -nb -lb
```

You can change also change the size of the kmer you're using, and prevent normalization:

```
$ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.npy -k 4 -nc -ns
```

#### pearson

To find Pearson correlations between kmer count profiles, run `pearson.py`. Running the program and options are similar to `kmers_counts.py`. Input files for `pearson.py` will always be the output files from one or more runs of `kmer_counts.py`. The default setting accept two numpy files and output a third numpy file.

```
$ python pearson.py /path/to/kc_out.npy /path/to/kc_out.npy -o /path/to/out.npy
```

The only other options besides the `-o` flag control binary versus .csv input and output. If you have a non binary input file (i.e. a .csv file) and also want a non binary output file, you can do:

```
$ python pearson.py /path/to/kc_out.npy /path/to/kc_out.npy -o /path/to/out.npy -nbi -nbo
```

### Advanced usage

A common task is to use the normalization vectors from a large .fa file to analyze specific lncRNAs of interest. Currently, this cannot be done from the command-line, but is still fairly straightforward from within python. This example will serve as a walkthrough on using SEEKR as a module and will take you through creating the normalization vectors (using 4mers) to creating the Pearson's similarity matrix.

```python
import numpy as np
from kmer_counts import BasicCounter
from pearson import pearson

v22_gencode = 'v22_gencode.fa'
gencode_counter = BasicCounter(v22_gencode, k=4)
gencode_counter.get_counts()
np.save('mean.npy', gencode_counter.mean)
np.save('std.npy', gencode_counter.std)

xist = 'xist.fa'
lncRNAs = 'other_lncs.fa'
xist_counter = BasicCounter(xist, '4mers_xist.npy', mean='mean.npy', std='std.npy', k=4)
lncs_counter = BasicCounter(lncRNAs, '4mers_lncs.npy', mean='mean.npy', std='std.npy', k=4)
xist_counter.make_count_file()
lncs_counter.make_count_file()

sim = pearson(xist_counter.counts, lncs_counter.counts, outfile='xist_vs_lncs.npy')

```

This will write five files to disk:

* `mean.npy`: Mean vector for gencode human lncRNAs. Once this has been saved, the first portion of the code doesn't need to be run again.
* `std.npy`: Standard deviation vector for gencode human lncRNAs.
* `4mers_xist.npy`: Normalized kmer profile for Xist.
* `5mers_lncs.npy`: Normalized kmer profile for other lncRNAs of interest.
* `xist_vs_lncs.npy`: Pearson's r values for all pairwise comparisons between Xist and the other lncRNAs.

### Help

For full documentation of the parameters and flags, you can run `kmer_counts.py`  or  `pearson.py` without any arguments.

```
$ python kmer_counts.py
```

## Issues

Any suggestions, questions, or problems can be directed to our [GitHub Issues page](https://github.com/CalabreseLab). #TODO update to full link after acceptance.
