# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

---

## [Current]

### Added

*

### Fixed

*

### Changed

*
---

## 1.5.1

### Fixed

* --log2 compatibility in domain pearson. ([See Issue](https://github.com/CalabreseLab/seekr/issues/19))
* Update examples. ([See Issue](https://github.com/CalabreseLab/seekr/issues/18))
  
## 1.5.0

### Added

* Unofficial support for different alphabets
* Additional log transformation option

### Fixed

* Length normalization divisor (length to [length -k +1])
  * Includes a small fix to the length normalization step of the core *k*-mer counting code. In v1.4.2, *k*-mer counts were normalized to the number of basepairs in each input sequence. In v1.4.3, *k*-mer counts are normalized to the number of *k*-mers in each input sequence ([length_of_sequence] -[*k*-mer_length]+1). This is a minor change for correctness that should not meaningfully affect results.


### Changed

* Behavior of --log2 argument
  * Includes a redesigned flag to indicate the method of *k*-mer standardization, and an additional option for *k*-mer standardization: --log2 [1,2,3] or -l [1,2,3]. These options correspond to log-transforming pre-standardization, post-standardization, or no log-transform,
respectively. 
  * `--log2 post` (Default standardization method). This is the same default standardization method used in SEEKR v1.4.2. For a given set of sequences, *k*-mers are counted, then length normalized (counts per kb of sequence), then z-scores for each *k*-mer are calculated, and then these z-scores are log2-tranformed. See [PMID 31097619](https://pubmed.ncbi.nlm.nih.gov/31097619/) for examples and an in-depth description of the rationale for using log2-transformed z-scores as a default.
  * `--log2 none` is the same as the no-log transformation option (-nl) of seekr v1.4.2. Here, after *k*-mers are counted and length-normalized, z-scores are calculated and used without any transformation. This was the approached used in our original SEEKR publication ([PMID 30224646](https://pubmed.ncbi.nlm.nih.gov/30224646/)).
  * `--log2 pre` is the additional/new option for standardization. In this approach, *k*-mers are counted across a set of sequences and then length normalized (counts per kb of sequence). For each *k*-mer that has a zero-count value in each sequence, a pseudo-count of 1 is added; this allows the *k*-mer count values to be log2 transformed. Z-scores are calculated after log2-transformation of *k*-mer counts, and these z-scores can then be used directly for comparisons. In effect, this standardization method is not much different from the default option of --log2 2 (users can compare for themselves). It is, however, a slightly cleaner heuristic. In the time since our original two publications using seekr, we have noted that *k*-mer counts in the mouse and human transcriptomes tend to follow a log-normal distribution; thus, we currently favor this method of standardization. 
* Includes small updates to “notes” and “Help” sections of README.md



---

## 1.4.0

### Added

* `seekr_pwms` is now callable from the command line
* `seekr_gen_rand_rnas` is live

### Fixed

* Updated README
* Let `seekr_visualize_distro` handle other matrices
### Changed

* In `seekr_domain_pearson`, change the way percentiles are calculated, to now be relative to a reference fasta.
* Improve error when passing a bad release to `seekr_download_gencode`.


---

## 1.3.0

### Added

* `seekr_canonical_gencode` command line script filters for -001 transcripts.
* Example integration script in the test directory.
* Add legacy option to continue using Louvain instead of Leiden.
* Travis CI automatic push testing.
* `seekr_visualize_distro` command makes distribution of r-values.
* `seekr_domain_pearson` command line script compares queries and domains in targets.

### Fixed

* Separate fasta.Downloader's url building from file downloading functionality.
* Convert arguments to integers appropriately in console_scripts.
* Add help strings to `seed` docs.
* 'None' is no longer a part of downloaded file names.
* Provide unique default path for dumping gml file if one isn't provided.

### Changed

* In `seekr_graph`, change '-l', '--limit' to '-t', '--threshold'.
* Set default community detection to Leiden.
* Remove testing of fasta file downloading.

---

## 1.2.0

### Added

* This CHANGELOG was added!
* Users can see all commands and examples by running "seekr".
* Users can download fasta files from Gencode with `seekr_download`.
* Log2 normalization is now possible and on by default.
* Users can find Louvain based transcript communities with `seekr_graph`.
* Package info is now in `__version__.py` (modeled after `requests`).
  This allows users to do things like `seekr.__version__` in a REPL.

### Fixed

* Removed examples of `k=7` from README.
* Stop building list twice while reading fasta data.

### Changed

* All commands now start with 'seekr_' (e.g. 'pearson' is now 'seekr_pearson')
* README has undergone a large re-write to reflect changes and new features.
* Console commands produce plain text files by default, instead of binary.
* Requirements are all described in setup.py
  requirements.txt has been removed.
