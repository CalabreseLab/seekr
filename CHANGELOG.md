# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

---

## [Current]

This is where you should add changes as you go along.
They will be shuttled into a version number just before the version is released.

### Added

* `seekr_canonical_gencode` command line script filters for -001 transcripts.
* Example integration script in the test directory.
* Add legacy option to continue using Louvain instead of Leiden.
* Travis CI automatic push testing.
* `seekr_visualize_distro` command makes distribution of r-values.

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
