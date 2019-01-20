# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

---

## [Current]

This is where you should add changes as you go along.
They will be shuttled into a version number just before the version is released.

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
