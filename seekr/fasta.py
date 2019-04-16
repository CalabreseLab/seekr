# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:13:42 2015

@author: jessime

Warnings
--------
* A lot of this code is dependent on GENCODE formatted fasta files.
* The C code in the shuffle function will segfault if the string is too big.
This will cause the Python Kernel to crash.

Issues
------
1. Filters are too similar. Shouldn't repeat that code. Just send in the filter condition.
2. Maker and Extractor share code that maybe should be in parent class
3. Do `if outfile: save \ always return`
"""

import os
import gzip
import pickle
import shutil
import ftplib
import requests
import urllib.request
import numpy as np

from contextlib import closing
from os.path import exists, join
from os import makedirs
from ushuffle import shuffle
from ushuffle import set_seed

from seekr.my_tqdm import my_tqdm, my_trange
from seekr.fasta_reader import Reader


class Maker:
    """Manipulates fasta files into new fasta files.

    Parameters
    ----------
    infasta : str (default=None)
        Name of input fasta file to be manipulated
    outfasta : str (default=None)
        Location to store filtered fasta file
    outnames : str (default=None)
        Location to stores saved transcript names

    Attributes
    ----------
    data : list
        Zipped list of names and seqs tuples
    names : list
        Header lines of fasta file
    seqs : list
        Sequences from fasta file
    """

    def __init__(self, infasta=None, outfasta=None, outnames=None):
        self.infasta = infasta
        if infasta is not None:
            self.data, self.names, self.seqs = Reader(infasta).get_data()
        self.outfasta = outfasta
        self.outnames = outnames

    def _name_dump(self, common_name_list):
        """Save new names if necessary

        Parameters
        ----------
        common_name_list : list
            Common style names of sequences which are being kept
        """
        if self.outnames is not None:
            pickle.dump(common_name_list, open(self.outnames, 'wb'))

    def filter1(self, zeros=1, unique_per_gene=False):
        """Filters a fasta file to only keep *1 transcript names.

        Parameters
        ----------
        zeros: int (default=1)
            The number of zeros preceding `1`.
            (e.g. zeros=1 will filter for `01` transcripts. zeros=2 will filter for `001`)
        unique_per_gene: bool (default=False)
            If True, when multiple 01 isoforms exist per gene, only keep the smaller name.
            (e.g. if XIST-001 and XIST-201 both exist, keep XIST-001)

        Returns
        -------
        common_name_list: List[str]
            The common style names of the sequences which are being kept
        """
        common_name_list = []
        with open(self.outfasta, 'w') as outfasta:
            current_gene = ''
            best_name_for_gene = ''
            best_seq_for_gene = ''
            for name, seq in self.data:
                name_data = name.split('|')
                common_name = name_data[4]
                gene = name_data[1]
                canonical = common_name[-(zeros+1):] == '0'*zeros+'1'
                if canonical:
                    if not unique_per_gene:
                        common_name_list.append(common_name)
                        outfasta.write(name+'\n')
                        outfasta.write(seq+'\n')
                        continue
                    new_gene = gene != current_gene
                    if not best_name_for_gene:
                        best_name_for_gene = name
                        best_seq_for_gene = seq
                        current_gene = gene
                    elif best_seq_for_gene and not new_gene:
                        new_is_better = common_name[-3:] < best_name_for_gene.split('|')[4][-3:]
                        if new_is_better:
                            best_name_for_gene = name
                            best_seq_for_gene = seq
                        warning = (f'Gene {gene} has at least two viable isoforms. '
                                   f'Keeping: {best_name_for_gene}')
                        print(warning)
                    elif new_gene:
                        common_name_list.append(best_name_for_gene.split('|')[4])
                        outfasta.write(best_name_for_gene+'\n')
                        outfasta.write(best_seq_for_gene + '\n')
                        current_gene = gene
                        best_name_for_gene = name
                        best_seq_for_gene = seq
            if best_seq_for_gene:
                common_name_list.append(common_name)
                outfasta.write(best_name_for_gene + '\n')
                outfasta.write(best_seq_for_gene + '\n')
        self._name_dump(common_name_list)
        return common_name_list

    def filter_size(self, size_lim, keep_all_below=True):
        """Filters fasta file based on size of each transcript

        Parameters
        ----------
            size_lim : int
                the size limit to keep above or below
            keep_all_below : bool (default=True)
                Keep sequences longer or shorter than specified length?

        Returns
        -------
        filtered_fasta : list
            Lines of the new fasta file
        """
        common_name_list = []
        filtered_fasta = []
        for name, seq in self.data:
            common_name = name.split('|')[4]
            length = len(seq)
            if ((length <= size_lim and keep_all_below) or
            (length >= size_lim and not keep_all_below)):
                if self.outnames:
                    common_name_list.append(common_name)
                filtered_fasta.append((name, seq))

        if self.outfasta is not None:
            with open(self.outfasta, 'w') as outfasta:
                for dt in filtered_fasta:
                    outfasta.write(dt+'\n')
        self._name_dump(common_name_list)
        return filtered_fasta

    def filter_name(self, keep_names):
        """Filters fasta file based on a list of common names

        Parameters
        ----------
        keep_names : list
            The common style names to be kept from the fasta file

        Returns
        -------
        filtered_fasta : list
            Lines of the new fasta file
        """
        common_names = [n.split('|')[4] for n in self.names]
        names_dict = dict(zip(common_names, self.data))
        filtered_fasta = []
        for n in keep_names:
            try:
                filtered_fasta.append(names_dict[n])
            except KeyError:
                print('{} is not in fasta file'.format(n))

        if self.outfasta is not None:
            with open(self.outfasta, 'w') as outfasta:
                for dt in filtered_fasta:
                    outfasta.write(dt[0]+'\n')
                    outfasta.write(dt[1]+'\n')
        return filtered_fasta

    def separate(self, outdir, filenames=None):
        """Separate a fasta file into individual fasta files

        Parameters
        ----------
        outdir : str
            Path of directory in which to save the new fasta files
        filenames : [str] (default=None)
            Optional list of names for files instead of integers.
        """
        if not exists(outdir):
            makedirs(outdir)

        if filenames is None:
            filenames = range(len(self.data))
        for fn, (name, seq) in zip(filenames, self.data):
            fn = '{}.fa'.format(fn)
            with open(join(outdir, fn), "w") as outfile:
                outfile.write(name+'\n')
                outfile.write(seq+'\n')


class RandomMaker(Maker):
    """Provides functions for creating randomly generated fasta files

    Parameters
    ----------
    infasta: str (default=None)
        Name of input fasta file to be manipulated
    outfasta: str (default=None)
        Location to store filtered fasta file
    k: int (default=1)
        The size of kmer to conserve between original and random sequences
    mutations: int
        Number of SNP mutations to make in sequence
    seed: int (default=None)
        Seed to use for reproducible shuffling with ushuffle
    individual: bool (default=True)
        Whether to conserve kmers of each sequence or the entire fasta file

    Attributes
    ----------
    error: float
        The mean error in kmer frequency betwen original counts and new counts
    weights: DataFrame
        Shows fractions of each kmer for old and new seqs, along with error
    new_fasta_seqs: list
        New, randomized sequences with conserved kmer frequencies

    Notes
    -----
    The kmers are conserved using the uShuffle algorithm described here:
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192
    The C code can be found here:
        https://github.com/dcjones/cbgb/tree/master/ushuffle
    The Cython code can be found here:
        https://github.com/guma44/ushuffle
    """

    def __init__(self, infasta=None, outfasta=None, k=1, mutations=0, seed=None, individual=True):
        super(RandomMaker, self).__init__(infasta=infasta, outfasta=outfasta)
        self.k = k
        self.mutations = mutations
        self.individual = individual
        self.weights = None
        self.new_fasta_seqs = None
        if seed is not None:
            set_seed(seed)
            np.random.seed(seed)

    def shuffle(self, seq):
        rand_seq = shuffle(seq.encode(), self.k).decode('utf-8')
        rand_array = np.array(list(rand_seq))
        indices = np.random.choice(len(seq), self.mutations, replace=False)
        new_bases = np.random.choice(list('AGTC'), self.mutations)
        rand_array[indices] = new_bases
        rand_seq = ''.join(rand_array)
        return rand_seq

    def get_random_seqs(self, seqs):
        """Generates random RNA sequences

        Parameters
        ----------
        seqs : [str]
            Original sequences to create random copies of

        Returns
        -------
        new_seqs : [str]
            Randomly generated sequences
        """
        new_seqs = [self.shuffle(s) for s in my_tqdm()(seqs)]
        return new_seqs

    def split(self, seq):
        """Return concatenated sequence to individual sequence lengths"""
        lengths = [len(s) for s in self.seqs]
        new_seqs = []
        for l in lengths:
            new_seqs.append(seq[:l])
            seq = seq[l:]
        return new_seqs

    def inject_seqs(self, new_seqs):
        """Use the original fasta file headers to put the random sequences in fasta format

        Parameters
        ----------
        new_seqs : list
            string elements, each a randomly generated RNA

        Returns
        -------
        new_fasta_seqs : list
            The lines of a new fasta file filled with random RNAs
        """
        new_fasta_seqs = []
        for name, new_seq in zip(self.names, new_seqs):
            new_fasta_seqs.append(name)
            new_fasta_seqs.append(new_seq)
        return new_fasta_seqs

    def save(self):
        if self.outfasta is not None:
            with open(self.outfasta, 'w') as outfasta:
                for line in self.new_fasta_seqs:
                    outfasta.write('{}\n'.format(line))

    def synthesize_random(self):
        """Make random RNAs based on natural RNAs.

        Returns
        -------
        new_fasta_seqs : list
            Lines of the new fasta file
        """
        if self.individual:
            new_seqs = self.get_random_seqs(self.seqs)
        else:
            new_seqs = self.get_random_seqs([''.join(self.seqs)])
            new_seqs = self.split(new_seqs[0])
        self.new_fasta_seqs = self.inject_seqs(new_seqs)
        self.save()
        return self.new_fasta_seqs


class Downloader:
    """Download fasta files from a variety of online sources

    #TODO Expand to other sources of fasta files.
    """

    def __init__(self):
        pass

    def find_current_release(self, species):
        """Scrape Genecode's site to find the latest release value.

        Parameters
        ----------
        species: str
            Name of species (human or mouse)
        """
        url = f'https://www.gencodegenes.org/{species}/'
        html = requests.get(url).text
        for line in html.splitlines():
            if '<title>' in line:
                title = line
                break
        release = title.split('Release')[1].strip().strip('</title>')
        return release

    def build_url(self, biotype, species, release):
        """Build the correct ftp URL to download from GENCODE.

        Parameters
        ----------
        biotype: str
            Name of Genocde set to download. Must be one of ('all', 'pc', 'lncRNA').
        species: str (default='human')
            Name of species. Must be one of: ('human' or 'mouse').
        release: str (default=None)
            Name of specific release to download (e.g. 'M5'). If None, download latest release.

        Returns
        -------
        url: str
            FTP file to download.
        release: str
            Name of specific release to download (e.g. 'M5'). Will not be None.
        """
        error_msg = "'biotype' must be in ('all', 'pc', 'lncRNA')."
        assert biotype in ('all', 'pc', 'lncRNA'), error_msg
        error_msg = "'species' must be either 'human' or 'mouse'."
        assert species in ('human', 'mouse'), error_msg
        biotype2prefix = {'all': '', 'pc': 'pc_', 'lncRNA': 'lncRNA_'}
        prefix = biotype2prefix[biotype]
        if release is None:
            release = self.find_current_release(species)
        if species == 'mouse':
            error_msg = "Mouse releases must begin with 'M'."
            assert release[0] == 'M', error_msg
        url_base = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_'
        url = f'{species}/release_{release}/gencode.v{release}.{prefix}transcripts.fa.gz'
        url = url_base + url
        return url, release

    def gunzip(self, gzip_path):
        """Unzip a gzipped file and remove orginal.

        Paramters
        ---------
        gzip_path: str
            Gzipped file location
        """
        out_path = gzip_path.strip('.gz')
        with gzip.open(gzip_path, 'rb') as in_file:
            with open(out_path, 'wb') as out_file:
                shutil.copyfileobj(in_file, out_file)
        os.remove(gzip_path)

    def get_gencode(self, biotype, species='human', release=None, out_path=None, unzip=True):
        """Download .fa.gz file from Gencode's site.

        Parameters
        ----------
        biotype: str
            Name of Genocde set to download. Must be one of ('all', 'pc', 'lncRNA').
        species: str (default='human')
            Name of species. Must be one of: ('human' or 'mouse').
        release: str (default=None)
            Name of specific release to download (e.g. 'M5'). If None, download latest release.
        out_path: str (default=None)
            Path to location for fasta file. Default will save by release name.
        unzip: bool (default=True)
            If False, do not gunzip fasta file after downloading
        """
        url, release = self.build_url(biotype, species, release)
        if out_path is not None:
            error_msg = "Even if unzipping, 'out_path' must end with '.gz'."
            assert out_path.endswith('.gz'), error_msg
        try:
            with closing(urllib.request.urlopen(url)) as r:
                if out_path is None:
                    out_path = f'v{release}_{biotype}.fa.gz'
                with open(out_path, 'wb') as out_file:
                    shutil.copyfileobj(r, out_file)
            if unzip:
                self.gunzip(out_path)
        except urllib.error.URLError as url_error:
            print('The file failed to download because:\n', url_error)
            cd_err = "<urlopen error ftp error: error_perm('550 Failed to change directory.',)>"
            if str(url_error) == cd_err:
                print('Did you pass a valid `--release` value (e.g. M14, 22)?')

