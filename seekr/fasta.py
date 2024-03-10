# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:13:42 2015

@author: jessime

Warnings
--------
* A lot of this code is dependent on GENCODE formatted fasta files.

"""

import os
import gzip
# import pickle
import shutil
import ftplib
import requests
import urllib.request

from contextlib import closing
# from os.path import exists, join
# from os import makedirs


# from seekr.my_tqdm import my_tqdm, my_trange
# from seekr.fasta_reader import Reader


class Downloader:
    """Download fasta and gtf files from gencode
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
        url = f"https://www.gencodegenes.org/{species}/"
        html = requests.get(url).text
        for line in html.splitlines():
            if "<title>" in line:
                title = line
                break
        release = title.split("Release")[1].strip().strip("</title>")
        return release

    def build_url(self, biotype, species, gtf, release):
        """Build the correct ftp URL to download from GENCODE.

        Parameters
        ----------
        biotype: str
            Name of Genocde set to download. Must be one of ('all', 'pc', 'lncRNA').
        species: str (default='human')
            Name of species. Must be one of: ('human' or 'mouse').
        gtf: bool (default=False)
            If True, download the most comprehensive GTF file of the same species and release.
        release: str (default=None)
            Name of specific release to download (e.g. 'M5'). If None, download latest release.

        Returns
        -------
        url: str
            FTP file to download.
        gtf_url: str
            URL to download GTF file.
        release: str
            Name of specific release to download (e.g. 'M5'). Will not be None.
        """
        error_msg = "'biotype' must be in ('all', 'pc', 'lncRNA')."
        assert biotype in ("all", "pc", "lncRNA"), error_msg
        error_msg = "'species' must be either 'human' or 'mouse'."
        assert species in ("human", "mouse"), error_msg
        biotype2prefix = {"all": "", "pc": "pc_", "lncRNA": "lncRNA_"}
        prefix = biotype2prefix[biotype]
        if release is None:
            release = self.find_current_release(species)
        if species == "mouse":
            error_msg = "Mouse releases must begin with 'M'."
            assert release[0] == "M", error_msg
        url_base = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_"
        url = f"{species}/release_{release}/gencode.v{release}.{prefix}transcripts.fa.gz"
        url = url_base + url
        if gtf:
            gtf_url = f"{species}/release_{release}/gencode.v{release}.chr_patch_hapl_scaff.annotation.gtf.gz"
            gtf_url = url_base + gtf_url
        else:
            gtf_url = None

        return url, gtf_url, release

    def gunzip(self, gzip_path):
        """Unzip a gzipped file and remove orginal.

        Paramters
        ---------
        gzip_path: str
            Gzipped file location
        """
        out_path = gzip_path.strip(".gz")
        with gzip.open(gzip_path, "rb") as in_file:
            with open(out_path, "wb") as out_file:
                shutil.copyfileobj(in_file, out_file)
        os.remove(gzip_path)

    def get_gencode(self, biotype, species="human", gtf=False, release=None, fasta_path=None, gtf_path=None, unzip=True):
        """Download .fa.gz file and/or .gtf file from Gencode's site.

        Parameters
        ----------
        biotype: str
            Name of Genocde set to download. Must be one of ('all', 'pc', 'lncRNA').
        species: str (default='human')
            Name of species. Must be one of: ('human' or 'mouse').
        gtf: bool (default=False)
            If True, download the most comprehensive GTF file of the same species and release.
            gencode.v{release}.chr_patch_hapl_scaff.annotation.gtf.gz
        release: str (default=None)
            Name of specific release to download (e.g. 'M5'). If None, download latest release.
        fasta_path: str (default=None)
            Path to location for fasta file. Default will save by release name.
        gtf_path: str (default=None)
            Path to location for gtf file. Default will save by release name.
        unzip: bool (default=True)
            If False, do not gunzip fasta file after downloading
        """
        url, gtf_url, release = self.build_url(biotype, species, gtf, release)
        
        if fasta_path is not None:
            error_msg = "Even if unzipping, 'fasta_path' must end with '.gz'."
            assert fasta_path.endswith(".gz"), error_msg
        
        if gtf_path is not None:
            error_msg = "Even if unzipping, 'gtf_path' must end with '.gz'."
            assert gtf_path.endswith(".gz"), error_msg

        try:
            with closing(urllib.request.urlopen(url)) as r:
                if fasta_path is None:
                    fasta_path = f"v{release}_{biotype}.fa.gz"
                with open(fasta_path, "wb") as out_file:
                    shutil.copyfileobj(r, out_file)
            if unzip:
                self.gunzip(fasta_path)

            if gtf:
                with closing(urllib.request.urlopen(gtf_url)) as r:
                    if gtf_path is None:
                        gtf_path = f"v{release}_{biotype}.chr_patch_hapl_scaff.annotation.gtf.gz"
                    with open(gtf_path, "wb") as out_file:
                        shutil.copyfileobj(r, out_file)
                if unzip:
                    self.gunzip(gtf_path)
                    
        except urllib.error.URLError as url_error:
            print("The file failed to download because:\n", url_error)
            cd_err = "<urlopen error ftp error: error_perm('550 Failed to change directory.',)>"
            if str(url_error) == cd_err:
                print("Did you pass a valid `--release` value (e.g. M14, 22)?")
