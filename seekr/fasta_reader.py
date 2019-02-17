# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 13:10:49 2016

@author: jessime
"""


class Reader():
    """Fixes any compatibility issues a fasta file might have with this code.

    Parameters
    ----------
    infasta : str (default=None)
        Name of input fasta file to be manipulated
    outfasta : str (default=None)
        location to store extracted data from infasta
    names : iter (default=None)
        Common style names to use in header lines

    Attributes
    ----------
    data : list
        Raw lines of the infasta file
        Note: This is different than the data attribute in other classes

    Examples
    --------
    Putting the sequence on one line instead of breaking it every 80 chars.
    Making sure the whole sequence is capitalized.
    Restructuring the name line to work with GENCODE's naming.
    """

    def __init__(self, infasta=None, outfasta=None, names=None):
        self.infasta = infasta
        self.outfasta = outfasta
        self.names = names

        self.data = None

    def _read_data(self):
        """Sets data to stripped lines from the fasta file
        """
        with open(self.infasta) as infasta:
            self.data = [l.strip() for l in infasta]

    def _upper_seq_per_line(self):
        """Sets data to upper case, single line sequences for each header
        """
        new_data = []
        seq = ''
        for i, line in enumerate(self.data):
            if line[0] == '>':
                if seq:
                    new_data.append(seq.upper())
                    seq = ''
                else:
                    assert i == 0, 'There may be a header without a sequence at line {}.'.format(i)
                new_data.append(line)
            else:
                seq += line
        new_data.append(seq.upper())
        self.data = new_data

    def get_lines(self):
        self._read_data()
        self._upper_seq_per_line()
        return self.data

    def get_seqs(self):
        clean_data = self.get_lines()
        seqs = clean_data[1::2]
        return seqs

    def get_headers(self):
        clean_data = self.get_lines()
        headers = clean_data[::2]
        return headers

    def get_data(self, tuples_only=False):
        clean_data = self.get_lines()
        headers = clean_data[::2]
        seqs = clean_data[1::2]
        tuples = zip(headers, seqs)
        if tuples_only:
            return tuples
        else:
            return tuples, headers, seqs

    def supply_basic_header(self):
        """Convert headerlines to GENCODE format with only common name and length"""
        new_fasta = []
        
        if self.names is None:
            self.names = iter(self.get_headers())
        for i, line in enumerate(self.data):
            if line[0] =='>':
                name = next(self.names).strip('>')
                length = len(self.data[i+1])
                new_fasta.append('>||||{}||{}|'.format(name, length))
            else:
                new_fasta.append(line)
        return new_fasta

    def save(self):
        """Write self.data to a new fasta file"""
        with open(self.outfasta, 'w') as outfasta:
            for line in self.data:
                outfasta.write(line+'\n')
