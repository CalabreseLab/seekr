import sys
import requests
import gzip
import pkg_resources

from pathlib import Path

import pytest

from seekr import fasta


class TestMaker:

    def test_filter01(self, tmpdir):
        infasta = 'tests/data/v22_pc_head.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        maker = fasta.Maker(infasta, str(Path(tmpdir, 'out.fa')))
        names = maker.filter1()
        assert names == ['OR4F5-001', 'FO538757.3-201', 'FO538757.2-201']

    def test_filter001(self, tmpdir):
        infasta = 'tests/data/v22_pc_head.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        maker = fasta.Maker(infasta, str(Path(tmpdir, 'out.fa')))
        names = maker.filter1(2)
        assert names == ['OR4F5-001']

    def test_filter01_1per(self, tmpdir, capsys):
        infasta = 'tests/data/example2.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        maker = fasta.Maker(infasta, str(Path(tmpdir, 'out.fa')))
        names = maker.filter1(unique_per_gene=True)
        assert names == ['JK-001', 'JK2-001']
        captured = capsys.readouterr()
        expected = ('Gene ENSG1 has at least two viable isoforms. '
                    'Keeping: >ENST1|ENSG1|OTTHUMG1|OTTHUMT1|JK-001|JK|918|CDS:1-918|\n'
                    'Gene ENSG2 has at least two viable isoforms. '
                    'Keeping: >ENST4|ENSG2|OTTHUMG2|OTTHUMT4|JK2-001|JK2|918|CDS:1-918|\n')
        assert captured.out == expected


class TestRandomMaker:

    @pytest.mark.skipif(sys.platform == 'darwin', reason='Random seed fails on Mac')
    def test_shuffle1(self):
        rand_maker = fasta.RandomMaker(seed=1)
        rand_seq = rand_maker.shuffle('AGTCAGTC')
        assert rand_seq == 'TTGGACAC'

    def test_shuffle2(self):
        rand_maker = fasta.RandomMaker(k=2, seed=1)
        rand_seq = rand_maker.shuffle('AGTCAGTCAGTCAGTC')
        assert rand_seq == 'AGTCAGTCAGTCAGTC'
        kmer_counts = {
            'AA': 0,
            'AG': 4,
            'AT': 0,
            'AC': 0,
            'GA': 0,
            'GG': 0,
            'GT': 4,
            'GC': 0,
            'TA': 0,
            'TG': 0,
            'TT': 0,
            'TC': 4,
            'CA': 3,
            'CG': 0,
            'CT': 0,
            'CC': 0
        }
        rand_kmer_counts = {k:rand_seq.count(k) for k in kmer_counts}
        assert kmer_counts == rand_kmer_counts

    def test_mutations(self):
        rand_maker = fasta.RandomMaker(seed=1)
        rand_maker.mutations = 2
        rand_seq = rand_maker.shuffle('AAAAAA')
        assert rand_seq == 'ACGAAA'

    @pytest.mark.skipif(sys.platform == 'darwin', reason='Random seed fails on Mac')
    def test_get_random_seqs(self):
        rand_maker = fasta.RandomMaker(k=2, seed=1)
        seqs = [
            'AGTCAGTCAGTCAGTC',
            'ATGATATATATATGAT',
            'ATGCATAGTTTTTTTTTCTGC'
        ]
        rand_seqs = rand_maker.get_random_seqs(seqs)
        expected = ['AGTCAGTCAGTCAGTC',
                    'ATGATGATATATATAT',
                    'ATTGCTTTGTTTAGCATTTTC']
        assert rand_seqs == expected

    def test_split(self):
        rand_maker = fasta.RandomMaker(k=2, seed=1)
        seqs = [
            'this sentence is 35 characters long',
            'this one is 14'
        ]
        rand_maker.seqs = seqs
        seq = 'TCATTAAGCGCGTCGGTCTCTGTGTACGTCATCTCCATTTTTTTTCGTG'
        rand_seqs = rand_maker.split(seq)
        expected = ['TCATTAAGCGCGTCGGTCTCTGTGTACGTCATCTC',
                    'CATTTTTTTTCGTG']
        assert rand_seqs == expected
        assert len(expected[0]) == 35
        assert len(expected[1]) == 14

    def test_inject_seqs(self):
        rand_maker = fasta.RandomMaker()
        rand_maker.names = ['>seq1', '>seq2']
        new_seqs = ['this is new', 'also new']
        new_fasta_seqs = rand_maker.inject_seqs(new_seqs)
        expected = ['>seq1', 'this is new', '>seq2', 'also new']
        assert new_fasta_seqs == expected

    @pytest.mark.skipif(sys.platform == 'darwin', reason='Random seed fails on Mac')
    def test_synthesize_random(self, tmpdir):
        infasta = 'tests/data/example.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        outfasta = str(tmpdir.join('rand.fa'))
        rand_maker = fasta.RandomMaker(infasta, outfasta, seed=1)
        rand_maker.synthesize_random()
        with open(outfasta) as outfasta:
            rand_fasta = ''.join(next(outfasta) for i in range(6))
        expected = ('>SEQ1\n'
                    'AAAAAA\n'
                    '>SEQ2\n'
                    'GTGTGTGTGTTG\n'
                    '>SEQ3\n'
                    'ATTACTTGGCACCGGA\n')
        assert rand_fasta == expected


class TestDownloader:

    def _expected_current_release(self, species):
        # Note, this information will change over time, so you can't hard code a value.
        # My strategy is to get the value a second way and ensure it's the same.
        url = f'https://www.gencodegenes.org/{species}/'
        html = requests.get(url).text
        for line in html.splitlines():
            if line.startswith('<h1>'):
                header = line
                break
        release = header.split()[1]
        return release

    def test_find_current_release(self):
        release = fasta.Downloader().find_current_release('human')
        expected = self._expected_current_release('human')
        assert expected == release

    def test_gunzip(self, tmpdir):
        test_path = Path(tmpdir, 'test.txt.gz')
        unzipped = Path(tmpdir, 'test.txt')
        content = 'Hello, World!'.encode()
        with gzip.open(test_path, 'wb') as test_file:
            test_file.write(content)
        fasta.Downloader().gunzip(str(test_path))
        assert not test_path.exists()
        assert unzipped.exists()
        assert unzipped.read_text() == content.decode()

    def test_build_url(self):
        downloader = fasta.Downloader()
        url, release = downloader.build_url('lncRNA', 'human', None)
        exp_release = self._expected_current_release('human')
        expected = ('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_'
                    f'human/release_{exp_release}/gencode.v{exp_release}.lncRNA_transcripts.fa.gz')
        assert url == expected

    def test_build_url_mouse(self):
        downloader = fasta.Downloader()
        url, release = downloader.build_url('all', 'mouse', None)
        exp_release = self._expected_current_release('mouse')
        expected = ('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_'
                    f'mouse/release_{exp_release}/gencode.v{exp_release}.transcripts.fa.gz')
        assert url == expected

    def test_build_url_pc22(self):
        downloader = fasta.Downloader()
        url, release = downloader.build_url('pc', 'human', '22')
        expected = ('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_'
                    'human/release_22/gencode.v22.pc_transcripts.fa.gz')
        assert url == expected

    # TODO Consider using mocking to help with testing downloads, or move to integration tests...
    # def test_get_gencode(self, tmpdir):
    #     downloader = fasta.Downloader()
    #     out_path = Path(tmpdir, 'human_lncs.fa.gz')
    #     unzipped = Path(tmpdir, 'human_lncs.fa')
    #     downloader.get_gencode('lncRNA', out_path=str(out_path))
    #     assert not out_path.exists()
    #     assert unzipped.exists()
    #     with unzipped.open() as in_file:
    #         line = next(in_file)
    #         assert line.startswith('>')
    #
    # def test_get_gencode_mouse(self, tmpdir):
    #     downloader = fasta.Downloader()
    #     out_path = Path(tmpdir, 'mouse_all.fa.gz')
    #     unzipped = Path(tmpdir, 'mouse_all.fa')
    #     downloader.get_gencode('all', 'mouse', out_path=str(out_path))
    #     assert not out_path.exists()
    #     assert unzipped.exists()
    #     with unzipped.open() as in_file:
    #         line = next(in_file)
    #         assert line.startswith('>')
    #
    # def test_get_gencode_pc22(self, tmpdir):
    #     downloader = fasta.Downloader()
    #     out_path = Path(tmpdir, 'human_pc_22.fa.gz')
    #     unzipped = Path(tmpdir, 'human_pc_22.fa')
    #     downloader.get_gencode('pc', release='22', out_path=str(out_path))
    #     assert not out_path.exists()
    #     assert unzipped.exists()
    #     with unzipped.open() as in_file:
    #         # Since this file is both "locked" and unzipped, let's check the line count.
    #         count = len(in_file.readlines())
    #         assert count == 187052
    #
    # def test_get_gencode_zipped(self, tmpdir):
    #     downloader = fasta.Downloader()
    #     out_path = Path(tmpdir, 'mouse_lncs.fa.gz')
    #     downloader.get_gencode('lncRNA', 'mouse', 'M5', out_path=str(out_path), unzip=False)
    #     assert out_path.exists()
