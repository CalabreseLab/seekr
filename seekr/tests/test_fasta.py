import requests
import gzip
import pkg_resources

from pathlib import Path

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


class TestRandomMaker:
    pass
    # TODO Add testing


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
