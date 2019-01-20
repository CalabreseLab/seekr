import requests
import gzip
from pathlib import Path

from seekr import fasta


class TestMaker:
    pass #TODO


class TestRandomMaker:
    pass


class TestDownloader:

    def test_find_current_release(self):
        # Note, this information will change over time, so you can't hard code a value.
        # My strategy is to get the value a second way and ensure it's the same.
        expected = fasta.Downloader().find_current_release('human')
        url = f'https://www.gencodegenes.org/human/'
        html = requests.get(url).text
        for line in html.splitlines():
            if line.startswith('<h1>'):
                header = line
                break
        release = header.split()[1]
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

    def test_get_gencode(self, tmpdir):
        downloader = fasta.Downloader()
        out_path = Path(tmpdir, 'human_lncs.fa.gz')
        unzipped = Path(tmpdir, 'human_lncs.fa')
        downloader.get_gencode('lncRNA', out_path=str(out_path))
        assert not out_path.exists()
        assert unzipped.exists()
        with unzipped.open() as in_file:
            line = next(in_file)
            assert line.startswith('>')

    def test_get_gencode_mouse(self, tmpdir):
        downloader = fasta.Downloader()
        out_path = Path(tmpdir, 'mouse_all.fa.gz')
        unzipped = Path(tmpdir, 'mouse_all.fa')
        downloader.get_gencode('all', 'mouse', out_path=str(out_path))
        assert not out_path.exists()
        assert unzipped.exists()
        with unzipped.open() as in_file:
            line = next(in_file)
            assert line.startswith('>')

    def test_get_gencode_pc22(self, tmpdir):
        downloader = fasta.Downloader()
        out_path = Path(tmpdir, 'human_pc_22.fa.gz')
        unzipped = Path(tmpdir, 'human_pc_22.fa')
        downloader.get_gencode('pc', release='22', out_path=str(out_path))
        assert not out_path.exists()
        assert unzipped.exists()
        with unzipped.open() as in_file:
            # Since this file is both "locked" and unzipped, let's check the line count.
            count = len(in_file.readlines())
            assert count == 187052

    def test_get_gencode_zipped(self, tmpdir):
        downloader = fasta.Downloader()
        out_path = Path(tmpdir, 'mouse_lncs.fa.gz')
        downloader.get_gencode('lncRNA', 'mouse', 'M5', out_path=str(out_path), unzip=False)
        assert out_path.exists()
