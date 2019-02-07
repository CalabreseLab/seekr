import pytest
import pkg_resources
import networkx
import numpy as np
import pandas as pd

from seekr import utils


class TestMaker:

    def _save_adj(self, out_dir, binary=True):
        adj = self._build_adj()
        if binary:
            out_path = str(out_dir.join('adj.npy'))
            np.save(out_path, adj.values)
        else:
            out_path = str(out_dir.join('adj.csv'))
            adj.to_csv(out_path)
        return out_path

    def _build_adj(self):
        kmers = 'tests/data/example_2mers.npy'
        kmers = pkg_resources.resource_filename('seekr', kmers)
        kmers = np.load(kmers)
        adj = np.corrcoef(kmers) * -1  # Flip signs for fewer negatives
        names = list(range(5))
        adj = pd.DataFrame(adj, names, names)
        return adj

    def test_get_adj_ndarray(self):
        expected = self._build_adj().values
        adj = utils.get_adj(expected)
        assert type(adj) == np.ndarray

    def test_get_adj_df(self):
        expected = self._build_adj()
        adj = utils.get_adj(expected)
        assert type(adj) == pd.DataFrame

    def test_get_adj_str_ndarray(self, tmpdir):
        adj_path = self._save_adj(tmpdir)
        adj = utils.get_adj(adj_path)
        assert type(adj) == np.ndarray

    def test_get_adj_str_df(self, tmpdir):
        adj_path = self._save_adj(tmpdir, binary=False)
        adj = utils.get_adj(adj_path)
        assert type(adj) == pd.DataFrame