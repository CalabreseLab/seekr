import pytest
import pkg_resources
import networkx
import numpy as np
import pandas as pd

from seekr import graph


class TestMaker:

    def _build_adj(self):
        kmers = 'tests/data/example_2mers.npy'
        kmers = pkg_resources.resource_filename('seekr', kmers)
        kmers = np.load(kmers)
        adj = np.corrcoef(kmers) * -1  # Flip signs for fewer negatives
        names = list(range(5))
        adj = pd.DataFrame(adj, names, names)
        return adj

    def _build_disconnected_graph(self):
        g = networkx.Graph()
        edges = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 4), (2, 5), (2, 6), (7, 8), (8, 9)]
        g.add_edges_from(edges)
        networkx.set_edge_attributes(g, 1, 'weight')
        return g

    def _save_adj(self, out_dir, binary=True):
        adj = self._build_adj()
        if binary:
            out_path = str(out_dir.join('adj.npy'))
            np.save(out_path, adj.values)
        else:
            out_path = str(out_dir.join('adj.csv'))
            adj.to_csv(out_path)
        return out_path

    def _get_maker_with_partition(self, out_path):
        gml_path = str(out_path.join('out_main_sub.gml'))
        adj = self._build_adj()
        maker = graph.Maker(adj, gml_path=gml_path, seed=0)
        maker.build()
        maker.save(True)
        maker.get_partition()
        return maker

    def test_apply_threshold(self):
        adj = self._build_adj()
        maker = graph.Maker(adj)
        maker.apply_threshold()
        assert np.alltrue(maker.adj.values.diagonal() == np.zeros(5))
        assert adj.values[1, 0] == 0
        assert adj.values[1, 2] != 0

    def test_apply_threshold_t1(self):
        adj = self._build_adj()
        maker = graph.Maker(adj, threshold=1)
        maker.apply_threshold()
        assert maker.adj.values.sum() == 0

    def test_apply_threshold_ndarray(self):
        adj = self._build_adj().values
        maker = graph.Maker(adj, threshold=1)
        maker.apply_threshold()
        assert maker.adj.sum() == 0

    def test_build(self):
        adj = self._build_adj()
        maker = graph.Maker(adj)
        maker.build()
        assert type(maker.graph) == networkx.Graph
        assert len(maker.graph) == 5
        assert len(maker.graph.edges()) == 9
        assert len(networkx.get_edge_attributes(maker.graph, 'weight')) == 9
        assert maker.adj is None
        assert maker.main_sub is not None

    def test_build_ndarray(self):
        adj = self._build_adj().values
        maker = graph.Maker(adj)
        maker.build()
        assert type(maker.graph) == networkx.Graph
        assert len(maker.graph) == 5
        assert len(maker.graph.edges()) == 9
        assert len(networkx.get_edge_attributes(maker.graph, 'weight')) == 9

    def test_build_no_clear_adj(self):
        adj = self._build_adj()
        maker = graph.Maker(adj)
        maker.build(clear_adj=False)
        assert maker.adj.equals(adj)

    def test_build_no_main_sub(self):
        adj = self._build_adj()
        maker = graph.Maker(adj)
        maker.build(main_sub=False)
        assert maker.main_sub is None

    def test_find_main_sub(self):
        g = self._build_disconnected_graph()
        maker = graph.Maker()
        maker.graph = g
        maker.find_main_sub()
        assert list(maker.main_sub) == list(range(7))

    def test_save(self, tmpdir):
        gml_path = str(tmpdir.join('out.gml'))
        adj = self._build_adj()
        maker = graph.Maker(adj, gml_path=gml_path)
        maker.build()
        maker.save()
        saved = networkx.read_gml(gml_path)
        expected = [str(n) for n in maker.graph.nodes()]
        assert list(saved.nodes()) == expected
        expected = [(str(n1), str(n2)) for n1, n2 in maker.graph.edges()]
        assert list(saved.edges()) == expected

    def test_save_main_sub(self, tmpdir):
        gml_path = str(tmpdir.join('out_main_sub.gml'))
        adj = self._build_adj()
        maker = graph.Maker(adj, gml_path=gml_path)
        maker.build()
        maker.save(main_sub=True)
        saved = networkx.read_gml(gml_path)
        expected = [str(n) for n in maker.main_sub.nodes()]
        assert list(saved.nodes()) == expected
        expected = [(str(n1), str(n2)) for n1, n2 in maker.main_sub.edges()]
        assert list(saved.edges()) == expected

    def test_get_partition(self, tmpdir):
        maker = self._get_maker_with_partition(tmpdir)
        assert np.isclose(maker.partition.modularity, -0.08024691358024699)
        assert maker.partition.membership == [1, 0, 1, 0, 0]

    def test_membership2attribute(self, tmpdir):
        maker = self._get_maker_with_partition(tmpdir)
        name2group = maker.membership2attribute()
        assert name2group == {'0': 1, '1': 0, '2': 1, '3': 0, '4': 0}
        assert name2group == networkx.get_node_attributes(maker.graph, 'Group')

    def test_membership2attribute_disconnected(self, tmpdir):
        gml_path = str(tmpdir.join('out_main_sub.gml'))
        g = self._build_disconnected_graph()
        maker = graph.Maker(gml_path=gml_path)
        maker.graph = g
        maker.find_main_sub()
        maker.save(True)
        maker.get_partition()
        name2group = maker.membership2attribute()
        expected = {0: 1,
                    1: 1,
                    2: 0,
                    3: 1,
                    4: 0,
                    5: 0,
                    6: 0,
                    7: 2,
                    8: 2,
                    9: 2}
        assert name2group == expected

    def test_membership2attribute_disconnected_ncomms1(self, tmpdir):
        gml_path = str(tmpdir.join('out_main_sub.gml'))
        g = self._build_disconnected_graph()
        maker = graph.Maker(gml_path=gml_path, n_comms=1)
        maker.graph = g
        maker.find_main_sub()
        maker.save(True)
        maker.get_partition()
        name2group = maker.membership2attribute()
        expected = {0: 1,
                    1: 1,
                    2: 0,
                    3: 1,
                    4: 0,
                    5: 0,
                    6: 0,
                    7: 1,
                    8: 1,
                    9: 1}
        assert name2group == expected

    def test_membership2attribute_disconnected_ncomms3(self, tmpdir):
        gml_path = str(tmpdir.join('out_main_sub.gml'))
        g = self._build_disconnected_graph()
        maker = graph.Maker(gml_path=gml_path, n_comms=3, gamma=10)
        maker.graph = g
        maker.find_main_sub()
        maker.save(True)
        maker.get_partition()
        name2group = maker.membership2attribute()
        expected = {0: 0,
                    1: 1,
                    2: 2,
                    3: 3,
                    4: 3,
                    5: 3,
                    6: 3,
                    7: 3,
                    8: 3,
                    9: 3}
        assert name2group == expected

    def test_make_gml_file(self, tmpdir):
        gml_path = str(tmpdir.join('out.gml'))
        csv_path = str(tmpdir.join('out.csv'))
        adj = self._build_adj()
        maker = graph.Maker(adj, gml_path=gml_path, csv_path=csv_path)
        maker.make_gml_csv_files()
        in_graph = networkx.read_gml(gml_path)
        assert list(in_graph.nodes()) == [str(i) for i in range(5)]
        assert len(networkx.get_node_attributes(in_graph, 'Group')) == 5
        df = pd.read_csv(csv_path, index_col=0)
        assert np.alltrue(df.index.values == np.arange(5))
        assert np.alltrue(df['Group'].values == np.array([1, 0, 1, 0, 0]))

