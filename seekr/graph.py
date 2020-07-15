#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import uuid
import numpy as np
import pandas as pd
import igraph
import networkx
import louvain
import leidenalg

from pathlib import Path

from seekr.utils import get_adj


class Maker:
    """Generate and manipulate network graphs using networkx and igraph.

    Parameters
    ----------
    adj: str | ndarray | DataFrame
        (Path to) adjacency matrix describing weighted edges between nodes.
    gml_path: str
        Path to output file for new network.
    csv_path: str
        Path to two column file storing transcript name and membership.
    leiden: bool (default=True)
        If set to False, use the Louvain algorithm for detection instead of Leiden.
    threshold: float (default=0)
        Value for thresholding adjacency matrix. Below this threshold, all edges are 0.
    gamma: float (default=1)
        Resolution parameter for community detection algorithm.
    n_comms: int (default=5)
        Number of communities to find. This does not count a null community.
    seed: int (default=None)
        Set the seed for reproducible results from community detection.

    Attributes
    ----------
    graph: networkx.Graph
        Primary graph build from adjacency matrix
    main_sub: networkx.Graph
        Largest connected subgraph of primary graph
    partition: RBConfigurationVertexPartition
        Contains the membership data for each node
    df: DataFrame
        Stores transcript name and membership
    tempfile_path: Path
        Location of temporary gml file if gml_path is not specified.
        Needed for converting networkx Graph to igraph Graph.
    """

    def __init__(self, adj=None, gml_path=None, csv_path=None, leiden=True,
                 threshold=0, gamma=1, n_comms=5, seed=None):
        self.adj = adj
        if adj is not None:
            self.adj = get_adj(adj)
        self.gml_path = gml_path
        self.csv_path = csv_path
        self.detector = leidenalg if leiden else louvain
        self.threshold = threshold
        self.gamma = gamma
        self.n_comms = n_comms
        self.seed = seed

        self.graph = None
        self.main_sub = None
        self.partition = None
        self.df = None
        self.tempfile_path = None
        if self.gml_path is None:
            self.tempfile_path = Path(str(uuid.uuid4()) + '.gml')

    def apply_threshold(self):
        """Remove low weighted edges from graph."""
        adj = self.adj
        if isinstance(adj, pd.DataFrame):
            adj = adj.values
        np.fill_diagonal(adj, 0)
        adj[adj < self.threshold] = 0

    def find_main_sub(self):
        """Calculate the largest connected subgraph

        Notes
        -----
        https://en.wikipedia.org/wiki/Connected_component_(graph_theory)
        """
        subgraphs = [self.graph.subgraph(c) for c in networkx.connected_components(self.graph)]
        graph_sizes = [sub.size() for sub in subgraphs]
        self.main_sub = subgraphs[graph_sizes.index(max(graph_sizes))]

    def save(self, main_sub=False):
        """Saves graph in the .gml format

        Parameters
        ----------
        main_sub: bool (default=False)
            If True, the graph stored in self.main_sub will be saved
        """
        graph = self.main_sub if main_sub else self.graph
        # This gives a way to dump the main_sub even when a gml_path hasn't been assigned.
        gml_path = self.gml_path if self.gml_path else self.tempfile_path
        networkx.write_gml(graph, gml_path)

    def get_partition(self, gml_path=None):
        """Run community detection on igraph Graph from .gml file.

        Parameters
        ----------
        gml_path: str (default=None)
            Alternative location for loading in .gml file. Defaults to self.gml_path.
        """
        paths = (gml_path, self.gml_path, self.tempfile_path)
        gml_path = next(path for path in paths if path is not None)
        graph = igraph.Graph.Read_GML(str(gml_path))
        extra_args = {}
        if self.seed is not None:
            if self.detector is louvain:
                self.detector.set_rng_seed(self.seed)
            else:
                extra_args['seed'] = self.seed
        self.partition = self.detector.find_partition(graph,
                                                      self.detector.RBConfigurationVertexPartition,
                                                      weights='weight',
                                                      resolution_parameter=self.gamma,
                                                      **extra_args)

    def membership2attribute(self):
        """Store communities in graph.

        The function maps values from the partition to graph nodes.
        Nodes found in communities with values larger than n_comms are put into the null community.
        Similarly, nodes not found in main_sub will also be added to the null community.

        Returns
        -------
        name2group: {str: int}
            Maps name of transcript to its community in the partition
        """
        n_comms = min(self.n_comms, len(set(self.partition.membership)))
        main2group = dict(zip(self.main_sub.nodes(), self.partition.membership))
        name2group = {}
        for node in self.graph.nodes():
            group = min(n_comms, main2group.get(node, n_comms))
            name2group[node] = group
        networkx.set_node_attributes(self.graph, name='Group', values=name2group)
        return name2group

    def build(self, clear_adj=True, main_sub=True):
        """Makes networkx graph from adjacency matrix

        Parameters
        ----------
        clear_adj: bool (default=True)
            If False, skip running self.adj=None to clear memory.
        main_sub: bool (default=True)
            If False, skip finding main commected subgraph.
        """
        self.apply_threshold()
        if isinstance(self.adj, np.ndarray):
            self.graph = networkx.from_numpy_array(self.adj)
        else:
            self.graph = networkx.from_pandas_adjacency(self.adj)
        networkx.relabel_nodes(self.graph, lambda n: str(n), False)
        if clear_adj:
            self.adj = None
        if main_sub:
            self.find_main_sub()

    def make_gml_csv_files(self):
        """Wrapper function for most common way to generate .gml file."""
        self.build()
        self.save(True)
        self.get_partition()
        name2group = self.membership2attribute()
        if self.gml_path is not None:  # This is outside of the save function intentionally.
            self.save()
        self.df = pd.DataFrame.from_dict(name2group, orient='index')
        self.df.columns = ['Group']
        if self.csv_path is not None:
            self.df.to_csv(self.csv_path)
        if self.tempfile_path:
            self.tempfile_path.unlink()
