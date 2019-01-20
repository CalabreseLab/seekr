#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""


import numpy as np
import pandas as pd
import igraph
import networkx
import louvain


class Maker:
    """Generate and manipulate network graphs using networkx and igraph.

    Parameters
    ----------
    adj: str | ndarray | DataFrame
        (Path to) adjacency matrix describing weighted edges between nodes.
    gml_path: str
        Path to output file for new network
    csv_path: str
        Path to two column file storing transcript name and membership
    limit: float (default=0)
        Value for thresholding adjacency matrix. Below this limit, all edges are 0.
    gamma: float (default=1)
        Resolution parameter for louvain algorithm
    n_comms: int (default=5)
        Number of communities to find. This does not count a null community.
    seed: int (default=None)
        Set the seed for reproducible results from louvain.

    Attributes
    ----------
    graph: networkx.Graph
        Primary graph build from adjacency matrix
    main_sub: networkx.Graph
        Largest connected subgraph of primary graph
    partition: louvain.RBConfigurationVertexPartition
        Contains the membership data for each node
    df: DataFrame
        Stores transcript name and membership
    """

    def __init__(self, adj=None, gml_path=None, csv_path=None,
                 limit=0, gamma=1, n_comms=5, seed=None):
        self.adj = adj
        if adj is not None:
            self.adj = self.get_adj(adj)
        self.gml_path = gml_path
        self.csv_path = csv_path
        self.limit = limit
        self.gamma = gamma
        self.n_comms = n_comms
        self.seed = seed

        self.graph = None
        self.main_sub = None
        self.partition = None
        self.df = None

    def get_adj(self, adj):
        """Load adjacency matrix from .csv or .npy file, if necessary.

        Parameters
        ----------
        adj: str | ndarray | DataFrame
            (Path to) adjacency matrix describing weighted edges between nodes.

        Returns
        -------
        adj: ndarray | DataFrame
            Adjacency matrix describing weighted edges between nodes.
        """
        adj_types = (str, pd.DataFrame, np.ndarray)
        err_msg = f'adj must be one of {adj_types}, not {type(adj)}.'
        assert type(adj) in adj_types, err_msg
        if isinstance(adj, str):
            try:
                adj = pd.read_csv(adj, index_col=0)
            except UnicodeDecodeError:
                adj = np.load(adj)
        return adj

    def threshold(self):
        """Remove low weighted edges from graph."""
        adj = self.adj
        if isinstance(adj, pd.DataFrame):
            adj = adj.values
        np.fill_diagonal(adj, 0)
        adj[adj < self.limit] = 0

    def find_main_sub(self):
        """Calculate the largest connected subgraph

        Notes
        -----
        https://en.wikipedia.org/wiki/Connected_component_(graph_theory)
        """
        subgraphs = list(networkx.connected_component_subgraphs(self.graph))
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
        networkx.write_gml(graph, self.gml_path)

    def get_partition(self, gml_path=None):
        """Run Louvain to call communities on the graph

        Parameters
        ----------
        gml_path: str (default=None)
            Alternative location for loading in .gml file. Defaults to self.gml_path.
        """
        if gml_path is None:
            gml_path = self.gml_path
        graph = igraph.Graph.Read_GML(gml_path)
        if self.seed is not None:
            louvain.set_rng_seed(self.seed)
        self.partition = louvain.find_partition(graph,
                                                louvain.RBConfigurationVertexPartition,
                                                weights='weight',
                                                resolution_parameter=self.gamma)

    def membership2attribute(self):
        """Store communities in graph.

        The function maps values from the Louvain partition to graph nodes.
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
        self.threshold()
        if isinstance(self.adj, np.ndarray):
            self.graph = networkx.from_numpy_array(self.adj)
        else:
            self.graph = networkx.from_pandas_adjacency(self.adj)
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
        self.save()
        self.df = pd.DataFrame.from_dict(name2group, orient='index')
        self.df.columns = ['Group']
        if self.csv_path is not None:
            self.df.to_csv(self.csv_path)
