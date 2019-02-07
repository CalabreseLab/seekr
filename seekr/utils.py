
import numpy as np
import pandas as pd


def get_adj(adj):
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
