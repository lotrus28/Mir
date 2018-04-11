# py 2.7
import copy
import networkx as nx
import numpy as np
import pandas as pd

from collections import deque

def add_symmetry(df):

    for i in df.index:
        for j in df.columns:
            if df.loc[i,j]!=0:
                df.loc[j,i] = -1*(df.loc[i,j])

    return(df)

def from_pd_df(df):

    gr = nx.DiGraph()
    for i in df.index:
        for j in df.columns:
            if df.loc[i, j] != 0:
                gr.add_edge(i, j, weight=df.loc[i, j])
    return(gr)

def bellman_ford(G, source, weight='weight'):
    """Compute shortest path lengths and predecessors on shortest paths
    in weighted graphs.

    The algorithm has a running time of O(mn) where n is the number of
    nodes and m is the number of edges.  It is slower than Dijkstra but
    can handle negative edge weights.

    Parameters
    ----------
    G : NetworkX graph
       The algorithm works for all types of graphs, including directed
       graphs and multigraphs.

    source: node label
       Starting node for path

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    Returns
    -------
    pred, dist : dictionaries
       Returns two dictionaries keyed by node to predecessor in the
       path and to the distance from the source respectively.

    Raises
    ------
    NetworkXUnbounded
       If the (di)graph contains a negative cost (di)cycle, the
       algorithm raises an exception to indicate the presence of the
       negative cost (di)cycle.  Note: any negative weight edge in an
       undirected graph is a negative cost cycle.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(5, create_using = nx.DiGraph())
    >>> pred, dist = nx.bellman_ford(G, 0)
    >>> sorted(pred.items())
    [(0, None), (1, 0), (2, 1), (3, 2), (4, 3)]
    >>> sorted(dist.items())
    [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]

    >>> from nose.tools import assert_raises
    >>> G = nx.cycle_graph(5, create_using = nx.DiGraph())
    >>> G[1][2]['weight'] = -7
    >>> assert_raises(nx.NetworkXUnbounded, nx.bellman_ford, G, 0)

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    The dictionaries returned only have keys for nodes reachable from
    the source.

    In the case where the (di)graph is not connected, if a component
    not containing the source contains a negative cost (di)cycle, it
    will not be detected.

    """
    if source not in G:
        raise KeyError("Node %s is not found in the graph" % source)

    for u, v, attr in G.selfloop_edges(data=True):
        if attr.get(weight, 1) < 0:
            raise nx.NetworkXUnbounded("Negative cost cycle detected.")

    dist = {source: 0}
    pred = {source: None}

    if len(G) == 1:
        return pred, dist

    return _bellman_ford_relaxation(G, pred, dist, [source], weight)

def _bellman_ford_relaxation(G, pred, dist, source, weight):
    """Relaxation loop for Bellman–Ford algorithm

    Parameters
    ----------
    G : NetworkX graph

    pred: dict
        Keyed by node to predecessor in the path

    dist: dict
        Keyed by node to the distance from the source

    source: list
        List of source nodes

    weight: string
       Edge data key corresponding to the edge weight

    Returns
    -------
    Returns two dictionaries keyed by node to predecessor in the
       path and to the distance from the source respectively.

    Raises
    ------
    NetworkXUnbounded
       If the (di)graph contains a negative cost (di)cycle, the
       algorithm raises an exception to indicate the presence of the
       negative cost (di)cycle.  Note: any negative weight edge in an
       undirected graph is a negative cost cycle
    """
    if G.is_multigraph():
        def get_weight(edge_dict):
            return min(eattr.get(weight, 1) for eattr in edge_dict.values())
    else:
        def get_weight(edge_dict):
            return edge_dict.get(weight, 1)

    G_succ = G.succ if G.is_directed() else G.adj
    inf = float('inf')
    n = len(G)

    count = {}
    q = deque(source)
    in_q = set(source)
    while q:
        u = q.popleft()
        in_q.remove(u)
        # Skip relaxations if the predecessor of u is in the queue.
        if pred[u] not in in_q:
            dist_u = dist[u]
            for v, e in G_succ[u].items():
                dist_v = dist_u + get_weight(e)
                if dist_v < dist.get(v, inf):
                    if v not in in_q:
                        q.append(v)
                        in_q.add(v)
                        count_v = count.get(v, 0) + 1
                        if count_v == n:
                            raise nx.NetworkXUnbounded(
                                "Negative cost cycle detected.")
                        count[v] = count_v
                    dist[v] = dist_v
                    pred[v] = u

    return pred, dist

def fill_zeros(df):

    temp = copy.copy(df)

    # Deal only with positive energy reactions
    temp[temp<0] = 0
    # Maximisation => MInimisation
    temp = -1*temp

    gr = from_pd_df(temp)

    for i in temp.index:

        try:
            p = bellman_ford(gr, i, weight = 'weight')

        except nx.NetworkXUnbounded:
            print('{} is involved in an infinite cycle'.format(i))

            return(None)

        for j in p[1]:
            if df.loc[i,j] == 0:
                # Add a small fee for direct reaction to aid symbiotic chains
                temp.loc[i,j] = p[1][j] + 0.05
                temp.loc[j, i] = -1*p[1][j]

    temp = -1*temp

    temp = add_symmetry(temp)

    return(temp)

names = [_ for _ in 'ABCDEFGHJK']

en_df = pd.DataFrame(np.zeros((10,10)), index=names, columns=names)

# cycle1
en_df.loc['A','F'] = 1.9
en_df.loc['F','B'] = 0.8
en_df.loc['B', 'G'] = 4.6
en_df.loc['G','C'] = 1.5
# Direct a -> c should be no more than indirect conversion
en_df.loc['A','C'] = 8.7

# cycle2
en_df.loc['H','D'] = 5
en_df.loc['H','E'] = 2
en_df.loc['E','J'] = 1.1

en_df.loc['K','E'] = 1.2

en_df.loc['H', 'G'] = 1

en_df = add_symmetry(en_df)
x = fill_zeros(en_df)

