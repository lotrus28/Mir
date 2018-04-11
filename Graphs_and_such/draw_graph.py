# This script draws a circular graph from a dematrx.txt

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def show_graph_with_labels(adjacency_matrix, mylabels):
    rows, cols = np.where(adjacency_matrix > 0)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.DiGraph()
    gr.add_edges_from(edges)
    nx.draw(gr, pos = nx.shell_layout(gr),  node_size=500, labels=mylabels, with_labels=True)
    plt.show()

x = pd.read_table('dematrix.txt', sep ='\t', header = None, index_col = None)

A = x.as_matrix()

labels = {x:x for x in range(A.shape[0])}

show_graph_with_labels(A, labels)