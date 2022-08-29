#!/usr/bin/env python3
"""
Leiden clustering script used in manuscript 
NEUROTRANSMITTER TRANSPORTER/RECEPTOR CO-EXPRESSION SHARES ORGANIZATIONAL TRAITS WITH BRAIN STRUCTURE AND FUNCTION
https://doi.org/10.1101/2022.08.26.505274

Negative-assymetric Leiden clustering is implemented using the leidenalg python package (https://github.com/vtraag/leidenalg), adapting the approach suggested
in the documentation (https://readthedocs.org/projects/leidenalg/downloads/pdf/latest/)
"""
import networkx as nx
import igraph as ig
import numpy as np
from scipy.stats import zscore
import pandas as pd
import leidenalg as ld
from sklearn.metrics.cluster import adjusted_rand_score


input_path = 'path/to/data/'
parcels=100

cortex=pd.read_csv(input_path + '{}Parcels7Networks_receptorprofiles.csv'.format(parcels), index_col=0)
cortex=cortex.apply(zscore)

subcortex = pd.read_csv(input_path + 'tian_subcortex.csv', index_col=0)
subcortex = subcortex.apply(zscore)

def cluster(inp, res, n_reps):
    coma = inp.transpose().corr('spearman').values
    G_nx = nx.Graph(coma)
    G = ig.Graph.from_networkx(G_nx)
    matrix = np.zeros(shape=(n_reps, len(coma)))
    for i in range(n_reps):
        G_pos = G.subgraph_edges(G.es.select(weight_gt=0), delete_vertices=False)
        G_neg = G.subgraph_edges(G.es.select(weight_lt=0), delete_vertices=False)
        G_neg.es['weight'] = [-w for w in G_neg.es['weight']]
        optimiser = ld.Optimiser()
        optimiser.consider_comms = ld.ALL_COMMS
        part_pos = ld.RBConfigurationVertexPartition(G_pos, weights='weight', resolution_parameter=res)
        part_neg = ld.RBConfigurationVertexPartition(G_neg, weights='weight', resolution_parameter=res)
        diff = optimiser.optimise_partition_multiplex([part_pos, part_neg], layer_weights=[1, -1])
        assigns = np.array(part_pos.membership).astype('float')
        assigns += 1
        matrix[i] += assigns
    stability = np.zeros((n_reps, n_reps))
    for i in range(n_reps):
        for it in range(i, n_reps):
            stability[i, it] += adjusted_rand_score(matrix[i], matrix[it])
    mirrored = stability + stability.T
    np.fill_diagonal(mirrored, 1)
    means = np.mean(mirrored, axis=0)
    # get index and return
    best_partition = matrix[list(means).index(max(means))]
    d = {'gamma': res, 'partition': best_partition, 'mean rand': max(means),
         ' mean variance': np.mean(np.var(mirrored, axis=1))}

    np.save(input_path + 'partition/{}_partition_{}.npy'.format(len(coma), res), d)







