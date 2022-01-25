"""
reads-in solar heritability results for the functional connectivity,
for left hemisphere and individual subfields
"""
import os
import h5py
import numpy as np
import pandas as pd


# define node strings in a list for 360 cortical parcellations
tot_node_num = 360
node_str = []
for i in range(1, tot_node_num+1):
    node_str.append('node_'+ str(i))
print(len(node_str))
print(node_str[0], '...', node_str[-1])

# read-in subiculum-cortex mean fc heritability (from SOLAR results)
topdir = '../solar/FC_LSUB/'
H_LSUB = np.zeros((1024, 360))
P_LSUB = np.zeros((1024, 360))

for low_i in range(0, 1024):
 
    lowdir = 'fc_' + str(low_i+1)
    fname = os.path.join(topdir, lowdir, 'fc_lsub_results_herit.txt')

    results = pd.read_csv(fname, index_col = 0, header = 0)
    results.index.name = 'node'

    df_results = pd.DataFrame(index = node_str, columns = ['H2r', 'rp'])
    # sorting results from node_1 to node_360
    for nodeID in range(1, tot_node_num+1):

        iA = results.index.get_loc(nodeID)
        iB = df_results.index.get_loc('node_'+ str(nodeID))

        df_results.iloc[iB]['H2r'] = results.iloc[iA]['H2r']
        df_results.iloc[iB]['rp']  = results.iloc[iA]['rp']

    H_LSUB[low_i, :] = np.array(df_results['H2r'], dtype='float64')
    P_LSUB[low_i, :] = np.array(df_results['rp'], dtype='float64')
    

# read-in CA-cortex mean fc heritability (from SOLAR results)
topdir = '../solar/FC_LCA/'
H_LCA = np.zeros((2048, 360))
P_LCA = np.zeros((2048, 360))

for low_i in range(0, 2048):
 
    lowdir = 'fc_' + str(low_i+1)
    fname = os.path.join(topdir, lowdir, 'fc_lca_results_herit.txt')

    results = pd.read_csv(fname, index_col = 0, header = 0)
    results.index.name = 'node'

    df_results = pd.DataFrame(index = node_str, columns = ['H2r', 'rp'])

    for nodeID in range(1, tot_node_num+1):
        
        iA = results.index.get_loc(nodeID)
        iB = df_results.index.get_loc('node_'+ str(nodeID))

        df_results.iloc[iB]['H2r'] = results.iloc[iA]['H2r']
        df_results.iloc[iB]['rp']  = results.iloc[iA]['rp']

    H_LCA[low_i, :] = np.array(df_results['H2r'], dtype='float64')
    P_LCA[low_i, :] = np.array(df_results['rp'], dtype='float64')


# read-in DG-cortex mean fc heritability (from SOLAR results)
topdir = '../solar/FC_LDG/'

H_LDG = np.zeros((1024, 360))
P_LDG = np.zeros((1024, 360))

for low_i in range(0, 1024):
 
    lowdir = 'fc_' + str(low_i+1)
    fname = os.path.join(topdir, lowdir, 'fc_ldg_results_herit.txt')

    results = pd.read_csv(fname, index_col = 0, header = 0)
    results.index.name = 'node'

    df_results = pd.DataFrame(index = node_str, columns = ['H2r', 'rp'])

    for nodeID in range(1, tot_node_num+1):

        iA = results.index.get_loc(nodeID)
        iB = df_results.index.get_loc('node_'+ str(nodeID))

        df_results.iloc[iB]['H2r'] = results.iloc[iA]['H2r']
        df_results.iloc[iB]['rp']  = results.iloc[iA]['rp']

    H_LDG[low_i, :] = np.array(df_results['H2r'], dtype='float64')
    P_LDG[low_i, :] = np.array(df_results['rp'], dtype='float64')

# concatenate all subfield data together and save
H = np.concatenate((H_LSUB, H_LCA, H_LDG), axis=0)
P = np.concatenate((P_LSUB, P_LCA, P_LDG), axis=0)
print(H.shape, P.shape, np.nanmin(H), np.nanmax(H))

h = h5py.File('../data/tout_group/Hmean709_FC_herit_left.h5', 'w')
h.create_dataset('h2r', data = H)
h.create_dataset('p', data = P)
h.close()

