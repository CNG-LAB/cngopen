"""
reads-in solar genetic correlation results for subfields-to-cortex scov
"""
import h5py
import numpy as np
import pandas as pd

# data dir
odir = '../data/tout_group/'

# read-in genetic correlations for SUB-cortex
matrix_RG_LSUB = np.zeros((360, 1024))
tot_node_num = 1024
for k in range(1, 360+1):
    fname = '../solar/scov_lsub/scov_lsub_%s/results_gencor.txt' % (k)
    # READ DATA
    df = pd.read_csv(fname,
                     index_col = 0,
                     header = 0)
    df.index.name = 'node'
    # REORDER INDICES
    node_str = []
    for i in range(1, tot_node_num+1):
        node_str.append( 'nodeH_%s_INORM' % str(i))
        
    df_ordered = pd.DataFrame(index = node_str,
                              columns = ['rE', 'rp', 'rG', 'p', 'p2'])
    # MATCH DATA IN NEW DF
    for nodeID in node_str:
        iA = df.index.get_loc(nodeID)
        iB = df_ordered.index.get_loc(nodeID)
        df_ordered.iloc[iB]['rE'] = df.iloc[iA]['rE']
        df_ordered.iloc[iB]['rp'] = df.iloc[iA]['rp']
        df_ordered.iloc[iB]['rG'] = df.iloc[iA]['rG']
        df_ordered.iloc[iB]['p']  = df.iloc[iA]['p']
        df_ordered.iloc[iB]['p2'] = df.iloc[iA]['p2']
    data_vox = np.array(df_ordered['rG'], dtype = 'float')
    matrix_RG_LSUB[(k-1), :] = data_vox
    
# read-in genetic correlations for CA-cortex
matrix_RG_LCA = np.zeros((360, 2048))
tot_node_num = 2048
for k in range(1, 360+1):
    fname = '../solar/scov_lca/scov_lca_%s/results_gencor.txt' % (k)
    # READ DATA
    df = pd.read_csv(fname,
                     index_col = 0,
                     header = 0)
    df.index.name = 'node'
    # REORDER INDICES
    node_str = []
    for i in range(1, tot_node_num+1):
        node_str.append( 'nodeH_%s_INORM' % str(i))
    df_ordered = pd.DataFrame(index = node_str,
                              columns = ['rE', 'rp', 'rG', 'p', 'p2'])
    # MATCH DATA IN NEW DF
    for nodeID in node_str:
        iA = df.index.get_loc(nodeID)
        iB = df_ordered.index.get_loc(nodeID)
        df_ordered.iloc[iB]['rE'] = df.iloc[iA]['rE']
        df_ordered.iloc[iB]['rp'] = df.iloc[iA]['rp']
        df_ordered.iloc[iB]['rG'] = df.iloc[iA]['rG']
        df_ordered.iloc[iB]['p']  = df.iloc[iA]['p']
        df_ordered.iloc[iB]['p2'] = df.iloc[iA]['p2']
    data_vox = np.array(df_ordered['rG'], dtype = 'float')
    matrix_RG_LCA[(k-1), :] = data_vox
    
# read-in genetic correlations for DG-cortex
matrix_RG_LDG = np.zeros((360, 1024))
tot_node_num = 1024
for k in range(1, 360+1):
    fname = '../solar/scov_ldg/scov_ldg_%s/results_gencor.txt' % (k)
    # READ DATA
    df = pd.read_csv(fname,
                     index_col = 0,
                     header = 0,)
    df.index.name = 'node'
    # REORDER INDICES
    node_str = []
    for i in range(1, tot_node_num+1):
        node_str.append( 'nodeH_%s_INORM' % str(i))
    df_ordered = pd.DataFrame(index = node_str,
                              columns = ['rE', 'rp', 'rG', 'p', 'p2'])
    # MATCH DATA IN NEW DF
    for nodeID in node_str:
        iA = df.index.get_loc(nodeID)
        iB = df_ordered.index.get_loc(nodeID)
        df_ordered.iloc[iB]['rE'] = df.iloc[iA]['rE']
        df_ordered.iloc[iB]['rp'] = df.iloc[iA]['rp']
        df_ordered.iloc[iB]['rG'] = df.iloc[iA]['rG']
        df_ordered.iloc[iB]['p']  = df.iloc[iA]['p']
        df_ordered.iloc[iB]['p2'] = df.iloc[iA]['p2']
    data_vox = np.array(df_ordered['rG'], dtype = 'float')
    matrix_RG_LDG[(k-1), :] = data_vox

matrix_subfields = np.concatenate((matrix_RG_LSUB, matrix_RG_LCA, 
                                   matrix_RG_LDG), axis=1)

print(matrix_subfields.shape, np.nanmin(matrix_subfields), np.nanmax(matrix_subfields))
# (360, 4096) -0.5069329 1.0
h = h5py.File('../data/tout_group/Hmean709gen_subfields.h5', 'w')
h.create_dataset('data', data=matrix_subfields)
h.close()




