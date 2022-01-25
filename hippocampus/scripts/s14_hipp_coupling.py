"""
computes the subfield coupling maps
"""
import os
import h5py
import numpy as np
from scipy.stats import pearsonr

ddir = '../data/tout_group/'

h1  = h5py.File(os.path.join(ddir, 'H709_mye_cortex.h5'), 'r');
mye_cortex = np.array(h1['data'])

h2  = h5py.File(os.path.join(ddir, 'H709_mye_subfields.h5'), 'r');
mye_subfields = np.array(h2['data'])

mye_LSUB = mye_subfields[0:1024,:]
mye_LCA  = mye_subfields[1024:1024+2048,:]
mye_LDG  = mye_subfields[1024+2048:1024+2048+1024,:]

scov_LSUB = np.corrcoef(mye_cortex, mye_LSUB)[360:,0:360]
scov_LCA  = np.corrcoef(mye_cortex, mye_LCA)[360:,0:360]
scov_LDG  = np.corrcoef(mye_cortex, mye_LDG)[360:,0:360]


h3 = h5py.File(os.path.join(ddir, 'Hmean709_FC_left.h5'), 'r')
fcon = np.array(h3['data'])

fcon_LSUB = fcon[0:1024,:]
fcon_LCA  = fcon[1024:1024+2048,:]
fcon_LDG  = fcon[1024+2048:1024+2048+1024,:]


r_lsub = np.zeros((1024,))
r_lca  = np.zeros((2048,))
r_ldg  = np.zeros((1024,))

i = 0
for i in range(0,1024):
    r_lsub[i] = pearsonr(fcon_LSUB[i,:], scov_LSUB[i,:])[0]

j = 0
for j in range(0,2048):
    r_lca[j] = pearsonr(fcon_LCA[j,:], scov_LCA[j,:])[0]

k = 0
for k in range(0,1024):
    r_ldg[k] = pearsonr(fcon_LDG[k], scov_LDG[k,:])[0]
    
    
data = np.concatenate((r_lsub, r_lca, r_ldg))

h = h5py.File('../data/tout_group/Hmean709_coupling_left.h5', 'w')
h.create_dataset('data', data = data)
h.close()
