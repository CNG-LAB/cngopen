"""
computes the gradients of hippocampus-to-cortex genetic correlations
"""
import os
import h5py
import numpy as np
from brainspace.gradient import GradientMaps

ddir = '../data/tout_group/'

# read-in hippocampus-to-cortex genetic correlations
k = h5py.File(os.path.join(ddir, 'Hmean709gen_subfields.h5'), 'r')
gen_corr = np.array(k['data']).T 
gen_corr[np.where(np.isnan(gen_corr))]=0

print(gen_corr.shape, np.nanmin(gen_corr), np.nanmax(gen_corr))
# (4096, 360) -0.5069329 1.0

# get the gradients
np.random.seed(0)
gm = GradientMaps(approach = 'dm', kernel='normalized_angle');
gm = gm.fit(gen_corr)

G1 = -1*gm.gradients_[:,0]
G2 = -1*gm.gradients_[:,1]

print(G1.min(), G1.max()) # -0.12774627084212994 0.10906496096980112
print(G2.min(), G2.max()) # -0.08851499714954007 0.07428325700001913

G1_gen_LSUB = G1[0:1024,]
G2_gen_LSUB = G2[0:1024,]

G1_gen_LCA = G1[1024:1024+2048,]
G2_gen_LCA = G2[1024:1024+2048,]

G1_gen_LDG = G1[1024+2048:1024+2048+1024,]
G2_gen_LDG = G2[1024+2048:1024+2048+1024,]

# concatenate subfield gradients and save
data = np.concatenate((G1.reshape(1,-1), G2.reshape(1,-1)))

h = h5py.File(os.path.join(ddir, 'Hmean709genGradients_left.h5'), 'w')
h.create_dataset('data', data = data)
h.close()
