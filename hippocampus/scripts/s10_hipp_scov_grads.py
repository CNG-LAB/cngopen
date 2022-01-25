"""
computes the hippocampus-to-cortex structural intensity covariance and 
its gradients 
"""
import os
import numpy as np
import h5py
from brainspace.gradient import GradientMaps

# data dirs
ddir = '../data/tout_group'

# get cortical T1w/T2w measures across subjects
h1 = h5py.File(os.path.join(ddir, 'H709_mye_cortex.h5'), 'r')

mye_cortex = np.array(h1['data'])
print(mye_cortex.shape, mye_cortex.min(), mye_cortex.max())
# (360, 709) 0.02919340319931507 4.150848388671875

# get subfield T1w/T2w measures across subjects
h2 = h5py.File(os.path.join(ddir, 'H709_mye_subfields.h5'), 'r')

mye_subfields = np.array(h2['data'])
print(mye_subfields.shape, mye_subfields.min(), mye_subfields.max())
# (4096, 709), 0.8110781, 33.222504

# compute the covariance matrix
mye_scov = np.corrcoef(mye_cortex, mye_subfields)

h2 = h5py.File('../data/tout_group/Hmean709scov_all.h5', 'w')
h2.create_dataset('data', data = mye_scov)
h2.close()

# cut the subfield-to-cortex covariance portion from the matrix
mye_scov = mye_scov[360:, 0:360]
print(mye_scov.shape, mye_scov.min(), mye_scov.max())
# (4096, 360) -0.03890701537805737 0.790964073398187

# get gradients
np.random.seed(0)
gm = GradientMaps(approach = 'dm', kernel='normalized_angle');
gm = gm.fit(mye_scov);   # (4096, 10)

G1 = gm.gradients_[:,0]      #(4096,)
G2 = -1*gm.gradients_[:,1]   #(4096,)
G3 = gm.gradients_[:,2]      #(4096,)

G1_LSUB = G1[0:1024]
G1_LCA = G1[1024:1024+2048]
G1_LDG = G1[1024+2048:1024+2048+1024]

G2_LSUB = G2[0:1024]
G2_LCA = G2[1024:1024+2048]
G2_LDG = G2[1024+2048:1024+2048+1024]

# concatenate subfield gradients and save
data = np.concatenate((G1.reshape(1,-1), G2.reshape(1,-1)))

print(G1.min(), G1.max()) # (-0.14878235398954343, 0.17707568637907664)
print(G2.min(), G2.max()) # (-0.10707012947806552, 0.04844059970187576)


h = h5py.File('../data/tout_group/Hmean709scovGradients_left.h5', 'w')
h.create_dataset('data', data = data)
h.close()