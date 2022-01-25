"""
computed the fused gradients of hippocampal subfields.
fused matrix is derived from fcon and scov measures
"""
import os
import h5py
import numpy as np
from brainspace.gradient import GradientMaps

gdir = '../data/tout_group/'

# get hippocampus-to-cortex functional connectivity matrix
fcfile = os.path.join(gdir, 'Hmean709_FC_left.h5')
f = h5py.File(fcfile, 'r')
fcon = np.array(f['data'])

print(fcon.shape, fcon.min(), fcon.max())
# (4096, 360) -0.005300521852874321 0.39153784016161197

# get hippocampus-to-cortex structural covariance matrix
scfile = os.path.join(gdir, 'Hmean709scov_all.h5')
s = h5py.File(scfile, 'r')
scov = np.array(s['data'])

scov = scov[360:, 0:360]
print(scov.shape, scov.min(), scov.max())
# (4096, 360) -0.03890701537805737 0.790964073398187

# define fusion function from the **Brainspace
def fusion(*args):
    from scipy.stats import rankdata
    from sklearn.preprocessing import minmax_scale

    max_rk = [None] * len(args)
    masks = [None] * len(args)
    for j, a in enumerate(args):
        m = masks[j] = a != 0
        a[m] = rankdata(a[m])
        max_rk[j] = a[m].max()

    max_rk = min(max_rk)
    for j, a in enumerate(args):
        m = masks[j]
        a[m] = minmax_scale(a[m], feature_range=(1, max_rk))

    return np.hstack(args)

# get the fused matrix
fcon[fcon < 0] = 0
scov[scov < 0] = 0
Hfused = fusion(fcon, scov)
print(Hfused.shape, Hfused.min(), Hfused.max())
# (4096, 720) 0.0 1473890.0

# compute the fused gradients
np.random.seed(0)
Hfused_gm = GradientMaps(n_components=3, approach = 'dm', 
                            kernel='normalized_angle')
Hfused_gm.fit(Hfused)

G1_fused = -1 * Hfused_gm.gradients_[:,0]
G2_fused = Hfused_gm.gradients_[:,1]

print(G1_fused.min(), G1_fused.max()) # -0.20110914442086306 0.1757402036390
print(G2_fused.min(), G2_fused.max()) # -0.08698611405029918 0.1112287977782
data = np.concatenate((G1_fused.reshape(1,-1), G2_fused.reshape(1,-1))).T

h = h5py.File('../data/tout_group/Hmean709fusedGradients_left.h5', 'w')
h.create_dataset('data', data = data)
h.close()
