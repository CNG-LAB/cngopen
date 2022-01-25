"""
computes the group-level gradients (G1, G2, ...) of hippocampal-cortical
functional connectivity (fc), separately for L & R hemispheres.
"""
import os
import h5py
import numpy as np
from brainspace.gradient import GradientMaps

# data dirs
ddir     = '../data/'   
conndir  = '../data/tout_hippoc/' 
odir     = '../data/tout_group/'

# get HCP - S900 subject list        
subjlist = '../data/subjectListS900_QC_gr.txt'
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
subjlist = joinedlist = mylist[:-1]
print('We have now %i subjects... ' % (len(subjlist)))  # 709 subjects

# Left hemisphere
# get hippocampus-to-cortex connectivity across subjects
j = 0
H = np.zeros((len(subjlist), 4096, 360 ));
for subjID in subjlist:   
    # get the hippocampus-to-cortex connectivity for each subject & append
    subjconn_left = os.path.join(conndir, 'HCP_' + subjID + '_left.h5')
    with h5py.File(subjconn_left, "r") as f:        
        subjdata_left = np.array(f['HCP_'+subjID])   # (4096, 360)
        H[j, : , :] = subjdata_left
        j +=1

Hmean = H.mean(axis=0)
print('We averaged hippoca-to-cortex conn. from %i subjects...' % (j))
print('HHHH ', Hmean.shape, np.ndim(Hmean), Hmean.min(), Hmean.max())

# get gradients
np.random.seed(0)
gm = GradientMaps(approach = 'dm', kernel='normalized_angle');
gm = gm.fit(Hmean);   # (4096, 10)

h = h5py.File(os.path.join(odir, 'Hmean709connGradients_left.h5'), 'w')
h.create_dataset('gradients_', data = gm.gradients_ )
h.create_dataset('lambdas_', data = gm.lambdas_ )
h.close()


# Right hemisphere
# get average hippocampus-to-cortex connectivity across subjects
j = 0
H = np.zeros((len(subjlist), 4096, 360 ));
for subjID in subjlist:   
    # get the hippocampus-to-cortex connectivity for each subject & append
    subjconn_left = os.path.join(conndir, 'HCP_' + subjID + '_right.h5')
    with h5py.File(subjconn_left, "r") as f:        
        subjdata_left = np.array(f['HCP_'+subjID])   # (4096, 360)
        H[j, : , :] = subjdata_left
        j +=1

Hmean = H.mean(axis=0)
print('We averaged hippoca-to-cortex conn. from %i subjects...' % (j))
print('HHHH ', Hmean.shape, np.ndim(Hmean), Hmean.min(), Hmean.max())

# get gradients
np.random.seed(0)
gm = GradientMaps(approach = 'dm', kernel='normalized_angle');
gm = gm.fit(Hmean);   # (4096, 10)

h = h5py.File(os.path.join(odir, 'Hmean709connGradients_right.h5'), 'w')
h.create_dataset('gradients_', data = gm.gradients_ )
h.create_dataset('lambdas_', data = gm.lambdas_ )
h.close()
