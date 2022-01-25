"""
computes individual-level functional connectivity gradients 
by aligning them to the group-level gradients
usage: $ python s02_hipp_fc_grads_indiv.py HCP_165840 LSUB
"""

import os, sys
import h5py
import numpy as np
from brainspace.gradient import GradientMaps

# data dir
ddir     = '../data/tout_hippoc/'   
odir     = '../data/tout_hippoc_grad/' 

# LEFT HIPPOCAMPUS
# read-in group-level fc gradient 
group_gradient_file = '../data/tout_group/Hmean709connGradients_left.h5'
with h5py.File(group_gradient_file, 'r') as g:
    group_gradients = np.array(g['gradients_'])  # (4096, 10)
print('Left: We had computed group-level gradients : ', group_gradients.shape)  

#subjID = 'HCP_165840'
subjID = sys.argv[1]   

# roi = 'LSUB'
roi = sys.argv[2]

# get the hippocampus-to-cortex connectivity for each subject & append
subjconn = os.path.join(ddir, subjID + '_left.h5')

with h5py.File(subjconn, "r") as f:        
    subjdata = np.array(f[subjID])   # (4096, 360)
    print(subjID, subjdata.shape)

if roi =='LSUB':
    subjconn_roi = subjdata[0:1024,:]
    group_gradients_roi = group_gradients[0:1024,:]
elif roi == 'LCA':
    subjconn_roi = subjdata[1024:1024+2048,:]
    group_gradients_roi = group_gradients[1024:1024+2048,:]
elif roi == 'LDG':
    subjconn_roi = subjdata[1024+2048:1024+2048+1024,:]
    group_gradients_roi = group_gradients[1024+2048:1024+2048+1024,:]
print(subjconn_roi.shape, group_gradients_roi.shape)      

# aligning individual gradient to group-level gradient
galign = GradientMaps(kernel = 'normalized_angle', 
                      approach = 'dm', 
                      alignment = 'procrustes')
galign.fit(subjconn_roi, reference = group_gradients_roi)

G1 = galign.gradients_[:,0]  
G2 = galign.gradients_[:,1]  
G3 = galign.gradients_[:,2]  
 

# RIGHT HIPPOCAMPUS  
if roi == 'RSUB' or roi == 'RCA' or roi == 'RDG':
  
    group_gradient_file = '../data/tout_group/Hmean709connGradients_right.h5'
    with h5py.File(group_gradient_file, 'r') as g:
        group_gradients = np.array(g['gradients_'])  # (4096, 10)
    print('Right: We had group-level gradients:', group_gradients.shape)  


    # get the hippocampus-to-cortex connectivity for each subject & append
    subjconn = os.path.join(ddir, subjID + '_right.h5')

    with h5py.File(subjconn, "r") as f:        
        subjdata = np.array(f[subjID])   # (4096, 360)
        print(subjID, subjdata.shape)

    if roi =='RSUB':
        subjconn_roi = subjdata[0:1024,:]
        group_gradients_roi = group_gradients[0:1024,:]
    elif roi == 'RCA':
        subjconn_roi = subjdata[1024:1024+2048,:]
        group_gradients_roi = group_gradients[1024:1024+2048,:]
    elif roi == 'RDG':
        subjconn_roi = subjdata[1024+2048:1024+2048+1024,:]
        group_gradients_roi = group_gradients[1024+2048:1024+2048+1024,:]
    print(subjconn_roi.shape, group_gradients_roi.shape)      

    # aligning individual gradient to group-level gradient
    galign = GradientMaps(kernel = 'normalized_angle', 
                          approach = 'dm', 
                          alignment = 'procrustes')
    galign.fit(subjconn_roi, reference = group_gradients_roi)

    G1 = galign.gradients_[:,0]  
    G2 = galign.gradients_[:,1]  
    G3 = galign.gradients_[:,2]  
     
# save
h = h5py.File(os.path.join(odir, subjID + '_G1_%s.h5' %(roi)), 'w')
h.create_dataset(subjID, data = G1); h.close()

h = h5py.File(os.path.join(odir, subjID + '_G2_%s.h5' %(roi)), 'w')
h.create_dataset(subjID, data = G2); h.close()

h = h5py.File(os.path.join(odir, subjID + '_G3_%s.h5' %(roi)), 'w')
h.create_dataset(subjID, data = G3); h.close()

