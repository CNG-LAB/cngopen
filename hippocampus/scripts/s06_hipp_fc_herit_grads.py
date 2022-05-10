"""
computes the gradients of fcon-heritability for subfields
"""
import os
import h5py
import numpy as np
from brainspace.gradient import GradientMaps

# data dir
ddir = '../data/tout_group'

# read-in heritability of mean hippocampal-cortical mean FC
ffile = h5py.File(os.path.join(ddir, 'Hmean709_FC_herit_left.h5'), 'r')
fc_herit_left = np.array(ffile['h2r'])
# check for nan
fc_herit_left[np.where(np.isnan(fc_herit_left))] = 0
print(fc_herit_left.shape)  # (4096, 360)


# read-in hippocampal-cortical FC gradients
ddir     = '../data/tout_hippoc/'      
group_gradient_file = '../data/tout_group/Hmean709connGradients_left.h5'
with h5py.File(group_gradient_file, 'r') as g:
    group_gradients = np.array(g['gradients_'])  # (4096, 10)
print('Left: We had computed group-level gradients : ', group_gradients.shape)  


# subiculum-cortex mean fc-heritability gradient
fc_herit_left_LSUB   = fc_herit_left[0:1024, :]
group_gradients_LSUB = group_gradients[0:1024,:]
# aligning individual gradient to group-level fcon gradient
np.random.seed(0)
galign = GradientMaps(kernel = 'normalized_angle', 
                      approach = 'dm', 
                      alignment = 'procrustes',
                      )

galign.fit(fc_herit_left_LSUB, reference = group_gradients_LSUB)

G1_LSUB = galign.aligned_[:,0]  
G2_LSUB = galign.aligned_[:,1]  
G3_LSUB = galign.aligned_[:,2]  

lambdas_LSUB = galign.lambdas_

# CA-cortex mean fc-heritability gradient
fc_herit_left_LCA    = fc_herit_left[1024:1024+2048, :]
group_gradients_LCA = group_gradients[1024:1024+2048,:]
np.random.seed(0)
galign = GradientMaps(kernel = 'normalized_angle', 
                      approach = 'dm', 
                      alignment = 'procrustes')
galign.fit(fc_herit_left_LCA, reference = group_gradients_LCA)
G1_LCA = galign.aligned_[:,0]  
G2_LCA = galign.aligned_[:,1]  
G3_LCA = galign.aligned_[:,2]  

lambdas_LCA = galign.lambdas_

# DG-cortex mean fc-heritability gradient
fc_herit_left_LDG    = fc_herit_left[1024+2048:1024+2048+1024, :]
group_gradients_LDG  = group_gradients[1024+2048:1024+2048+1024,:]
np.random.seed(0)
galign = GradientMaps(kernel = 'normalized_angle', 
                      approach = 'dm', 
                      alignment = 'procrustes')

galign.fit(fc_herit_left_LDG, reference = group_gradients_LDG)
G1_LDG = galign.aligned_[:,0]  
G2_LDG = galign.aligned_[:,1]  
G3_LDG = galign.aligned_[:,2]  

lambdas_LDG = galign.lambdas_

# concatenate and save
G1 = np.concatenate((G1_LSUB, G1_LCA, G1_LDG), axis=0)
G2 = np.concatenate((G2_LSUB, G2_LCA, G2_LDG), axis=0)
G3 = np.concatenate((G3_LSUB, G3_LCA, G3_LDG), axis=0)

data = np.concatenate((G1.reshape(-1,1), G2.reshape(-1,1), 
                       G3.reshape(-1,1)), axis=1)
print(data.shape) # (4096, 3)

lambdas = np.concatenate((lambdas_LSUB.reshape(-1,1), 
                          lambdas_LCA.reshape(-1,1),
                          lambdas_LDG.reshape(-1,1)), axis=1)
print(lambdas.shape) # (10, 3)

h = h5py.File('../data/tout_group/Hmean709_FC_herit_gradients_left.h5', 'w')
h.create_dataset('gradients', data = data)
h.create_dataset('lambdas', data = lambdas)
h.close()
