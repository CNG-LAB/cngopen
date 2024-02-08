#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:46:29 2023

@author: ajohn
script to compute TC functional connectivity gradients at group level

1. import fc matrix
2. compute gradients
3. project them onto thalamus volume
Additionally
4. split them into 10 bins
5. calculate the mean connectivity per !!!structural!!!! gradient bin 
6. calculate the mean connectivity per !!!functional!!!! gradient bin 
"""

import numpy as np
from brainspace.gradient import GradientMaps
from brainspace.gradient import  compute_affinity 
import nibabel as nb
import pandas as pd

#1. import fc matrix and mask

# import refined thalamus mask as reference
thala_ref_lh_path="/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()

# import functional connectivity matrix
fc_l=np.load("/data/p_02666/Project1_thalamus/functional_connectivity/fc_l.npy")

#2.compute gradients
gm_l = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
gm_l.fit(fc_l)

# save first gradient inversed to have same directionaly as struc conn gradient
gm_l.gradients_[:,0]=gm_l.gradients_[:,0]*(-1)

np.save("/data/p_02666/Project1_thalamus/functional_connectivity/gradients_left.npy", gm_l.gradients_)
np.save("/data/p_02666/Project1_thalamus/functional_connectivity/lambdas_lh.npy", gm_l.lambdas_)


for g in range(3):
# project gradient values onto thalamus 
# iterate over first 3 gradients
    image_tmp=np.zeros(thala_ref_lh.shape)       # create empty image same size as thalamus mask
    idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
    #iterate over indices
    for i, idx in enumerate(idx_l):  

        image_tmp[idx]=gm_l.gradients_[i, g]        # fill in gradients value in thalamus mask     
        # save as nifti, with same header informations as thalamus mask
        tha_mask_l_=nb.load(thala_ref_lh_path)
        clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)
        nb.save(clipped_img, '/data/p_02666/Project1_thalamus/functional_connectivity/gradient{}_left_tha.nii.gz'.format(g+1))

#compute affinity matrix
affinity_matrix = compute_affinity(fc_l, kernel="normalized_angle", sparsity=0.9, pre_sparsify=True, non_negative=True, gamma=None)
np.save("/data/p_02666/Project1_thalamus/functional_connectivity/affinity.npy", affinity_matrix)


# BINNING the gradient into 10 bins
b=10
for g in range(3):
    bins=pd.qcut(gm_l.gradients_[:,g], b)
    values=bins.codes

    #load thalamus_mask
    tha_mask_l=nb.load(thala_ref_lh_path).get_fdata()
    image_tmp=np.zeros(tha_mask_l.shape)            # create empty image same size as thalamus mask 
    idx_l=zip(*np.where(tha_mask_l==1))             # collect indices of voxels where mask = 1 
    
    #iterate over indices
    for i, idx in enumerate(idx_l):
        image_tmp[idx]=values[i]+1                  #fill in bin label (start with 1 not zero)
    
    tha_mask_l_=nb.load(thala_ref_lh_path)
    clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/functional_connectivity/gradient{}_left_tha_{}bins.nii.gz'.format(g+1,b))
    

"""right"""
#1. import fc matrix and mask

# import refined thalamus mask as reference
thala_ref_rh_path="/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()

# import functional connectivity matrix
fc_r=np.load("/data/p_02666/Project1_thalamus/functional_connectivity/fc_r.npy")

#2. COMPUTE GRADIENT
gm_r = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
gm_r.fit(fc_r)

# save first gradient inversed to have same directionaly as struc conn gradient
gm_r.gradients_[:,0]=gm_r.gradients_[:,0]*(-1)

np.save("/data/p_02666/Project1_thalamus/functional_connectivity/gradients_right.npy", gm_r.gradients_)

for g in range(3):
# project gradient values onto thalamus 
# iterate over first 3 gradients
    image_tmp=np.zeros(thala_ref_rh.shape)       # create empty image same size as thalamus mask
    idx_r=zip(*np.where(thala_ref_rh==1))           # collect indices of voxels where mask = 1 
    #iterate over indices
    for i, idx in enumerate(idx_r):         
        image_tmp[idx]=gm_r.gradients_[i, g]        # fill in gradients value in thalamus mask     
        # save as nifti, with same header informations as thalamus mask
        tha_mask_r_=nb.load(thala_ref_rh_path)
        clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)
        nb.save(clipped_img, '/data/p_02666/Project1_thalamus/functional_connectivity/gradient{}_right_tha.nii.gz'.format(g+1))

b=10
for g in range(3):
    bins=pd.qcut(gm_r.gradients_[:,g], b)
    values=bins.codes

    #load thalamus_mask
    tha_mask_r=nb.load(thala_ref_rh_path).get_fdata()
    image_tmp=np.zeros(tha_mask_r.shape)           # create empty image same size as thalamus mask 
    idx_r=zip(*np.where(tha_mask_r==1))             # collect indices of voxels where mask = 1 
    
    #iterate over indices
    for i, idx in enumerate(idx_r):
        image_tmp[idx]=values[i]+1                  #fill in bin label (start with 1 not zero)
    
    tha_mask_r_=nb.load(thala_ref_rh_path)
    clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/functional_connectivity/gradient{}_right_tha_{}bins.nii.gz'.format(g+1,b))
    
