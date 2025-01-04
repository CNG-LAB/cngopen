#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 12:41:42 2022

@author: ajohn
script to compute gradients from structural connectivity matrix at group level
#left and right hemisphere separately (hard coded)

!! set treshhold for gradient computation 
(first version we used 0.9 (default), for final version we changed to 0.75 (because structural connectome is sparse per se))

"""
import numpy as np
from brainspace.gradient import GradientMaps
from brainspace.gradient import compute_affinity
import nibabel as nb

# set threshold for computing the gradient
sparsity=0.75   

""" left hemisphere """

# paths to save results
gradients = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_lh.npy"
affinity = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/affinity_lh.npy"
lambdas = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/lambdas_lh.npy"


# import refined thalamus mask as reference
thala_ref_lh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()

# import structural connectivity matrix
conn_matrix_l=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_lh_avg.npy")


# COMPUTE GRADIENT
gm_l = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
gm_l.fit(conn_matrix_l, sparsity=sparsity)
np.save(gradients, gm_l.gradients_)     
np.save(lambdas, gm_l.lambdas_)         


# project gradient values onto thalamus 
# iterate over first 3 gradients
for g in range(3):
    image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
    idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
    #iterate over indices
    for i, idx in enumerate(idx_l):         
        image_tmp[idx]=gm_l.gradients_[i, g]        # fill in gradients value in thalamus mask     
    
    # save as nifti, with same header informations as thalamus mask
    tha_mask_l_=nb.load(thala_ref_lh_path)
    clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)
    
    nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradient{}_left_tha.nii.gz'.format(g+1))


#compute affinity matrix
affinity_matrix = compute_affinity(conn_matrix_l, kernel="normalized_angle", sparsity=sparsity, pre_sparsify=True, non_negative=True, gamma=None)
np.save(affinity,affinity_matrix)


""" right hemisphere """

# paths to save results
gradients = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_rh.npy"
affinity = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/affinity_rh.npy"
lambdas = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/lambdas_rh.npy"


# import refined thalamus mask as reference
thala_ref_rh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()

# import structural connectivity matrix
conn_matrix_r=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_rh_avg.npy")


# COMPUTE GRADIENT
gm_r = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
gm_r.fit(conn_matrix_r, sparsity=sparsity)
np.save(gradients, gm_r.gradients_)
np.save(lambdas, gm_r.lambdas_)

# project gradient values onto thalamus 
# iterate over first 3 gradients
for g in range(3):
    image_tmp=np.zeros(thala_ref_rh.shape)          # create empty image same size as thalamus mask
    idx_r=zip(*np.where(thala_ref_rh==1))           # collect indices of voxels where mask = 1 
    #iterate over indices
    for i, idx in enumerate(idx_r):         
        image_tmp[idx]=gm_r.gradients_[i, g]        # fill in gradients value in thalamus mask     
    
    # save as nifti, with same header informations as thalamus mask
    tha_mask_r_=nb.load(thala_ref_rh_path)
    clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)
    
    nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradient{}_right_tha.nii.gz'.format(g+1))

#compute affinity matrix
affinity_matrix = compute_affinity(conn_matrix_r, kernel="normalized_angle", sparsity=sparsity, pre_sparsify=True, non_negative=True, gamma=None)
np.save(affinity,affinity_matrix)

