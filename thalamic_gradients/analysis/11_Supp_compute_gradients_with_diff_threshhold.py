#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on 25.1.2024

@author: ajohn
script to compute structural connectivity gradients with different threshholds 
#left and right hemisphere separately (still hard coded)

"""
import numpy as np
from brainspace.gradient import GradientMaps
import nibabel as nb


th=0.25
thresh=25

""" left hemisphere """


# import refined thalamus mask as reference
thala_ref_lh_path="/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()

# import structural connectivity matrix
conn_matrix_l=np.load("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/struc_conn_matrix_lh_avg.npy")


# COMPUTE GRADIENT
gm_l = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
gm_l.fit(conn_matrix_l, sparsity=th)



""" right hemisphere """

# import refined thalamus mask as reference
thala_ref_rh_path="/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()

# import structural connectivity matrix
conn_matrix_r=np.load("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/struc_conn_matrix_rh_avg.npy")


# COMPUTE GRADIENT
gm_r = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
gm_r.fit(conn_matrix_r, sparsity=th)


""" project gradient values onto thalamus  """


# iterate over first 2 gradients
for g in range(2):
    
    image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
    idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
    #iterate over indices
    for i, idx in enumerate(idx_l):         
        image_tmp[idx]=gm_l.gradients_[i, g]        # fill in gradients value in thalamus mask     
    
    idx_r=zip(*np.where(thala_ref_rh==1)) 
    #iterate over indices
    for i, idx in enumerate(idx_r):         
        image_tmp[idx]=gm_r.gradients_[i, g]        # fill in gradients value in thalamus mask     
    
    
    # save as nifti, with same header informations as thalamus mask
    tha_mask_l_=nb.load(thala_ref_lh_path)
    clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradient_threshholds/gradient{}_tha_{}.nii.gz'.format(g+1,thresh))

