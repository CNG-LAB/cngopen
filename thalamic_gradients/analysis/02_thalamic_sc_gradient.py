#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 12:41:42 2022

@author: ajohn
script to compute gradients from structural connectivity matrix at group level
bins the gradient into 10 bins
#left and right hemisphere separately (hard coded)

!! set treshhold for gradient computation 
(first version we used 0.9 (default), for final version we changed to 0.75 (because structural connectome is sparse per se))

"""
import numpy as np
from brainspace.gradient import GradientMaps
from brainspace.gradient import compute_affinity
import nibabel as nb
import pandas as pd


""" left hemisphere """

# paths to save results
gradients = "/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradients_lh.npy"
affinity = "/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/affinity_lh.npy"
lambdas = "/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/lambdas_lh.npy"


# import refined thalamus mask as reference
thala_ref_lh_path="/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()

# import structural connectivity matrix
conn_matrix_l=np.load("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/struc_conn_matrix_lh_avg.npy")


sparsity=0.75   # set threshhold for computing the gradient


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
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradient{}_left_tha.nii.gz'.format(g+1))


#compute affinity matrix
affinity_matrix = compute_affinity(conn_matrix_l, kernel="normalized_angle", sparsity=sparsity, pre_sparsify=True, non_negative=True, gamma=None)
np.save(affinity,affinity_matrix)


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
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradient{}_left_tha_{}bins.nii.gz'.format(g+1,b))
    
    # save mean projections per bin 
    mean_projections=np.zeros((b,200))              # 10 bins, 200 parcels
    bins_nr=np.arange(1,b+1,1)
    for i, bin_nr in enumerate(bins_nr):
        idx_mean=np.where(bins.codes==i)                 # return indices of each bin
        mean=np.mean(conn_matrix_l[idx_mean,:], axis=1)  # build mean over these rows
        mean_projections[i,0:100]=mean
    np.save("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/mean_projections_per_bin/gradient{}_left_thala_{}bins_mean.npy".format(g+1,b), mean_projections)




""" right hemisphere """

# paths to save results
gradients = "/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradients_rh.npy"
affinity = "/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/affinity_rh.npy"
lambdas = "/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/lambdas_rh.npy"


# import refined thalamus mask as reference
thala_ref_rh_path="/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()

# import structural connectivity matrix
conn_matrix_r=np.load("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/struc_conn_matrix_rh_avg.npy")


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
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradient{}_right_tha.nii.gz'.format(g+1))

#compute affinity matrix
affinity_matrix = compute_affinity(conn_matrix_r, kernel="normalized_angle", sparsity=sparsity, pre_sparsify=True, non_negative=True, gamma=None)
np.save(affinity,affinity_matrix)


# BINNING the gradient into 10 bins
b=10
for g in range(3):
    bins=pd.qcut(gm_r.gradients_[:,g], b)
    values=bins.codes

    #load thalamus_mask
    tha_mask_r=nb.load(thala_ref_rh_path).get_fdata()
    image_tmp=np.zeros(tha_mask_r.shape)            # create empty image same size as thalamus mask 
    idx_r=zip(*np.where(tha_mask_r==1))             # collect indices of voxels where mask = 1 
    
    #iterate over indices
    for i, idx in enumerate(idx_r):
        image_tmp[idx]=values[i]+1                  #fill in bin label (start with 1 not zero)
    
    tha_mask_r_=nb.load(thala_ref_rh_path)
    clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)
    
    nb.save(clipped_img, '/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradient{}_right_tha_{}bins.nii.gz'.format(g+1,b))
    
    # save mean projections per bin 
    mean_projections=np.zeros((b,200))              # 10 bins, 200 parcels
    bins_nr=np.arange(1,b+1,1)
    for i, bin_nr in enumerate(bins_nr):
        idx_mean=np.where(bins.codes==i)                 # return indices of each bin
        mean=np.mean(conn_matrix_r[idx_mean,:], axis=1)  # build mean over these rows
        mean_projections[i,100:200]=mean
    np.save("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/mean_projections_per_bin/gradient{}_right_thala_{}bins_mean.npy".format(g+1,b), mean_projections)




