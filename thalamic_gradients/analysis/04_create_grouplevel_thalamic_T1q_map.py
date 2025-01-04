#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 13:11:48 2023

@author: ajohn
script to create grouplevel qT1 map for left and right hemisphere

"""

import numpy as np
import nibabel as nb


### left hemisphere

thala_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz").get_fdata()
idx_l=np.where(thala_ref_lh==1)

#calculate grouplevel T1q map
T1q_stack=np.zeros((91, 109, 91,50))
sublist = np.arange(1,51,1)            # subjectlist
for s,sub in enumerate(sublist):
    # load in mprage_T1map
    T1q_img = nb.load("/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC0{:02d}/fnirt_mp2rage_T1map_to_MNI2mm.nii.gz".format(sub)).get_fdata()
    T1q_stack[:,:,:,s]=T1q_img
T1q_mean=np.mean(T1q_stack, axis=3)
T1q=T1q_mean[idx_l] 

#save as txt
np.save("/Project1_thalamus_gradients/data/structural_covariance/subjectlevel_T1q_values_lh.npy",T1q_stack)

np.savetxt("/Project1_thalamus_gradients/data/structural_covariance/grouplevel_T1q_values_lh.txt",T1q, delimiter=' ')

#save as nifti
image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
#iterate over indices
for i, idx in enumerate(idx_l):         
    image_tmp[idx]=T1q[i]       # fill in T1q in thalamus mask     
    
# save as nifti, with same header informations as thalamus mask
thala_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz")
clipped_img = nb.Nifti1Image(image_tmp, thala_ref_lh.affine, thala_ref_lh.header)
    
nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_covariance/grouplevel_T1q_map_lh.nii.gz')


### right hemispherre

thala_ref_rh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz").get_fdata()
idx_r=np.where(thala_ref_rh==1)

#calculate grouplevel T1q map
T1q_stack=np.zeros((91, 109, 91,50))
sublist = np.arange(1,51,1)            # subjectlist
for s,sub in enumerate(sublist):
    # load in mprage_T1map
    T1q_img = nb.load("/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC0{:02d}/fnirt_mp2rage_T1map_to_MNI2mm.nii.gz".format(sub)).get_fdata()
    T1q_stack[:,:,:,s]=T1q_img
T1q_mean=np.mean(T1q_stack, axis=3)
T1q=T1q_mean[idx_r] 

np.save("/Project1_thalamus_gradients/data/structural_covariance/subjectlevel_T1q_values_rh.npy",T1q_stack)
#save as txt
np.savetxt("/Project1_thalamus_gradients/data/structural_covariance/grouplevel_T1q_values_rh.txt",T1q, delimiter=' ')

#save as nifti
image_tmp=np.zeros(thala_ref_rh.shape)          # create empty image same size as thalamus mask
idx_r=zip(*np.where(thala_ref_rh==1))           # collect indices of voxels where mask = 1 
#iterate over indices
for i, idx in enumerate(idx_r):         
    image_tmp[idx]=T1q[i]       # fill in gradients value in thalamus mask     
    
# save as nifti, with same header informations as thalamus mask
thala_ref_rh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz")
clipped_img = nb.Nifti1Image(image_tmp, thala_ref_rh.affine, thala_ref_rh.header)
    
nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_covariance/grouplevel_T1q_map_rh.nii.gz')









### for shemata in fig 2 also save some single subjects

#sub0
T1q_sub=T1q_stack[:,:,:,0]
thala_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz").get_fdata()
idx_l=np.where(thala_ref_lh==1)
T1q_sub=T1q_sub[idx_l] 

#save as nifti
idx_l=zip(*np.where(thala_ref_lh==1)) 
image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
#iterate over indices
for i, idx in enumerate(idx_l):         
    image_tmp[idx]=T1q_sub[i]       # fill in T1q in thalamus mask     
# save as nifti, with same header informations as thalamus mask
thala_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz")
clipped_img = nb.Nifti1Image(image_tmp, thala_ref_lh.affine, thala_ref_lh.header)
nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_covariance/T1q_map_lh_sub01.nii.gz')