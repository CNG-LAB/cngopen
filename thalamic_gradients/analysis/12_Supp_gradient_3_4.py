#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 15:43:14 2024

@author: ajohn
REVISION R2 Q2
save gradient 3 and 4 as nifti

#left and right hemisphere separately (hard coded)
"""
import numpy as np
import nibabel as nb


grad=3 # which gradient (as absolute number not index)
g=grad-1

""" left hemisphere """

# import refined thalamus mask as reference
thala_ref_lh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()

gradients=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_lh.npy")    

# project gradient values onto thalamus 

image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
#iterate over indices
for i, idx in enumerate(idx_l):         
    image_tmp[idx]=gradients[i, g]        # fill in gradients value in thalamus mask     

# save as nifti, with same header informations as thalamus mask
tha_mask_l_=nb.load(thala_ref_lh_path)
clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)

nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradient{}_left_tha.nii.gz'.format(g+1))



""" right hemisphere """

gradients_r=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_rh.npy")   

# import refined thalamus mask as reference
thala_ref_rh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()


image_tmp=np.zeros(thala_ref_rh.shape)          # create empty image same size as thalamus mask
idx_r=zip(*np.where(thala_ref_rh==1))           # collect indices of voxels where mask = 1 
#iterate over indices
for i, idx in enumerate(idx_r):         
    image_tmp[idx]=gradients_r[i, g]        # fill in gradients value in thalamus mask     

# save as nifti, with same header informations as thalamus mask
tha_mask_r_=nb.load(thala_ref_rh_path)
clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)

nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradient{}_right_tha.nii.gz'.format(g+1))



'''both hemispheres together'''

image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
#iterate over indices
for i, idx in enumerate(idx_l):         
    image_tmp[idx]=gradients[i, g]        # fill in gradients value in thalamus mask     

idx_r=zip(*np.where(thala_ref_rh==1)) 
#iterate over indices
for i, idx in enumerate(idx_r):         
    image_tmp[idx]=gradients_r[i, g]        # fill in gradients value in thalamus mask     


# save as nifti, with same header informations as thalamus mask
tha_mask_l_=nb.load(thala_ref_lh_path)
clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)

nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradient{}_tha.nii.gz'.format(g+1))

