#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:35:43 2022

@author: ajohn

1. import probabilisic tractography data
2. build structural connectivity matrix M x N matrix (M - thalamic voxels, N - 200 Schaefer parcels) 
   - iterate over thalamic voxels per schaefer parcel and put as columns inside matrix (number of samples seeded from that voxel reaching the relevant target mask)

data set: mica mics
schaefer parcellation: 200
left hemisphere and right hemisphere separately (hard coded)
"""

#import modules
import numpy as np
import nibabel as nb

""" left hemisphere """

sublist=np.arange(1,51,1) #subject list
schaefer_parcels_lh=np.arange(1,101,1) # left hemisphere

# paths to save results
save_to="/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_lh_sub.npy"
avg_save_to="/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_lh_avg.npy"

# import refined thalamus mask as reference
thala_ref_lh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()

# collect indices of voxels where mask = 1 
idx_l=np.where(thala_ref_lh==1)        
#counts number of voxels in thalamic_mask              
vox_lh=np.count_nonzero(thala_ref_lh) # 1068 voxels

# create array to stack connectivity matrices of all subjects
conn_matrix_stack=np.zeros((vox_lh,100,50))    

# CREATE CONNECTIVITY MATRIX based on probtrackx_output
# iterate over subjects
for s, sub in enumerate(sublist):
    print("sub_HC0{:02d}".format(sub))
    # iterate over each output file (one nifti per parcel) 
    for i, parcel in enumerate(schaefer_parcels_lh): 
        path="/mica-mics/probtrackx_out/sub_HC0{:02d}/ses-01/dwi/left_hem/seeds_to_space-MNI_atlas-schaefer-200_parcel_{}.nii.gz".format(sub,parcel)
        thala_to_parcel=nb.load(path).get_fdata()       # load output file as array
        values=thala_to_parcel[idx_l]                   # returns all values at location where thalamus mask = 1
        conn_matrix_stack[:,i,s]=values                 # puts the values as colum inside array, stack subjects
np.save(save_to, conn_matrix_stack)
#conn_matrix_stack=np.load(save_to)

# AVERAGE connectivity matrix across subjects
conn_matrix_l = np.mean(conn_matrix_stack, axis=2) 

# NORMALIZE columnwise 
for i in range(100):                                    # iterate over columns ( ~ parcels)
    max_val=np.max(conn_matrix_l[:,i])                  # seach for max per column 
    conn_matrix_l[:,i]=conn_matrix_l[:,i]/max_val       # divide every value in column by max value of this columns (to set it to 1)

np.save(avg_save_to, conn_matrix_l)



""" right hemisphere """

sublist=np.arange(1,51,1) #subject list
schaefer_parcels_rh=np.arange(101,201,1)   # right hemisphere

# paths to save results
save_to="/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_rh_sub.npy"
avg_save_to="/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_rh_avg.npy"

# import refined thalamus mask as reference
thala_ref_rh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()

# collect indices of voxels where mask = 1 
idx_r=np.where(thala_ref_rh==1)        
#counts number of voxels in thalamic_mask              
vox_rh=np.count_nonzero(thala_ref_rh)  # 1029 voxels

# create array to stack connectivity matrices of all subjects
conn_matrix_stack=np.zeros((vox_rh,100,50))    

# CREATE CONNECTIVITY MATRIX based on probtrackx_output
# iterate over subjects
for s, sub in enumerate(sublist):
    print("sub_HC0{:02d}".format(sub))
    # iterate over each output file (one nifti per parcel) 
    for i, parcel in enumerate(schaefer_parcels_rh): 
        path="/mica-mics/probtrackx_out/sub_HC0{:02d}/ses-01/dwi/right_hem/seeds_to_space-MNI_atlas-schaefer-200_parcel_{}.nii.gz".format(sub,parcel)
        thala_to_parcel=nb.load(path).get_fdata()       # load output file as array
        values=thala_to_parcel[idx_r]                   # returns all values at location where thalamus mask = 1
        conn_matrix_stack[:,i,s]=values                 # puts the values as colum inside array, stack subjects
np.save(save_to, conn_matrix_stack)
#conn_matrix_stack=np.load(save_to)


# AVERAGE connectivity matrix across subjects
conn_matrix_r = np.mean(conn_matrix_stack, axis=2) 

# NORMALIZE columnwise 
for i in range(100):                                    # iterate over columns ( ~ parcels)
    max_val=np.max(conn_matrix_r[:,i])                  # seach for max per column 
    conn_matrix_r[:,i]=conn_matrix_r[:,i]/max_val       # divide every value in column by max value of this columns (to set it to 1)

np.save(avg_save_to, conn_matrix_r)
