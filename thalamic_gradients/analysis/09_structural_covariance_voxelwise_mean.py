#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 13:30:24 2022

@author: ajohn
script to create voxelwise structural covariance matrix
1. use mean T1q across layers 
2. correlating between voxelwise thalamic T1q and cortex parcel T1q
3. links this to structural connectivity gradients by correlating columns (parcels) of structural covariance with G1 and G2 gradient loadings 
#left and right hemisphere separately (still hard coded)
change gradient manually in script
"""


import numpy as np
import nibabel as nb


sublist=np.arange(1,51,1)  


'''left thalamus'''
thala_ref_lh=nb.load("/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz").get_fdata()
nr_v=np.count_nonzero(thala_ref_lh)
allsub_T1q=np.zeros((nr_v,50)) 
idx_l=np.where(thala_ref_lh==1)


#iterate over mp2rage in mni space of each subject
for s,sub in enumerate(sublist):
    print(sub)
    #load in mprage_T1map -> 2mm 1mm?
    T1q=nb.load("/data/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC0{:02d}/fnirt_mp2rage_T1map_to_MNI2mm.nii.gz".format(sub)).get_fdata()
    values=T1q[idx_l]     #returns all values at location where thalamus mask is 1
    allsub_T1q[:,s]=values 
    
    
###correlate voxelwise####

#correlate between thalamic T1q voxels and schaefer parcels 
intensity_profiles=np.load("/data/p_02666/Project1_thalamus/structural_covariance/cortex_intensity_profiles/intensity_profiles_mean_schaefer200_allsub.npy")
parcels=np.arange(0,100,1)  #left hemisphere
    
#corr_coefficients -> returns Rvalues (Pearson)
corr_vx=np.zeros((nr_v,100))
vox_nr=np.arange(0,nr_v,1) ##?

#iterate over voxel
for v, vox in enumerate(vox_nr):
    voxel=allsub_T1q[v,:]
    
    for p,parc in enumerate(parcels):
        x=np.corrcoef(voxel,intensity_profiles[p,:])[0,1]
        corr_vx[v,p]=x   
            
            
np.save("/data/p_02666/Project1_thalamus/structural_covariance/voxelwise_struc_cov/struc_cov_meanint_lh.npy",corr_vx)


## correlate each structural covariance column with DTI gradient

#load gradient vector
grad=np.load("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradients_lh.npy")
R_values=np.zeros((100,1))


#iterate over columns
for p,parc in enumerate(parcels):
    R=np.corrcoef(corr_vx[:,p],grad[:,0])[0,1]     #change gradient 
    #collect R values per parcel
    R_values[p]=R

np.save("/data/p_02666/Project1_thalamus/structural_covariance/voxelwise_struc_cov/meanint_corr_grad1_colums_lh.npy", R_values)


###########################################################################################################################################


'''right thalamus'''
thala_ref_rh=nb.load("/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz").get_fdata()
nr_v=np.count_nonzero(thala_ref_rh)
allsub_T1q=np.zeros((nr_v,50)) 
idx_r=np.where(thala_ref_rh==1)

#iterate over mp2rage in mni space of each subject
for s,sub in enumerate(sublist):
    print(sub)
    #load in mprage_T1map -> 2mm 1mm?
    T1q=nb.load("/data/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC0{:02d}/fnirt_mp2rage_T1map_to_MNI2mm.nii.gz".format(sub)).get_fdata()
    values=T1q[idx_r]     #returns all values at location where thalamus mask is 1
    allsub_T1q[:,s]=values 


###correlate voxelwise####

#correlate between gradient T1q voxels and schaefer parcels  -> per layer
intensity_profiles=np.load("/data/p_02666/Project1_thalamus/structural_covariance/cortex_intensity_profiles/intensity_profiles_mean_schaefer200_allsub.npy")
parcels=np.arange(100,200,1)  #right hemisphere
    
#corr_coefficients -> returns Rvalues (Pearson)
corr_vx=np.zeros((nr_v,100))
vox_nr=np.arange(0,nr_v,1) ##?

#iterate over voxel
for v, vox in enumerate(vox_nr):
    voxel=allsub_T1q[v,:]
    
    for p,parc in enumerate(parcels):
        #x=stats.spearmanr(voxel,intensity_profiles[parc,:])[0]  #to try with spearman
        x=np.corrcoef(voxel,intensity_profiles[parc,:])[0,1]  #pearson
        corr_vx[v,p]=x   
            
np.save("/data/p_02666/Project1_thalamus/structural_covariance/voxelwise_struc_cov/struc_cov_meanint_rh.npy",corr_vx)

## correlate each structural covariance column with DTI gradient

#load gradient vector
grad=np.load("/data/p_02666/Project1_thalamus/structural_connectivity/parcels_200/gradients_rh.npy")
R_values=np.zeros((100,1))


for p,parc in enumerate(parcels):
    R=np.corrcoef(corr_vx[:,p],grad[:,0])[0,1]  #change gradient
    #collect R values per parcel
    R_values[p]=R

np.save("/data/p_02666/Project1_thalamus/structural_covariance/voxelwise_struc_cov/meanint_corr_grad1_colums_rh.npy", R_values)

