#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 13:44:25 2023

@author: ajohn

main question: we have differences in our structural covariance patterns between left and right and we want to check whether there is a left-right bias in the files
1. use mean T1q across layers 
2. correlating between left voxelwise thalamic T1q and wholebrain cortex parcel T1q
3. links this to left structural connectivity gradients by correlating columns (parcels) of structural covariance with G1 and G2 gradient loadings 
4. same for right hem
"""
#import modules
import numpy as np
import nibabel as nb

sublist=np.arange(1,51,1)  
"""
1) calculate structural covariance of left thalamus with whole brain
"""

### left thalamus -> whole brain
thala_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz").get_fdata()
nr_v=np.count_nonzero(thala_ref_lh)
allsub_T1q_l=np.zeros((nr_v,50)) 
idx_l=np.where(thala_ref_lh==1)

#iterate over mp2rage in mni space of each subject
for s,sub in enumerate(sublist):
    print(sub)
    #load in mprage_T1map -> 2mm 1mm?
    T1q=nb.load("/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC0{:02d}/fnirt_mp2rage_T1map_to_MNI2mm.nii.gz".format(sub)).get_fdata()
    values=T1q[idx_l]     #returns all values at location where thalamus mask is 1
    allsub_T1q_l[:,s]=values 
    
       
###correlate voxelwise####

#correlate between gradient T1q voxels and schaefer parcels  -> per layer
intensity_profiles=np.load("/Project1_thalamus_gradients/data/structural_covariance/cortex_intensity_profiles/intensity_profiles_mean_schaefer200_allsub.npy")
parcels=np.arange(0,200,1)  #both hemisphere
    
#corr_coefficients -> retunrs Rvalues (Pearson)
corr_vx=np.zeros((nr_v,200))
vox_nr=np.arange(0,nr_v,1) ##?

for v, vox in enumerate(vox_nr):
    voxel=allsub_T1q_l[v,:]
    
    for p,parc in enumerate(parcels):
        x=np.corrcoef(voxel,intensity_profiles[p,:])[0,1]
        corr_vx[v,p]=x   
            
#->scov
#load gradient vector
grad=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_lh.npy")
# 1. Gradient
R_values=np.zeros((200,1))

#iterate over columns
for p,parc in enumerate(parcels):
    R=np.corrcoef(corr_vx[:,p],grad[:,0])[0,1]
    #collect R values per parcel
    R_values[p]=R

np.save("/Project1_thalamus_gradients/data/structural_covariance/QC/R_values_left_th_whole_cortex_g1.npy", R_values)

# 2. Gradient
R_values=np.zeros((200,1))

#iterate over columns
for p,parc in enumerate(parcels):
    R=np.corrcoef(corr_vx[:,p],grad[:,1])[0,1]
    #collect R values per parcel
    R_values[p]=R

np.save("/Project1_thalamus_gradients/data/structural_covariance/QC/R_values_left_th_whole_cortex_g2.npy", R_values)

#####################################################################
### right thalamus -> whole brain
thala_ref_rh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz").get_fdata()
nr_v=np.count_nonzero(thala_ref_rh)
allsub_T1q=np.zeros((nr_v,50)) 
idx_r=np.where(thala_ref_rh==1)

#iterate over mp2rage in mni space of each subject
for s,sub in enumerate(sublist):
    print(sub)
    #load in mprage_T1map -> 2mm 1mm?
    T1q=nb.load("/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC0{:02d}/fnirt_mp2rage_T1map_to_MNI2mm.nii.gz".format(sub)).get_fdata()
    values=T1q[idx_r]     #returns all values at location where thalamus mask is 1
    allsub_T1q[:,s]=values 

###correlate voxelwise####

#correlate between gradient T1q voxels and schaefer parcels  -> per layer
intensity_profiles=np.load("/Project1_thalamus_gradients/data/structural_covariance/cortex_intensity_profiles/intensity_profiles_mean_schaefer200_allsub.npy")
parcels=np.arange(0,200,1)  #left hemisphere
    
#corr_coefficients -> retunrs Rvalues (Pearson)
corr_vx=np.zeros((nr_v,200))
vox_nr=np.arange(0,nr_v,1) ##?

#iterate over voxel
for v, vox in enumerate(vox_nr):
    voxel=allsub_T1q[v,:]
    
    for p,parc in enumerate(parcels):
        x=np.corrcoef(voxel,intensity_profiles[p,:])[0,1]
        corr_vx[v,p]=x   
#->scov

#load gradient vector
grad=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_rh.npy")
# 1. Gradient
R_values=np.zeros((200,1))

#iterate over columns
for p,parc in enumerate(parcels):
    R=np.corrcoef(corr_vx[:,p],grad[:,0])[0,1]
    #collect R values per parcel
    R_values[p]=R       

np.save("/Project1_thalamus_gradients/data/structural_covariance/QC/R_values_right_th_whole_cortex_g1.npy", R_values)

# 2. Gradient
R_values=np.zeros((200,1))

#iterate over columns
for p,parc in enumerate(parcels):
    R=np.corrcoef(corr_vx[:,p],grad[:,1])[0,1]
    #collect R values per parcel
    R_values[p]=R       

np.save("/Project1_thalamus_gradients/data/structural_covariance/QC/R_values_right_th_whole_cortex_g2.npy", R_values)








