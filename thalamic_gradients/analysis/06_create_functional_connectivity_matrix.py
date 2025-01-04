#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 5.7.2023

@author: ajohn

1. import timeseries data cortex in schaefer 200 parcellation 
2. import thalamus voxel timeseries
3. correlate timeseries to create functional connectivity matrix and fisher zscore rows
4. create grouplevel fc matrix

data set: mica mics
schaefer parcellation: 200
left hemisphere and right hemisphere separately (still hard coded)
"""


#import modules
import numpy as np
import nibabel as nb
import scipy.stats

sublist = np.arange(1,51,1)
schaefer_ts_stack = np.zeros((50,695,200))

#1. import timeseries schaefer parcellation 
for s,sub in enumerate(sublist):
    schaefer = np.loadtxt("/mica-mics/derivatives/micapipe/sub-HC0{:02d}/ses-01/func/sub-HC0{:02d}_ses-01_space-fsnative_atlas-schaefer200_desc-timeseries.txt".format(sub,sub),delimiter=",")
    print(schaefer.shape)

    #remove subcortex and MW
    columns_to_delete=[i for i in range(14)]+[14,115]   #0..13 is subcortex, 14 left MW, 115 right MW
    schaefer_ts=np.delete(schaefer, columns_to_delete,1)
    
    
    #the first 4 subject have a 100 timepoints longer timeseries -> use always the first 695 timepoints 
    if schaefer_ts.shape[0]>695:
        schaefer_ts=schaefer_ts[:695,:]
    #save ts from all subjects in one stack
    schaefer_ts_stack[s,:,:]=schaefer_ts
    

#2. import thalamus timeseries
# we use the refined thalamus mask in MNI space and the preprocessed ts warped to MNI space 
thalamus_ts_l_stack=np.zeros((50,1068,695))
thalamus_ts_r_stack=np.zeros((50,1029,695))

#import refined thalamus mask 
thalamus_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz").get_fdata()
thalamus_ref_rh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz").get_fdata()


for s,sub in enumerate(sublist):
    print(sub)
    #import preprocessed warped ts
    brain_ts=nb.load("/mica-mics/singularity_out_proc_rsfmri_fromjessica/ants_output_to_mni/sub-HC0{:02d}_ses-01_singleecho_clean_in_MNI.nii.gz".format(sub)).get_fdata()
    #get timeseries from thalamus voxels
    thalamus_ts_l=brain_ts[thalamus_ref_lh==1]    #left
    if thalamus_ts_l.shape[0]>695:
        thalamus_ts_l=thalamus_ts_l[:,:695]
    thalamus_ts_r=brain_ts[thalamus_ref_rh==1]    #right
    if thalamus_ts_r.shape[0]>695:
        thalamus_ts_r=thalamus_ts_r[:,:695]
    #save ts from all subjects in one stack
    thalamus_ts_l_stack[s,:,:] =thalamus_ts_l          # left
    thalamus_ts_r_stack[s,:,:] =thalamus_ts_r          # right

np.save("/Project1_thalamus_gradients/data/functional_connectivity/thalamus_timeseries_l_stack.npy", thalamus_ts_l_stack)    
np.save("/Project1_thalamus_gradients/data/functional_connectivity/thalamus_timeseries_r_stack.npy", thalamus_ts_r_stack)  

#3. correlate thalamus voxel ts with schaefer parcel ts
#left
fc_l_stack=np.zeros((50,1068,100))

for s,sub in enumerate(sublist):
    print(sub)
    #iterate over thalamus voxels
    for v in range(1068):
        #iterate over cortex parcels
        for p in range(100):
            coef, p_val = scipy.stats.pearsonr(thalamus_ts_l_stack[s,v,:], schaefer_ts_stack[s,:,p])
            fc_l_stack[s,v,p]=coef

#fisher z-transform each subjects rows                
for s,sub in enumerate(sublist):
    for v in range(1068):
       fc_l_stack[s,v,:]=np.arctanh(fc_l_stack[s,v,:])
    
np.save("/Project1_thalamus_gradients/data/functional_connectivity/fc_l_stack.npy" ,fc_l_stack)

#right
fc_r_stack=np.zeros((50,1029,100))

for s,sub in enumerate(sublist):
    print(sub)
    #iterate over thalamus voxels
    for v in range(1029):
        #iterate over cortex parcels
        for p in range(100):
            p_r=p+100    #because right hemisphere starts from index 100
            coef, p_val = scipy.stats.pearsonr(thalamus_ts_r_stack[s,v,:], schaefer_ts_stack[s,:,p_r])
            fc_r_stack[s,v,p]=coef
            
#z-transform each subjects rows            
for s,sub in enumerate(sublist):
    for v in range(1029):
       fc_r_stack[s,v,:]=np.arctanh(fc_r_stack[s,v,:])

np.save("/Project1_thalamus_gradients/data/functional_connectivity/fc_r_stack.npy" ,fc_r_stack)

#4. compute grouplevel fc matrices
fc_l= np.mean(fc_l_stack, axis=0)
np.save("/Project1_thalamus_gradients/data/functional_connectivity/fc_l.npy" ,fc_l)

fc_r= np.mean(fc_r_stack, axis=0)
np.save("/Project1_thalamus_gradients/data/functional_connectivity/fc_r.npy" ,fc_r)
