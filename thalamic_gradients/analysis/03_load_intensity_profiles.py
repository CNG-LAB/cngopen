#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 11:50:00 2022

@author: ajohn

- data: micamics processed with micapipe (-MPC)
- load micapipes intensity profiles (Schaefer_parcellation 200) of cortex
- stack all subjects
- remove pial and gm/wm layer as well as medial wall parcels 

"""
import numpy as np

sub_list=np.arange(1,51,1)
allsub_int_profiles=np.empty((14,202,50))


for s, sub in enumerate(sub_list):
    
    path = "/mica-mics/singularity_out/SUB_HC0{:02d}/micapipe/sub-HC0{:02d}/ses-01/anat/surfaces/micro_profiles/".format(sub,sub)
    int_profiles=path+"sub-HC0{:02d}_ses-01_space-fsnative_atlas-schaefer-200_desc-intensity_profiles.txt".format(sub)
    
    # Load the intensity profiles
    mtx_int = np.loadtxt(int_profiles, dtype=np.float, delimiter=' ')
    allsub_int_profiles[:,:,s]=mtx_int
    
# matrix has shape (14,202,50) -> delete first and last row (pial and gm/wm) and  column 0 and 201 (medial wall)
allsub_int_profiles=np.delete(allsub_int_profiles, (0,13), axis=0)
allsub_int_profiles=np.delete(allsub_int_profiles, (0,101), axis=1)  # shape (12,200,50)

np.save("/Project1_thalamus_gradients/structural_covariance/cortex_intensity_profiles/intensity_profiles_schaefer200_allsub.npy", allsub_int_profiles)


# build mean across layers
mean_int=np.mean(allsub_int_profiles, axis=0)

np.save("/Project1_thalamus_gradients/structural_covariance/cortex_intensity_profiles/intensity_profiles_mean_schaefer200_allsub.npy",mean_int)