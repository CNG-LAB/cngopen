#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revision Supplementary Figure 5

@author: ajohn 20.08.2024
script to compute correlation between grouplevel T1q map and subject T1q maps
#left and right hemisphere separately (still hard coded)
"""
import numpy as np
import nibabel as nb


""" left hemisphere """
# import individual qT1 maps and reference to extract qT1 values 
T1q_sub_lh=np.load("/Project1_thalamus_gradients/data/structural_covariance/subjectlevel_T1q_values_lh.npy")
thala_ref_lh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz").get_fdata()
idx_l=np.where(thala_ref_lh==1)

# import grouplevel qT1 maps
group_T1q=np.loadtxt("/Project1_thalamus_gradients/data/structural_covariance/grouplevel_T1q_values_lh.txt")

########## correlation between individual qT1_maps and grouplevel qT1_maps

correlation=np.zeros(50)
for s in range(50):
    print(s)
    T1q_sub=T1q_sub_lh[:,:,:,s]
    T1q_sub=T1q_sub[idx_l] 
    r = np.corrcoef(group_T1q,T1q_sub)
    correlation[s] = r[0,1] 


np.save("/Project1_thalamus_gradients/data/structural_covariance/Revision_individual_qT1maps/corr_individual_group_qT1maps_lh.npy",correlation)


""" right hemisphere """

# import individual qT1 maps and reference to extract qT1 values 
T1q_sub_rh=np.load("/Project1_thalamus_gradients/data//data/p_02666/Project1_thalamus_gradients/structural_covariance/subjectlevel_T1q_values_rh.npy")
thala_ref_rh=nb.load("/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz").get_fdata()
idx_r=np.where(thala_ref_rh==1)

# import grouplevel qT1 maps
group_T1q=np.loadtxt("/Project1_thalamus_gradients/data/structural_covariance/grouplevel_T1q_values_rh.txt")

########## correlation between individual qT1_maps and grouplevel qT1_maps

correlation=np.zeros(50)
for s in range(50):
    print(s)
    T1q_sub=T1q_sub_rh[:,:,:,s]
    T1q_sub=T1q_sub[idx_r] 
    r = np.corrcoef(group_T1q,T1q_sub)
    correlation[s] = r[0,1] 


np.save("/Project1_thalamus_gradients/data/structural_covariance/Revision_individual_qT1maps/corr_individual_group_qT1maps_rh.npy",correlation)
