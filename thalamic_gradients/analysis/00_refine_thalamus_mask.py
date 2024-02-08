#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:14:04 2022
- script refines thalamus mask by manually thresholding

@author: ajohn
"""

import nibabel as nb
import numpy as np


""" left hemisphere """

img=nb.load("/data/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC001/fnirt_uni_T1map_to_MNI2mm_iout.nii.gz").get_fdata()
thala=nb.load("/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh.nii.gz").get_fdata()


idx_l=zip(*np.where(thala==1))

for i, idx in enumerate(idx_l):
    if img[idx] < 1800000:       #manuell set threshhold
        thala[idx]=0
    

img=nb.load("/data/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC001/fnirt_uni_T1map_to_MNI2mm_iout.nii.gz")
clipped_img = nb.Nifti1Image(thala, img.affine, img.header)
nb.save(clipped_img, '/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz')





""" right hemisphere """

img=nb.load("/data/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC001/fnirt_uni_T1map_to_MNI2mm_iout.nii.gz").get_fdata()
thala=nb.load("/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh.nii.gz").get_fdata()


idx_l=zip(*np.where(thala==1))

for i, idx in enumerate(idx_l):
    if img[idx] < 1800000:       #manuell set threshhold
        thala[idx]=0
    

img=nb.load("/data/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_HC001/fnirt_uni_T1map_to_MNI2mm_iout.nii.gz")
clipped_img = nb.Nifti1Image(thala, img.affine, img.header)
nb.save(clipped_img, '/data/p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz')