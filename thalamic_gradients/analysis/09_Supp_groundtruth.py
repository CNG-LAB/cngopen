#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 14:53:36 2023
@author: ajohn

#check projections of structural connectivity, functional connectivity and structural covariance from different nuclei
#use thomas atlas and calculate structural connectivity for nuclei where we know the projections
#project then the mean structural connectivity of this nucleus on surface 
#in the Paper supplements we show AV, VLP, MD
# same for fc, scov

"""
import numpy as np
import nibabel as nb


nucleus_label=7         # 2=AV 6=VLP 12=MD  ##### 11=CM  12=MD 8=Pulvinar  4=VA  7=VPL  10=MGN

#import thalamus_mask left and right
thala_ref_lh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()
idx_l=np.where(thala_ref_lh==1)

thala_ref_rh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()
idx_r=np.where(thala_ref_rh==1)

#import thomas atlas 
thomas_matrix_l= nb.load("/Project1_thalamus_gradients/data/Atlas/THOMAS_download_from_zenodo/res_to_2mm/thomas_left_2mm.nii.gz").get_fdata()
thomas_matrix_r= nb.load("/Project1_thalamus_gradients/data/Atlas/THOMAS_download_from_zenodo/res_to_2mm/thomas_right_2mm.nii.gz").get_fdata()
label_l=thomas_matrix_l[idx_l]   #vector same length as mask containing the label numbers 
label_r=thomas_matrix_r[idx_r]   #vector same length as mask containing the label numbers 

#
sc_l=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_lh_avg.npy")
sc_r=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_rh_avg.npy")

#select nucleus
idx_nucleus_l=np.where(label_l==nucleus_label)  
idx_nucleus_r=np.where(label_r==nucleus_label)  

# select voxels of struc conn matrix belonging to this nucleus and average rows 
nucleus_conn_l=np.mean(sc_l[idx_nucleus_l],0)
nucleus_conn_r=np.mean(sc_r[idx_nucleus_r],0)
np.save("/Project1_thalamus_gradients/data/structural_connectivity/projections_thomasnuclei/sc_nucleus_{}_l".format(nucleus_label),nucleus_conn_l)
np.save("/Project1_thalamus_gradients/data/structural_connectivity/projections_thomasnuclei/sc_nucleus_{}_r".format(nucleus_label),nucleus_conn_r)


"""""""""""""""""""""""""""""""""""""""""""""""""""
FUNCTIONAL CONN
proof ground truth of functional connectivity 
#use thomas atlas and calculate functional connectivity for nuclei where we know the projections
#project then the mean functional connectivity of this nucleus on surface 
#in the Paper supplements we show AV, VLP, MD
"""

#import func conn matrix
fc_l=np.load("/Project1_thalamus_gradients/data/functional_connectivity/fc_l.npy")
fc_r=np.load("/Project1_thalamus_gradients/data/functional_connectivity/fc_r.npy")

#select nucleus
idx_nucleus_l=np.where(label_l==nucleus_label)  
idx_nucleus_r=np.where(label_r==nucleus_label)  

# select voxel of func conn matrix belonging to this nucleus and average rows 
nucleus_conn_l=np.mean(fc_l[idx_nucleus_l],0)
nucleus_conn_r=np.mean(fc_r[idx_nucleus_r],0)
np.save("/Project1_thalamus_gradients/data/functional_connectivity/projections_thomasnuclei/fc_nucleus_{}_l".format(nucleus_label),nucleus_conn_l)
np.save("/Project1_thalamus_gradients/data/functional_connectivity/projections_thomasnuclei/fc_nucleus_{}_r".format(nucleus_label),nucleus_conn_r)



"""""""""""""""""""""""""""""""""""""""""""""""""""
STRUCTURAL COVARIANCE
proof ground truth of structural covariance
#use thomas atlas and calculate structural covariance for nuclei where we know the projections
#project then the mean structural covariance of this nucleus on surface 
#in the Paper supplements we show AV, VLP, MD
"""

#import func conn matrix
scov_l=np.load("/Project1_thalamus_gradients/data/structural_covariance/voxelwise_struc_cov/struc_cov_lh.npy")
scov_r=np.load("/Project1_thalamus_gradients/data/structural_covariance/voxelwise_struc_cov/struc_cov_rh.npy")

#select nucleus
idx_nucleus_l=np.where(label_l==nucleus_label)  
idx_nucleus_r=np.where(label_r==nucleus_label)  

# select voxel of scov matrix belonging to this nucleus and average rows 
nucleus_conn_l=np.mean(scov_l[idx_nucleus_l],0)
nucleus_conn_r=np.mean(scov_r[idx_nucleus_r],0)
np.save("/Project1_thalamus_gradients/data/structural_covariance/projections_thomasnuclei/scov_l_nucleus_{}".format(nucleus_label),nucleus_conn_l)
np.save("/Project1_thalamus_gradients/data/structural_covariance/projections_thomasnuclei/scov_r_nucleus_{}".format(nucleus_label),nucleus_conn_r)


