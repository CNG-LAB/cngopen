#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revision Supplementary Figure 5
@author: ajohn 20.08.2024

script to compute gradient from structural connectivity matrix on subjectlevel
#left and right hemisphere separately (still hard coded)
#we use procrustes alignment to align the individual gradients to the grouplevel gradient (alignment = True)
#correlate individual level to grouplevel gradient -> save r values
"""
import numpy as np
from brainspace.gradient import GradientMaps
import nibabel as nb


alignment=True  #alignment of individual gradients to grouplevel using procrustes

""" left hemisphere """

# import refined thalamus mask as reference
thala_ref_lh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh_refined.nii.gz"
thala_ref_lh=nb.load(thala_ref_lh_path).get_fdata()
group_grad=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_lh.npy")

# import structural connectivity matrix stack (50 subjects)
conn_matrix_sub=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_lh_sub.npy")
sub_grad=np.zeros((1068,10,50))
sub_grad_aligned=np.zeros((1068,10,50))

for i in np.arange(0,50):
    
    sub=i+1
    # paths to save results
    gradients = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_gradients_lh.npy".format(sub)
    affinity = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_affinity_lh.npy".format(sub)
    lambdas = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_lambdas_lh.npy".format(sub)
    
    sparsity=0.75   # set threshhold for computing the gradient
    
    if alignment == False:
        
        # COMPUTE GRADIENT
        gm_l = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
        gm_l.fit(conn_matrix_sub[:,:,i], sparsity=sparsity)
        np.save(gradients, gm_l.gradients_) 
        np.save(lambdas, gm_l.lambdas_)
        sub_grad[:,:,i]=gm_l.gradients_
        
        # project gradient values onto thalamus 
        # iterate over first 2 gradients
        for g in range(2):
            image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
            idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
            #iterate over indices
            for i, idx in enumerate(idx_l):         
                image_tmp[idx]=gm_l.gradients_[i, g]        # fill in gradients value in thalamus mask     
            
            # save as nifti, with same header informations as thalamus mask
            tha_mask_l_=nb.load(thala_ref_lh_path)
            clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)
            
            nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_gradient{}_left_tha_.nii.gz'.format(sub,g+1))
    
    elif alignment == True:        
    # use procrustes alignment to align individual gradients to the grouplevel 

        align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                           kernel='normalized_angle', alignment='procrustes')  
        #align all individual gradients to grouplevel
        align.fit(conn_matrix_sub[:,:,i],reference=group_grad, sparsity=0.75)
        grad_aligned = align.aligned_
        sub_grad_aligned[:,:,i]=grad_aligned
        
        # project gradient values onto thalamus 
        # iterate over first 2 gradients        
        for g in range(2):
            image_tmp=np.zeros(thala_ref_lh.shape)          # create empty image same size as thalamus mask
            idx_l=zip(*np.where(thala_ref_lh==1))           # collect indices of voxels where mask = 1 
            #iterate over indices
            for i, idx in enumerate(idx_l):         
                image_tmp[idx]=align.aligned_[i, g]        # fill in gradients value in thalamus mask     
            
            # save as nifti, with same header informations as thalamus mask
            tha_mask_l_=nb.load(thala_ref_lh_path)
            clipped_img = nb.Nifti1Image(image_tmp, tha_mask_l_.affine, tha_mask_l_.header)
            
            nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_gradient{}_left_tha_aligned.nii.gz'.format(sub,g+1))

if alignment == False:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/all_individual_gradients_lh.npy", sub_grad)
elif alignment == True:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/all_individual_gradients_lh_aligned.npy", sub_grad_aligned)
########## correlation between individual gradients and grouplevel gradient

#gradient1 
correlation1=np.zeros(50)
for s in range(50):
    if alignment == False:
        r = np.corrcoef(group_grad[:,0],sub_grad[:,0,s])
    elif alignment == True:
        r = np.corrcoef(group_grad[:,0],sub_grad_aligned[:,0,s])   
    correlation1[s] = r[0,1] 

#gradient2
    correlation2=np.zeros(50)
    for s in range(50):
        if alignment == False:
            r = np.corrcoef(group_grad[:,1],sub_grad[:,1,s])
        elif alignment == True:
            r = np.corrcoef(group_grad[:,1],sub_grad_aligned[:,1,s])
        correlation2[s] = r[0,1] 

if alignment == False:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_1_lh.npy",correlation1)
    np.save("Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_2_lh.npy",correlation2)

elif alignment == True:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_1_lh_aligned.npy",correlation1)
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_2_lh_aligned.npy",correlation2)



""" right hemisphere """

# import refined thalamus mask as reference
thala_ref_rh_path="/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh_refined.nii.gz"
thala_ref_rh=nb.load(thala_ref_rh_path).get_fdata()
group_grad_r=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/gradients_rh.npy")

# import structural connectivity matrix stack (50 subjects)
conn_matrix_sub_r=np.load("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/struc_conn_matrix_rh_sub.npy")
sub_grad_r=np.zeros((1029,10,50))
sub_grad_aligned_r=np.zeros((1029,10,50))

for i in np.arange(0,50):
    
    sub=i+1
    # paths to save results
    gradients = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_gradients_rh.npy".format(sub)
    affinity = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_affinity_rh.npy".format(sub)
    lambdas = "/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_lambdas_rh.npy".format(sub)
    
    sparsity=0.75   # set threshhold for computing the gradient
    
    if alignment == False:
    
    # COMPUTE GRADIENT
        gm_r = GradientMaps(n_components=10, random_state=0, approach='dm', kernel='normalized_angle')
        gm_r.fit(conn_matrix_sub_r[:,:,i], sparsity=sparsity)
        np.save(gradients, gm_r.gradients_)
        np.save(lambdas, gm_r.lambdas_)
        sub_grad_r[:,:,i]=gm_r.gradients_

        # project gradient values onto thalamus 
        # iterate over first 3 gradients
        for g in range(2):
            image_tmp=np.zeros(thala_ref_rh.shape)          # create empty image same size as thalamus mask
            idx_r=zip(*np.where(thala_ref_rh==1))           # collect indices of voxels where mask = 1 
            #iterate over indices
            for i, idx in enumerate(idx_r):         
                image_tmp[idx]=gm_r.gradients_[i, g]        # fill in gradients value in thalamus mask     
        
            # save as nifti, with same header informations as thalamus mask
            tha_mask_r_=nb.load(thala_ref_rh_path)
            clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)
        
        nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_gradient{}_right_tha_.nii.gz'.format(sub,g+1))

    elif alignment == True:        
        # use procrustes alignment to align individual gradients to the grouplevel 

        align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                               kernel='normalized_angle', alignment='procrustes')  
        #align all individual gradients to grouplevel
        align.fit(conn_matrix_sub_r[:,:,i],reference=group_grad_r, sparsity=sparsity)
        grad_aligned = align.aligned_
        sub_grad_aligned_r[:,:,i]=grad_aligned
            
        # project gradient values onto thalamus 
        # iterate over first 2 gradients        
        for g in range(2):
            image_tmp=np.zeros(thala_ref_rh.shape)          # create empty image same size as thalamus mask
            idx_r=zip(*np.where(thala_ref_rh==1))           # collect indices of voxels where mask = 1 
            #iterate over indices
            for i, idx in enumerate(idx_r):         
                image_tmp[idx]=align.aligned_[i, g]        # fill in gradients value in thalamus mask     
                
            # save as nifti, with same header informations as thalamus mask
            tha_mask_r_=nb.load(thala_ref_rh_path)
            clipped_img = nb.Nifti1Image(image_tmp, tha_mask_r_.affine, tha_mask_r_.header)
                
            nb.save(clipped_img, '/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/sub-{}_gradient{}_right_tha_aligned.nii.gz'.format(sub,g+1))

if alignment == False:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/all_individual_gradients_rh.npy", sub_grad_r)
elif alignment == True:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/all_individual_gradients_rh_aligned.npy", sub_grad_aligned_r)

########## correlation between individual gradients and grouplevel gradient

#gradient1 
correlation1=np.zeros(50)
for s in range(50):
    if alignment == False:
        r = np.corrcoef(group_grad_r[:,0],sub_grad_r[:,0,s])
    elif alignment == True:
        r = np.corrcoef(group_grad_r[:,0],sub_grad_aligned_r[:,0,s])   
    correlation1[s] = r[0,1] 

#gradient2
    correlation2=np.zeros(50)
    for s in range(50):
        if alignment == False:
            r = np.corrcoef(group_grad_r[:,1],sub_grad_r[:,1,s])
        elif alignment == True:
            r = np.corrcoef(group_grad_r[:,1],sub_grad_aligned_r[:,1,s])
        correlation2[s] = r[0,1] 

if alignment == False:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_1_rh.npy",correlation1)
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_2_rh.npy",correlation2)

elif alignment == True:
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_1_rh_aligned.npy",correlation1)
    np.save("/Project1_thalamus_gradients/data/structural_connectivity/parcels_200/individual_gradients/corr_individual_group_gradients_2_rh_aligned.npy",correlation2)

