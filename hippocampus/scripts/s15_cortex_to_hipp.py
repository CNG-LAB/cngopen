"""
calculates the connectivity for cortex-to-hippocampus fashion
usage: $ python s15_cortex_to_hipp.py 100610
"""
import os, sys
import h5py
import numpy as np

# definde data directories
ddir     = '../data/';            
odir     = '../data/tout_cortex'
glassdir = os.path.join(ddir, 'glasserTimeseries/');    # cortex t-series
hippdir  = os.path.join(ddir, 'smoothTimeseries/');     # hippocampus t-series

# it is HCP data, so there will be 4 scans
scans = ['rfMRI_REST1_LR', 'rfMRI_REST1_RL', 
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL']

# hippocampus segmentations
roi_lsub = 'L_SUB'
roi_lca  = 'L_CA'
roi_ldg  = 'L_DG'
roi_rsub = 'R_SUB'
roi_rca  = 'R_CA'
roi_rdg  = 'R_DG'

glasservertexnum = 360;
sum_LSUB = np.zeros((glasservertexnum, 1))
sum_LCA  = np.zeros((glasservertexnum, 1))
sum_LDG  = np.zeros((glasservertexnum, 1))
sum_RSUB = np.zeros((glasservertexnum, 1))
sum_RCA  = np.zeros((glasservertexnum, 1))
sum_RDG  = np.zeros((glasservertexnum, 1))

#subjID = '100610'
subjID = sys.argv[1]

# get data files
subj_glass_file = os.path.join(glassdir, 'HCP_' + subjID +
                               '_glasserTimeseries.mat')
subj_hipp_file = os.path.join(hippdir, 'HCP_' + subjID +
                              '_smoothTimeseries.mat');

#  HDF reader for matlab v7.3 files
f_subj_glass = h5py.File(subj_glass_file, 'r')         
f_subj_hipp  = h5py.File(subj_hipp_file, 'r')           

# LSUB
for scan in scans:
    roi = roi_lsub
    # arrays
    subj_glass = np.array(f_subj_glass[scan]).T              # (1200, 360)
    subj_hipp = np.array(f_subj_hipp[scan][roi]).T           # (1200, 1024)
    subj_hippav = subj_hipp.mean(axis=1)                     # (1200,)
    subj_hippav = subj_hippav.reshape(len(subj_hippav), 1)   # (1200, 1)
    # connectivity from cortex to hippocampus segmentation
    subj_corr = np.corrcoef(subj_glass, subj_hippav, 
                            rowvar=False)[-1:,:-1].T         # (360, 1)
    # summing it up across roi & scans (Fisher r2z)
    sum_LSUB = sum_LSUB + np.arctanh(subj_corr)

# get average cortex connectivity
sum_LSUB = sum_LSUB / len(scans); 

h = h5py.File(os.path.join(odir, subjID + '_cortex_LSUB.h5'), 'w')
h.create_dataset(subjID, data = sum_LSUB)
h.close()

# LCA
for scan in scans:
    # arrays
    subj_glass = np.array(f_subj_glass[scan]).T              # (1200, 360)
    subj_hipp = np.array(f_subj_hipp[scan][roi_lca]).T       # (1200, 2048)
    subj_hippav = subj_hipp.mean(axis=1)                     # (1200,)
    subj_hippav = subj_hippav.reshape(len(subj_hippav), 1)   # (1200, 1)
    # connectivity from cortex to hippocampus segmentation
    subj_corr = np.corrcoef(subj_glass, subj_hippav, 
                            rowvar=False)[-1:,:-1].T         # (360, 1)
    # summing it up across roi & scans (Fisher r2z)
    sum_LCA = sum_LCA + np.arctanh(subj_corr)

sum_LCA  = sum_LCA / len(scans) 

h = h5py.File(os.path.join(odir, subjID + '_cortex_LCA.h5'), 'w')
h.create_dataset(subjID, data = sum_LCA)
h.close()

# LDG
for scan in scans:
    # arrays
    subj_glass = np.array(f_subj_glass[scan]).T              # (1200, 360)
    subj_hipp = np.array(f_subj_hipp[scan][roi_ldg]).T       # (1200, 1024)
    subj_hippav = subj_hipp.mean(axis=1)                     # (1200, )
    subj_hippav = subj_hippav.reshape(len(subj_hippav), 1)   # (1200, 1)
    # connectivity from cortex to hippocampus segmentation
    subj_corr = np.corrcoef(subj_glass, subj_hippav, 
                            rowvar=False)[-1:,:-1].T         # (360, 1)
    # summing it up across roi & scans (Fisher r2z)
    sum_LDG = sum_LDG + np.arctanh(subj_corr)
    
sum_LDG  = sum_LDG / len(scans) ; 

h = h5py.File(os.path.join(odir, subjID + '_cortex_LDG.h5'), 'w')
h.create_dataset(subjID, data = sum_LDG)
h.close()

# RSUB
for scan in scans:
    # arrays
    subj_glass = np.array(f_subj_glass[scan]).T              # (1200, 360)
    subj_hipp = np.array(f_subj_hipp[scan][roi_rsub]).T      # (1200, 1024)
    subj_hippav = subj_hipp.mean(axis=1)                     # (1200,)
    subj_hippav = subj_hippav.reshape(len(subj_hippav), 1)   # (1200, 1)
    # connectivity from cortex to hippocampus segmentation
    subj_corr = np.corrcoef(subj_glass, subj_hippav, 
                            rowvar=False)[-1:,:-1].T         # (360, 1)
    # summing it up across roi & scans (Fisher r2z)
    sum_RSUB = sum_RSUB + np.arctanh(subj_corr)

# get average cortex connectivity
sum_RSUB = sum_RSUB / len(scans); 

h = h5py.File(os.path.join(odir, subjID + '_cortex_RSUB.h5'), 'w')
h.create_dataset(subjID, data = sum_RSUB)
h.close()

# RCA
for scan in scans:
    # arrays
    subj_glass = np.array(f_subj_glass[scan]).T              # (1200, 360)
    subj_hipp = np.array(f_subj_hipp[scan][roi_rca]).T       # (1200, 2048)
    subj_hippav = subj_hipp.mean(axis=1)                     # (1200,)
    subj_hippav = subj_hippav.reshape(len(subj_hippav), 1)   # (1200, 1)
    # connectivity from cortex to hippocampus segmentation
    subj_corr = np.corrcoef(subj_glass, subj_hippav, 
                            rowvar=False)[-1:,:-1].T         # (360, 1)
    # summing it up across roi & scans (Fisher r2z)
    sum_RCA = sum_RCA + np.arctanh(subj_corr)

sum_RCA  = sum_RCA / len(scans) 

h = h5py.File(os.path.join(odir, subjID + '_cortex_RCA.h5'), 'w')
h.create_dataset(subjID, data = sum_RCA)
h.close()

# RDG
for scan in scans:
    # arrays
    subj_glass = np.array(f_subj_glass[scan]).T              # (1200, 360)
    subj_hipp = np.array(f_subj_hipp[scan][roi_rdg]).T       # (1200, 1024)
    subj_hippav = subj_hipp.mean(axis=1)                     # (1200, )
    subj_hippav = subj_hippav.reshape(len(subj_hippav), 1)   # (1200, 1)
    # connectivity from cortex to hippocampus segmentation
    subj_corr = np.corrcoef(subj_glass, subj_hippav, 
                            rowvar=False)[-1:,:-1].T         # (360, 1)
    # summing it up across roi & scans (Fisher r2z)
    sum_RDG = sum_RDG + np.arctanh(subj_corr)
    
sum_RDG  = sum_RDG / len(scans) ; 

h = h5py.File(os.path.join(odir, subjID + '_cortex_RDG.h5'), 'w')
h.create_dataset(subjID, data = sum_RDG)
h.close()
