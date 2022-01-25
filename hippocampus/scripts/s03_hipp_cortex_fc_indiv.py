"""
computes the individual hippocampus-to-cortex FC matrix (4096, 360),
separately for the left and right hemispheres
usage: $ python s03_hipp_cortex_fc_indiv.py HCP_100610
"""
import os, sys
import numpy as np
import h5py
from scipy.stats import pearsonr

# data dir
ddir     = '../data/'  
odir     = '../data/tout_hippoc'
                
glassdir = os.path.join(ddir, 'glasserTimeseries/');   # cortex t-series
hippdir  = os.path.join(ddir, 'smoothTimeseries/');    # hippocampus t-series

# it is HCP data, so there will be 4 scans
scans = ['rfMRI_REST1_LR', 'rfMRI_REST1_RL', 
         'rfMRI_REST2_LR', 'rfMRI_REST2_RL'];

# hippocampus segmentations
rois = ['L_SUB', 'L_CA', 'L_DG', 'R_SUB', 'R_CA', 'R_DG']

#subjID = 'HCP_100610'
subjID = sys.argv[1]

# get file full paths for each subject
subj_glass_file = os.path.join(glassdir, 
                               subjID + '_glasserTimeseries.mat')

subj_hipp_file = os.path.join(hippdir,
                              subjID + '_smoothTimeseries.mat')

#  HDF reader for matlab v7.3 files
f_subj_glass = h5py.File(subj_glass_file, 'r')
f_subj_hipp  = h5py.File(subj_hipp_file, 'r')

# zeros array for hippocampus-to-cortex FC
corr_hleft_all = np.zeros((4096, 360))
corr_hright_all = np.zeros((4096, 360))

for scan in scans:
    for roi in rois:

        subj_glass = np.array(f_subj_glass[scan]).T    # [1200x360]
        subj_hipp = np.array(f_subj_hipp[scan][roi]).T  

        # subj_hipp:
        #     L_SUB: [1200×1024]
        #     L_CA:  [1200×2048]
        #     L_DG:  [1200×1024]
        #     R_SUB: [1200×1024]
        #     R_CA:  [1200×2048]
        #     R_DG:  [1200×1024]

        # focus on left hippocampus
        if roi == 'L_SUB' or roi == 'L_CA' or roi == 'L_DG':
            # concatenate L_SUB, L_CA, L_DG
            if roi == 'L_SUB':
                hleft = subj_hipp
            elif roi == 'L_CA':
                hleft = np.concatenate((hleft, subj_hipp), axis=1)
            elif roi == 'L_DG':
                hleft = np.concatenate((hleft, subj_hipp), axis=1)

        # focus on right hippocampus
        if roi == 'R_SUB' or roi == 'R_CA' or roi == 'R_DG':
            # concatenate R_SUB, R_CA, R_DG
            if roi == 'R_SUB':
                hright = subj_hipp
            elif roi == 'R_CA':
                hright = np.concatenate((hright, subj_hipp), axis=1)
            elif roi == 'R_DG':
                hright = np.concatenate((hright, subj_hipp), axis=1)
            
    # hleft is the concetantion of L_SUB, L_CA, L_DG --> [1200x4096]    
    # column-wise pearsonr between hleft and subj_glass --> [4096x360]
    corr_hleft = np.zeros((hleft.shape[1], subj_glass.shape[1])) 
    for hcol in range(hleft.shape[1]): 
        for gcol in range(subj_glass.shape[1]):
            corr_hleft[hcol, gcol] = pearsonr(hleft[:,hcol], 
                                              subj_glass[:,gcol])[0]

    # hright is the concetantion of R_SUB, R_CA, R_DG --> [1200x4096]    
    # column-wise pearsonr between rright and subj_glass --> [4096x360]
    corr_hright = np.zeros((hright.shape[1], subj_glass.shape[1])) 
    for hcol in range(hright.shape[1]): 
        for gcol in range(subj_glass.shape[1]):
            corr_hright[hcol, gcol] = pearsonr(hright[:,hcol], 
                                               subj_glass[:,gcol])[0]
      
    # add hippocampus-to-cortex connectivity along scans (fisher RtoZ)
    corr_hleft_all += np.arctanh(corr_hleft)
    corr_hright_all += np.arctanh(corr_hright)

# average hippocampus-to-cortex connectivity across scans    
corr_hleft_all = corr_hleft_all / len(scans)  # [4096x360]  
corr_hright_all = corr_hright_all / len(scans)  # [4096x360]  

print(subjID, corr_hleft_all.shape, corr_hright_all.shape)

h = h5py.File(os.path.join(odir, subjID + '_left.h5'), 'w')
h.create_dataset(subjID, data = corr_hleft_all)
h.close()

h = h5py.File(os.path.join(odir, subjID + '_right.h5'), 'w')
h.create_dataset(subjID, data = corr_hright_all)
h.close()
