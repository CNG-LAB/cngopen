"""
computes the individual-level correlations between T1w/T2w and G2 (fc) maps
"""
import os
import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# get the list of subject id's
subjlist = '../data/subjectListS900_QC_gr.txt'
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
subjlist = joinedlist = mylist[:-1]
len(subjlist)

# create a dataframe with subfield columns
mycols = ['tG2_LSUB', 'tG2_LCA', 'tG2_LDG']
dcorr = pd.DataFrame(index = subjlist, columns = mycols)
dcorr.index.name = 'id'

# for each subject, compute the correlation between T1w/T2w and G2 (fc)
for subjID in subjlist:

    tdir = '../data/tout_T1wT2w_msm50/'

    vol2hipp_LSUB  = os.path.join(tdir, 'HCP_%s_t1t2_sub_left.h5' % (subjID))
    h_LSUB  = h5py.File(vol2hipp_LSUB, 'r')
    t_LSUB  = np.array(h_LSUB[subjID])

    vol2hipp_LCA  = os.path.join(tdir, 'HCP_%s_t1t2_ca_left.h5' % (subjID))
    h_LCA   = h5py.File(vol2hipp_LCA, 'r')
    t_LCA   = np.array(h_LCA[subjID])

    vol2hipp_LDG  = os.path.join(tdir, 'HCP_%s_t1t2_dg_left.h5' % (subjID))
    h_LDG   = h5py.File(vol2hipp_LDG, 'r')
    t_LDG   = np.array(h_LDG[subjID])

    gdir = '../data/tout_hippoc_grad_flipped_msm50/'
    
    gfile_LSUB = h5py.File(os.path.join(gdir,'HCP_' +subjID+'_G2_LSUB.h5'),'r')
    g2_LSUB = np.array(gfile_LSUB[subjID])  
    gfile_LSUB.close()

    gfile_LCA = h5py.File(os.path.join(gdir, 'HCP_' +subjID+ '_G2_LCA.h5'),'r')
    g2_LCA = np.array(gfile_LCA[subjID])  
    gfile_LCA.close()

    gfile_LDG = h5py.File(os.path.join(gdir, 'HCP_' +subjID+ '_G2_LDG.h5'),'r')
    g2_LDG = np.array(gfile_LDG[subjID])  
    gfile_LDG.close()
  
    iC = dcorr.index.get_loc(subjID)

    dcorr.iloc[iC]['tG2_LSUB'] = pearsonr(np.log(t_LSUB), g2_LSUB)[0]
    dcorr.iloc[iC]['tG2_LCA']  = pearsonr(np.log(t_LCA), g2_LCA)[0]
    dcorr.iloc[iC]['tG2_LDG']  = pearsonr(np.log(t_LDG), g2_LDG)[0]
    
# save the dataframe   
dcorr.to_csv('../data/tout_group/Hmean709_t1wt2_corr_G2.csv')

