"""
reads in the individual-level T1w/T2w ratios from blades dir,
saves it all together for all subjects,
also saves the mean T1w/T2w ratio for each subject in a data frame
"""
import os
import h5py
import numpy as np
import pandas as pd
import nibabel as nb

subjlist = '../data/subjectListS900_QC_gr.txt'
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
subjlist = joinedlist = mylist[:-1]
len(subjlist)

# create empty dataframe for the mean T1w/T2w intensity values 
# for each subject and each subfield
DF_indiv_mean_t1wt2w  = pd.DataFrame(index = subjlist, 
                                     columns = ['t1-t2-LSUB', 
                                                't1-t2-LCA', 
                                                't1-t2-LDG'])

# create empty arrays for the individual-level T1w/T2w intensity values 
all_LSUB = np.zeros((709,1024))
all_LCA  = np.zeros((709,2048))
all_LDG  = np.zeros((709,1024))

# compute the mean T1w/T2w ratios per subject and per subfield
# and fill out the empty arrays
workdir   = '../data/bladesMniGifti_T1wT2w/'

j = 0
for subjid in subjlist:

    # filenames for t1w-t2w ratios (resampled along hippocampus)
    vol2hipp_LSUB = os.path.join(workdir, 
                                 'HCP_%s_L_SUB_skelFinal.shape.gii' % (subjid))  
    vol2hipp_LCA  = os.path.join(workdir, 
                                 'HCP_%s_L_CA_skelFinal.shape.gii' % (subjid)) 
    vol2hipp_LDG  = os.path.join(workdir, 
                                 'HCP_%s_L_DG_skelFinal.shape.gii' % (subjid))  

    t1wt2w_LSUB = nb.load(vol2hipp_LSUB).agg_data()
    t1wt2w_LCA  = nb.load(vol2hipp_LCA).agg_data()
    t1wt2w_LDG  = nb.load(vol2hipp_LDG).agg_data()    
    
    DF_indiv_mean_t1wt2w.at[subjid, 't1-t2-LSUB'] = t1wt2w_LSUB.mean()     
    DF_indiv_mean_t1wt2w.at[subjid, 't1-t2-LCA' ] = t1wt2w_LCA.mean()     
    DF_indiv_mean_t1wt2w.at[subjid, 't1-t2-LDG' ] = t1wt2w_LDG.mean()     

    all_LSUB[j,:] = t1wt2w_LSUB
    all_LCA[j,:]  = t1wt2w_LCA
    all_LDG[j,:]  = t1wt2w_LDG
    
    j += 1    
        
print(j) # 709
        
# save mean T1w/T2w ratios per subject    
DF_indiv_mean_t1wt2w.to_csv('../data/tout_group/Hmean709_t1wt2_lsub_lca_ldg.csv')

# save T1w/T2w ratios per subject
h1 = h5py.File('../data/tout_group/H709_t1wt2w_lsub.h5', 'w')
h1.create_dataset('data', data = all_LSUB); h1.close()
h2 = h5py.File('../data/tout_group/H709_t1wt2w_lca.h5', 'w')
h2.create_dataset('data', data = all_LCA); h2.close()
h3 = h5py.File('../data/tout_group/H709_t1wt2w_ldg.h5', 'w')
h3.create_dataset('data', data = all_LDG); h3.close()
    