"""
gets the T1w/T2w ratios at cortical and subfield level,
for all subjects and saves it as matrix
"""
import os
import numpy as np
import nibabel as nb
import h5py
from numpy import genfromtxt
from brainspace.utils.parcellation import reduce_by_labels

# get HCP - S900 subject list        
subjlist = '../data/subjectListS900_QC_gr.txt'
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
subjlist = joinedlist = mylist[:-1]

mysubjects = []
for subj in subjlist:
    mysubjects.append(subj)
print(len(mysubjects))

# Glasser Labels
labeling_file = '../data/tout_group/glasser.csv' 
labeling = genfromtxt(labeling_file)
print(labeling.shape, labeling.min(), labeling.max())
# (64984,), 0.0, 360.0

# data dirs
ddir = '../data/hcp_smoothedmyelin/'
tdir = '../data/tout_T1wT2w_msm50/'
odir = '../data/tout_group/'

# get myelin measures for cortex across all subjects (360, 709)
count_f = 0
for subj in mysubjects:
    dataup = os.path.join(ddir, '%s/MNINonLinear/' %(subj))
    
    fLeft  = os.path.join(dataup, 'fsaverage_LR32k/%s.L.SmoothedMyelinMap.32k_fs_LR.func.gii' %(subj))
    fRight = os.path.join(dataup, 'fsaverage_LR32k/%s.R.SmoothedMyelinMap.32k_fs_LR.func.gii'  %(subj))

    if os.path.isfile(fLeft) and os.path.isfile(fRight):
 
        dLeft  = nb.load(fLeft).agg_data()   # 32k
        dRight = nb.load(fRight).agg_data()  # 32k
        # concatenate Left & Right hemis
        dLR = np.concatenate((dLeft, dRight), axis=0) # 64k
        # mapping it to Glasser parcellation    
        dLR_glass = reduce_by_labels(dLR, labels=labeling)[1:]
        
        if count_f == 0:
            mye_cortex = dLR_glass
            mye_cortex = mye_cortex.reshape(-1, 1)
        else:
            mye_cortex = np.concatenate((mye_cortex, dLR_glass.reshape(-1,1)), axis=1)
        count_f += 1         

print(count_f, mye_cortex.shape, mye_cortex.min(), mye_cortex.max())
# 709, (360, 709), 0.02919340319931507, 4.150848388671875)
h1 = h5py.File(os.path.join(odir, 'H709_mye_cortex.h5'), 'w')
h1.create_dataset('data', data = mye_cortex); h1.close()


# get myelin measures for subfields across all subjects (4096, 709)
count_t = 0
for subj in mysubjects:
    dataup = os.path.join(ddir, '%s/MNINonLinear/' %(subj))
    fLeft  = os.path.join(dataup, 'fsaverage_LR32k/%s.L.SmoothedMyelinMap.32k_fs_LR.func.gii' %(subj))
    fRight = os.path.join(dataup, 'fsaverage_LR32k/%s.R.SmoothedMyelinMap.32k_fs_LR.func.gii' %(subj))
   
    if os.path.isfile(fLeft) and os.path.isfile(fRight):
        
        vol2hipp_LSUB = os.path.join(tdir, 'HCP_%s_t1t2_sub_left.h5' % (subj))
        vol2hipp_LCA = os.path.join(tdir, 'HCP_%s_t1t2_ca_left.h5' % (subj))
        vol2hipp_LDG = os.path.join(tdir, 'HCP_%s_t1t2_dg_left.h5' % (subj))
        
        h_LSUB = h5py.File(vol2hipp_LSUB, 'r')
        h_LCA = h5py.File(vol2hipp_LCA, 'r')
        h_LDG = h5py.File(vol2hipp_LDG, 'r')
        
        dataLSUB = np.array(h_LSUB[subj])
        dataLCA = np.array(h_LCA[subj])
        dataLDG = np.array(h_LDG[subj])

        if count_t == 0:
            mye_LSUB = dataLSUB.reshape(-1,1)
            mye_LCA = dataLCA.reshape(-1,1)
            mye_LDG = dataLDG.reshape(-1,1)
        
        else:
            mye_LSUB = np.concatenate((mye_LSUB, dataLSUB.reshape(-1,1)), axis=1)
            mye_LCA = np.concatenate((mye_LCA, dataLCA.reshape(-1,1)), axis=1)
            mye_LDG = np.concatenate((mye_LDG, dataLDG.reshape(-1,1)), axis=1)

        count_t += 1

# concatenate myelin across subfields
mye_subfields = np.concatenate((mye_LSUB, mye_LCA, mye_LDG))
print(mye_subfields.shape, mye_subfields.min(), mye_subfields.max())
# (4096, 709), 0.8110781, 33.222504
h2 = h5py.File(os.path.join(odir, 'H709_mye_subfields.h5'), 'w')
h2.create_dataset('data', data = mye_subfields); h2.close()
