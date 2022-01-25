"""
computes and aligns group-level isocortex gradients
"""
import os
import h5py
import numpy as np
from brainspace.gradient import GradientMaps

# data dirs
ddir     = '../data/'   
conndir  = '../data/tout_hippoc/' 
odir     = '../data/tout_group/'

# get HCP - S900 subject list        
subjlist = '../data/subjectListS900_QC_gr.txt'
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
subjlist = joinedlist = mylist[:-1]
print('We have now %i subjects... ' % (len(subjlist)))  # 709 subjects

# get isocortex X isocortex connectivity

glassdir = '../data/glasserTimeseries/'

# it is HCP data, so there will be 4 scans
scans = ['rfMRI_REST1_LR', 'rfMRI_REST1_RL', 
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL']


Gcon_all = np.zeros((360,360))

## subjID = '100307'

i = 0
for subjID in subjlist:
    i += 1
    # get data files
    subj_glass_file = os.path.join(glassdir, 'HCP_' + subjID +
                                   '_glasserTimeseries.mat')

    #  HDF reader for matlab v7.3 files
    f_subj_glass = h5py.File(subj_glass_file, 'r')         

    for scan in scans:
        subj_glass  = np.array(f_subj_glass[scan])             # (360, 1200)
        Gcon = np.corrcoef(subj_glass)                         # (360, 360) 
        Gcon_all += Gcon

Gcon_all = Gcon_all / (len(subjlist) * len(scans))

print('We averaged isocortex-to-isocortex conn. from %i subjects...' % (i))
print('HHHH ', Gcon_all.shape, np.ndim(Gcon_all), Gcon_all.min(), Gcon_all.max())

h = h5py.File(os.path.join(odir, 'Hmean709isocortex_fcon.h5'), 'w')
h.create_dataset('data', data = Gcon_all)
h.close()


# get gradients
gm = GradientMaps(approach = 'dm', kernel='cosine');
gm = gm.fit(Gcon_all);   # (360, 360)


h = h5py.File(os.path.join(odir, 'Hmean709isocortex_gradients.h5'), 'w')
h.create_dataset('gradients_', data = gm.gradients_ )
h.create_dataset('lambdas_', data = gm.lambdas_ )
h.close()

