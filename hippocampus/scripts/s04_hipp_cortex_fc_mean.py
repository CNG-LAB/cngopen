"""
computes the mean hippocampal-cortical functional connectivity (fc) matrix,
for the left hemisphere subfields
"""
import os
import h5py
import numpy as np

# data dirs
ddir     = '../data/'   
conndir  = '../data/tout_hippoc/' 
odir     = '../data/tout_group/'

# get HCP - S900 subject list        
subjlist = '../data/subjectListS900_QC_gr.txt'
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
subjlist = joinedlist = mylist[:-1]
print('We have now %i subjects... ' % (len(subjlist)))  # 709


fc_left = np.zeros((4096, 360))

j = 0
for subjID in subjlist:

    fname = os.path.join(conndir, 'HCP_' + subjID + '_left.h5')
    f = h5py.File(fname, 'r')
    f = np.array(f['HCP_' + subjID])

    fc_left = fc_left + f
    j += 1
    
fc_left = fc_left / j

h = h5py.File('../data/tout_group/Hmean709_FC_left.h5', 'w')
h.create_dataset('data', data = fc_left)
h.close()

print(fc_left.min(), fc_left.max(), fc_left.shape, j)
# -0.005300521852874321, 0.39153784016161197, (4096, 360), 709