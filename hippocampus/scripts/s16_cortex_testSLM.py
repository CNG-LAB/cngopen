"""
SLM test for the cortex-to-hippocampus connectivity for individual subfields
usage: $ python s16_cortex_testSLM.py LSUB
"""
import os, sys
import h5py
import numpy as np
from numpy import genfromtxt 

# definde data directories
ddir      = '../data/'                                          #  data dir
cordir    = '../data/tout_cortex/'
odir      = '../data/tout_group'

# final subject list after QC         
subjlist = os.path.join(ddir, 'subjectListS900_QC_gr.txt');     # 709 subjects
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
mylist =  mylist[:-1] 
totnum = len(mylist)

labeling_file = '../data/tout_group/glasser.csv' 
mylabel = genfromtxt(labeling_file)
print('We have now %i subjects... ' % totnum)

# subfield = 'LSUB'
subfield = sys.argv[1]

# here we go
C360_all = np.zeros((len(mylist), 360))
i = 0
for subjID in mylist:
    subjsub= os.path.join(cordir, subjID + '_cortex_%s.h5' % (subfield))
    with h5py.File(subjsub, "r") as f:        
        subjdata = np.array(f[subjID])    
    C360_all[i, :] = subjdata.T
    i +=1   

print(C360_all.shape, C360_all.mean(axis=0).max())

# labeling from 360 to 64k points
C64k_all = np.zeros((len(mylist), 64984))
for i in range(0, len(mylist)):
    for j in range(1,360+1):
        C64k_all[i, np.where(mylabel == j)] = C360_all[i,(j-1)]

print(C64k_all.shape, C64k_all.mean(axis=0).max())

from brainspace.datasets import load_conte69
from brainspace.mesh import mesh_elements

# load poly data for 64k surface (for the test & plotting)
surf_lh, surf_rh = load_conte69()

# write surface coordinates and triangles in a dictionary
lh_coord = np.array(mesh_elements.get_points(surf_lh)).T
rh_coord = np.array(mesh_elements.get_points(surf_rh)).T
lh_tri = np.array(mesh_elements.get_cells(surf_lh))
rh_tri = np.array(mesh_elements.get_cells(surf_rh))

D = {}
D['coord'] = np.concatenate((lh_coord, rh_coord), axis=1)         # (3, 64984)     
D['tri'] = np.concatenate((lh_tri, rh_tri + lh_coord.shape[1]))   # (129960, 3)


# run slm
from brainstat.stats.terms import FixedEffect
from brainstat.stats.SLM import SLM

Y         = C64k_all
contrast  = np.ones((len(mylist),1))
term_     = FixedEffect(contrast)
model_    = 1 + term_

slm = SLM(model_, contrast = contrast)
slm.fit(Y)
Tvals = slm.t
Tvals.shape

h = h5py.File(os.path.join(odir, 'Tvals_cortex709_%s.h5' % (subfield)), 'w')
h.create_dataset('data', data = Tvals)
h.close()

