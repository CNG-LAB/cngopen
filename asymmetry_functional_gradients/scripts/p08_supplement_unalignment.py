"""
computes the unaligned gradients of functional connectivity for 4 different fashions:
LL & RR intra-hemispheric and LR & RL inter-hemispheric fashion.

change the input of HCP or UKB
"""
import os
import pandas as pd
import numpy as np
from brainspace.gradient import GradientMaps

datadir = '../data/data_results/supplementary/ukb/FC/'

graddir = '../data/data_results/supplementary/ukb/gradient/unalign/'

path = os.path.join(datadir, 'LL')
path_list = os.listdir(path)
path_list.sort()

n = len(path_list)
group_grad = np.array(pd.read_csv('../data/data_results/supplementary/ukb/gradient/group_grad_LL.csv', header=None))

# get the gradients for each subject and each FC (LL, RR, LR, RL) and
# align them to the group-level LL gradients, aligned_ = aligned gradients, gradients_ = non-aligned gradients

for i in path_list:
  # FC LL
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment = 'procrustes')  
  fc_LL = np.array(pd.read_csv(os.path.join(datadir, 'LL', i), header=None))
  align.fit(fc_LL, reference=group_grad)
  grad_LL = align.gradients_
  for j in range(10):
    if np.corrcoef(grad_LL[:,j], group_grad[:,j])[0][1] < 0:
      grad_LL[:,j] = -grad_LL[:,j]
  np.savetxt(os.path.join(graddir, 'LL', i), grad_LL, delimiter = ',')

  # FC RR
  align = GradientMaps(n_components=10, random_state=0, approach='dm',
                       kernel='normalized_angle', alignment = 'procrustes')
  fc_RR = np.array(pd.read_csv(os.path.join(datadir, 'RR', i), header=None))
  align.fit(fc_RR, reference=group_grad)
  grad_RR = align.gradients_
  for j in range(10):
    if np.corrcoef(grad_RR[:,j], group_grad[:,j])[0][1] < 0:
      grad_RR[:,j] = -grad_RR[:,j]
  np.savetxt(os.path.join(graddir, 'RR', i), grad_RR, delimiter = ',')  

  # FC LR 
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment = 'procrustes')
  fc_LR = np.array(pd.read_csv(os.path.join(datadir, 'LR', i), header=None))
  align.fit(fc_LR, reference=group_grad)
  grad_LR = align.gradients_
  for j in range(10):
    if np.corrcoef(grad_LR[:,j], group_grad[:,j])[0][1] < 0:
      grad_LR[:,j] = -grad_LR[:,j]
  np.savetxt(os.path.join(graddir, 'LR', i), grad_LR, delimiter = ',')
  
  # FC RL  
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment = 'procrustes')
  
  fc_RL = np.array(pd.read_csv(os.path.join(datadir, 'RL', i), header=None))
  align.fit(fc_RL, reference=group_grad)
  grad_RL = align.gradients_
  for j in range(10):
    if np.corrcoef(grad_RL[:,j], group_grad[:,j])[0][1] < 0:
      grad_RL[:,j] = -grad_RL[:,j]
  np.savetxt(os.path.join(graddir, 'RL', i), grad_RL, delimiter = ',')

  # intra asymmetry
  grad_intra = grad_LL - grad_RR
  np.savetxt(os.path.join(graddir, 'intra', i), grad_intra, delimiter = ',')

  # inter asymmetry
  grad_inter = grad_LR - grad_RL
  np.savetxt(os.path.join(graddir, 'inter', i), grad_inter, delimiter = ',')
  print('finish   ' + i)
