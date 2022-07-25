"""
computes the gradients of functional connectivity for:
LL & RR intra-hemispheric and LR & RL inter-hemispheric fashion. Here we combined 
LL and RR together, LR and RL together, to investigate potential normalization differences
resulting from computing separate gradients for left and right hemisphere connectomes.

change the input of HCP or UKB
"""
import os
import pandas as pd
import numpy as np
from brainspace.gradient import GradientMaps
import statistics
from scipy import stats

datadir = '../data/data_results/FC/'

graddir = '../data/data_results/gradient/same_model/'

# LLRR
fc_1 = np.array(pd.read_csv('../data/data_results/FC/LL_groupmean.csv', header=None))
fc_2 = np.array(pd.read_csv('../data/data_results/FC/RR_groupmean.csv', header=None))
gm = GradientMaps(n_components=10, random_state=0, approach='dm', 
                  kernel='normalized_angle')
gm.fit(np.vstack((fc_1, fc_2)))
np.savetxt(graddir+'group_grad_LLRR.csv', gm.gradients_, delimiter = ',')
np.savetxt(graddir+'group_grad_LLRR_lambdas.csv', gm.lambdas_, delimiter = ',')

# LRRL
fc_1 = np.array(pd.read_csv('../data/data_results/FC/LR_groupmean.csv', header=None))
fc_2 = np.array(pd.read_csv('../data/data_results/FC/RL_groupmean.csv', header=None))
gm = GradientMaps(n_components=10, random_state=0, approach='dm', 
                  kernel='normalized_angle')
gm.fit(np.vstack((fc_1, fc_2)))
np.savetxt(graddir+'group_grad_LRRL.csv', gm.gradients_, delimiter = ',')
np.savetxt(graddir+'group_grad_LRRL_lambdas.csv', gm.lambdas_, delimiter = ',')

# individual level
path = os.path.join(datadir, 'LL')
path_list = os.listdir(path)
path_list.sort()

n = len(path_list)

group_grad_LLRR = np.array(pd.read_csv(graddir+'group_grad_LLRR.csv', header=None))
group_grad_LRRL = np.array(pd.read_csv(graddir+'group_grad_LRRL.csv', header=None))

# get the gradients for each subject and each FC (LL, RR, LR, RL) and
# align them to the group-level LL gradients
for i in path_list:
  # FC LLRR
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LL = np.array(pd.read_csv(os.path.join(datadir, 'LL', i), header=None))
  fc_RR = np.array(pd.read_csv(os.path.join(datadir, 'RR', i), header=None))
  align.fit(np.vstack((fc_LL, fc_RR)), reference=group_grad_LLRR)
  grad_LLRR = align.aligned_
  for j in range(10):
    if np.corrcoef(grad_LLRR[:,j], group_grad_LLRR[:,j])[0][1] < 0:
      grad_LLRR[:,j] = -grad_LLRR[:,j]
  np.savetxt(os.path.join(graddir, 'LLRR', i), grad_LLRR, delimiter = ',')
  np.savetxt(os.path.join(graddir, 'intra', i), grad_LLRR[:180]-grad_LLRR[180:], delimiter = ',')

  # FC LRRL 
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LR = np.array(pd.read_csv(os.path.join(datadir, 'LR', i), header=None))
  fc_RL = np.array(pd.read_csv(os.path.join(datadir, 'RL', i), header=None))
  align.fit(np.vstack((fc_LR, fc_RL)), reference=group_grad_LRRL)
  grad_LRRL = align.aligned_
  for j in range(10):
    if np.corrcoef(grad_LRRL[:,j], group_grad_LRRL[:,j])[0][1] < 0:
      grad_LRRL[:,j] = -grad_LRRL[:,j]
  np.savetxt(os.path.join(graddir, 'LRRL', i), grad_LRRL, delimiter = ',')
  np.savetxt(os.path.join(graddir, 'inter', i), grad_LRRL[:180]-grad_LRRL[180:], delimiter = ',')
  print('finish   ' + i)

cadir= '../data'

ca_l = np.array(pd.read_csv(os.path.join(cadir, 'ca_glasser_network.csv'),
                                         header=None))[:,0][:180]
ca_r = np.array(pd.read_csv(os.path.join(cadir, 'ca_glasser_network.csv'),
                                         header=None))[:,0][180:]

for n in range(len(path_list)):
  ll = np.array(pd.read_csv(os.path.join(graddir,'LLRR',path_list[n]),header=None))[:180]
  rr = np.array(pd.read_csv(os.path.join(graddir,'LLRR',path_list[n]),header=None))[180:]
  lr = np.array(pd.read_csv(os.path.join(graddir,'LRRL',path_list[n]),header=None))[:180]
  rl = np.array(pd.read_csv(os.path.join(graddir,'LRRL',path_list[n]),header=None))[180:]
  
  # intra-hemisphere
  intra = [None] * 3  
  for i in range(3):
    intra[i] = (statistics.mean(ll[:,i][np.where(ca_l==1)])-statistics.mean(rr[:,i][np.where(ca_r==1)]),
                statistics.mean(ll[:,i][np.where(ca_l==2)])-statistics.mean(rr[:,i][np.where(ca_r==2)]),
                statistics.mean(ll[:,i][np.where(ca_l==3)])-statistics.mean(rr[:,i][np.where(ca_r==3)]),
                statistics.mean(ll[:,i][np.where(ca_l==4)])-statistics.mean(rr[:,i][np.where(ca_r==4)]),
                statistics.mean(ll[:,i][np.where(ca_l==5)])-statistics.mean(rr[:,i][np.where(ca_r==5)]),
                statistics.mean(ll[:,i][np.where(ca_l==6)])-statistics.mean(rr[:,i][np.where(ca_r==6)]),
                statistics.mean(ll[:,i][np.where(ca_l==7)])-statistics.mean(rr[:,i][np.where(ca_r==7)]),
                statistics.mean(ll[:,i][np.where(ca_l==8)])-statistics.mean(rr[:,i][np.where(ca_r==8)]),
                statistics.mean(ll[:,i][np.where(ca_l==9)])-statistics.mean(rr[:,i][np.where(ca_r==9)]),
                statistics.mean(ll[:,i][np.where(ca_l==10)])-statistics.mean(rr[:,i][np.where(ca_r==10)]),
                statistics.mean(ll[:,i][np.where(ca_l==11)])-statistics.mean(rr[:,i][np.where(ca_r==11)]),
                statistics.mean(ll[:,i][np.where(ca_l==12)])-statistics.mean(rr[:,i][np.where(ca_r==12)]))

  np.savetxt(os.path.join(graddir, 'ca/intra', path_list[n]), 
             np.array(intra).T, delimiter = ',')

  # inter-hemisphere  
  inter = [None] * 3  
  for i in range(3):
    inter[i] = (statistics.mean(lr[:,i][np.where(ca_l==1)])-statistics.mean(rl[:,i][np.where(ca_r==1)]),
                statistics.mean(lr[:,i][np.where(ca_l==2)])-statistics.mean(rl[:,i][np.where(ca_r==2)]),
                statistics.mean(lr[:,i][np.where(ca_l==3)])-statistics.mean(rl[:,i][np.where(ca_r==3)]),
                statistics.mean(lr[:,i][np.where(ca_l==4)])-statistics.mean(rl[:,i][np.where(ca_r==4)]),
                statistics.mean(lr[:,i][np.where(ca_l==5)])-statistics.mean(rl[:,i][np.where(ca_r==5)]),
                statistics.mean(lr[:,i][np.where(ca_l==6)])-statistics.mean(rl[:,i][np.where(ca_r==6)]),
                statistics.mean(lr[:,i][np.where(ca_l==7)])-statistics.mean(rl[:,i][np.where(ca_r==7)]),
                statistics.mean(lr[:,i][np.where(ca_l==8)])-statistics.mean(rl[:,i][np.where(ca_r==8)]),
                statistics.mean(lr[:,i][np.where(ca_l==9)])-statistics.mean(rl[:,i][np.where(ca_r==9)]),
                statistics.mean(lr[:,i][np.where(ca_l==10)])-statistics.mean(rl[:,i][np.where(ca_r==10)]),
                statistics.mean(lr[:,i][np.where(ca_l==11)])-statistics.mean(rl[:,i][np.where(ca_r==11)]),
                statistics.mean(lr[:,i][np.where(ca_l==12)])-statistics.mean(rl[:,i][np.where(ca_r==12)]))

  np.savetxt(os.path.join(graddir, 'ca/inter', path_list[n]), 
             np.array(inter).T, delimiter = ',')

# Stats: implement one-sample t-test on intra-hemispheric (LL, RR) and 
# inter-hemispheric (LR, RL) gradient score differences
# for each parcel (1,...,180) and for each gradient (1,2,3) separately
 
intra_g1 = [None] * len(path_list) # diff between LL and RR gradient 1
intra_g2 = [None] * len(path_list)
intra_g3 = [None] * len(path_list)
inter_g1 = [None] * len(path_list) # diff between LR and RL gradient 1
inter_g2 = [None] * len(path_list)
inter_g3 = [None] * len(path_list)    


for n in range(len(path_list)):
  intra_g1[n] = np.array(pd.read_csv(os.path.join(graddir, 'intra', path_list[n]), header=None))[:,0]
  intra_g2[n] = np.array(pd.read_csv(os.path.join(graddir, 'intra', path_list[n]), header=None))[:,1]
  intra_g3[n] = np.array(pd.read_csv(os.path.join(graddir, 'intra', path_list[n]), header=None))[:,2]
  inter_g1[n] = np.array(pd.read_csv(os.path.join(graddir, 'inter', path_list[n]), header=None))[:,0]
  inter_g2[n] = np.array(pd.read_csv(os.path.join(graddir, 'inter', path_list[n]), header=None))[:,1]
  inter_g3[n] = np.array(pd.read_csv(os.path.join(graddir, 'inter', path_list[n]), header=None))[:,2]

intra_g1_stats = [None] * 180
intra_g2_stats = [None] * 180
intra_g3_stats = [None] * 180
inter_g1_stats = [None] * 180
inter_g2_stats = [None] * 180
inter_g3_stats = [None] * 180

for i in range(180):
    intra_g1_stats[i] = stats.ttest_1samp(np.array(intra_g1)[:,i],0)
    intra_g2_stats[i] = stats.ttest_1samp(np.array(intra_g2)[:,i],0)
    intra_g3_stats[i] = stats.ttest_1samp(np.array(intra_g3)[:,i],0)
    inter_g1_stats[i] = stats.ttest_1samp(np.array(inter_g1)[:,i],0)
    inter_g2_stats[i] = stats.ttest_1samp(np.array(inter_g2)[:,i],0)
    inter_g3_stats[i] = stats.ttest_1samp(np.array(inter_g3)[:,i],0)

pd.DataFrame(np.array(intra_g1_stats)).to_csv(os.path.join(graddir, 'intra_g1_stats.csv' ))
pd.DataFrame(np.array(intra_g2_stats)).to_csv(os.path.join(graddir, 'intra_g2_stats.csv' ))
pd.DataFrame(np.array(intra_g3_stats)).to_csv(os.path.join(graddir, 'intra_g3_stats.csv' ))
pd.DataFrame(np.array(inter_g1_stats)).to_csv(os.path.join(graddir, 'inter_g1_stats.csv' ))
pd.DataFrame(np.array(inter_g2_stats)).to_csv(os.path.join(graddir, 'inter_g2_stats.csv' ))
pd.DataFrame(np.array(inter_g3_stats)).to_csv(os.path.join(graddir, 'inter_g3_stats.csv' ))

# FDR correction on the p-values
def fdr(p_vals):
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

pd.DataFrame(np.vstack((fdr(np.array(intra_g1_stats)[:,1]),
                        fdr(np.array(intra_g2_stats)[:,1]),
                        fdr(np.array(intra_g3_stats)[:,1]),
                        fdr(np.array(inter_g1_stats)[:,1]),
                        fdr(np.array(inter_g2_stats)[:,1]),
                        fdr(np.array(inter_g3_stats)[:,1]))).T).\
                        to_csv(os.path.join(graddir, 'g_stats_fdr.csv'),  
                               header = ['intra_g1_fdr','intra_g2_fdr',
                                         'intra_g3_fdr','inter_g1_fdr',
                                         'inter_g2_fdr','inter_g3_fdr'],
                                         index=None)
