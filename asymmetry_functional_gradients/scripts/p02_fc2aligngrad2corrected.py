"""
computes the gradients of functional connectivity for 4 different fashions:
LL & RR intra-hemispheric and LR & RL inter-hemispheric fashion.
using Cole-Anticevic network parcellations, we compute the gradient scores
in each network and for all 4 fashions. we implement one sample t-test
on intra-hemispheric and inter-hemispheric gradient score differences in
each network (LL-RR differences and LR-RL differences).
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
from brainspace.gradient import GradientMaps
import statistics

codedir = os.path.dirname(os.path.abspath(__file__))

datadir = os.path.join(os.path.dirname(codedir), 
                       'data/data_results/FC/')

graddir = os.path.join(os.path.dirname(codedir), 
                       'data/data_results/gradient/')

path = os.path.join(datadir, 'LL')
path_list = os.listdir(path)
path_list.sort()

# compute the group-level functional connectivity (FC) matrix 
# for intra-hemispheric (LL, RR) and inter-hemispheric (LR, RL) fashions
n = len(path_list)
matrix_fc_LL = [None] * len(path_list)
matrix_fc_RR = [None] * len(path_list)
matrix_fc_RL = [None] * len(path_list)
matrix_fc_LR = [None] * len(path_list)
matrix_fc_LLRR = [None] * len(path_list)
matrix_fc_LRRL = [None] * len(path_list)
total_fc_LL = 0
total_fc_RR = 0
total_fc_RL = 0
total_fc_LR = 0
total_fc_LLRR = 0
total_fc_LRRL = 0
total_llrr = 0
total_lrrl = 0

for i in range(n):
  matrix_fc_LL[i] = np.array(pd.read_csv(os.path.join(datadir, 'LL', 
                              path_list[i]), header=None))
  total_fc_LL += matrix_fc_LL[i]
  matrix_fc_RR[i] = np.array(pd.read_csv(os.path.join(datadir, 'RR', 
                              path_list[i]), header=None))
  total_fc_RR += matrix_fc_RR[i]
  matrix_fc_LR[i] = np.array(pd.read_csv(os.path.join(datadir, 'LR', 
                              path_list[i]), header=None))
  total_fc_LR += matrix_fc_LR[i]
  matrix_fc_RL[i] = np.array(pd.read_csv(os.path.join(datadir, 'RL', 
                              path_list[i]), header=None))
  total_fc_RL += matrix_fc_RL[i]
  matrix_fc_LLRR[i] = np.array(pd.read_csv(os.path.join(datadir, 'LL-RR', 
                                path_list[i]), header=None))
  total_fc_LLRR += matrix_fc_LLRR[i]
  matrix_fc_LRRL[i] = np.array(pd.read_csv(os.path.join(datadir, 'LR-RL', 
                                path_list[i]), header=None))
  total_fc_LRRL += matrix_fc_LRRL[i]
  
mean_fc_LL = total_fc_LL/n
mean_fc_RR = total_fc_RR/n
mean_fc_RL = total_fc_RL/n
mean_fc_LR = total_fc_LR/n
mean_fc_LLRR = total_fc_LLRR/n
mean_fc_LRRL = total_fc_LRRL/n

np.savetxt(os.path.join(datadir, 'LL_groupmean.csv'), mean_fc_LL, 
           delimiter = ',')
np.savetxt(os.path.join(datadir, 'RR_groupmean.csv'), mean_fc_RR, 
           delimiter = ',')
np.savetxt(os.path.join(datadir, 'LR_groupmean.csv'), mean_fc_LR, 
           delimiter = ',')
np.savetxt(os.path.join(datadir, 'RL_groupmean.csv'), mean_fc_RL, 
           delimiter = ',')
np.savetxt(os.path.join(datadir, 'LLRR_groupmean.csv'), mean_fc_LLRR, 
           delimiter = ',')
np.savetxt(os.path.join(datadir, 'LRRL_groupmean.csv'), mean_fc_LRRL, 
           delimiter = ',')


print('Done! fc matrix computation across all subjects finished...')

# get the gradients of group-level FC for left-left (LL) connections
gm = GradientMaps(n_components=10, random_state=0, approach='dm', 
                  kernel='normalized_angle')
LL = np.array(pd.read_csv(os.path.join(datadir, 'LL_groupmean.csv'), 
                          header=None))
gm.fit(LL)
group_grad_LL = gm.gradients_

np.savetxt(os.path.join(graddir, 'group_grad_LL.csv'), group_grad_LL, 
           delimiter = ',')
np.savetxt(os.path.join(graddir, 'group_grad_LL_lambdas.csv'), gm.lambdas_, 
           delimiter = ',')

# get the gradients for each subject and each FC (LL, RR, LR, RL) and
# align them to the group-level LL gradients
for i in path_list:
  # FC LL
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LL = np.array(pd.read_csv(os.path.join(datadir, 'LL', i), header=None))
  align.fit(fc_LL, reference=group_grad_LL)
  grad_LL = align.gradients_
  np.savetxt(os.path.join(graddir, 'LL', i), grad_LL, delimiter = ',')

  # FC RR
  align = GradientMaps(n_components=10, random_state=0, approach='dm',
                       kernel='normalized_angle', alignment='procrustes')
  fc_RR = np.array(pd.read_csv(os.path.join(datadir, 'RR', i), header=None))
  align.fit(fc_RR,reference=group_grad_LL)
  grad_RR = align.gradients_
  np.savetxt(os.path.join(graddir, 'RR', i), grad_RR, delimiter = ',')  

  # FC LR 
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  fc_LR = np.array(pd.read_csv(os.path.join(datadir, 'LR', i), header=None))
  align.fit(fc_LR,reference=group_grad_LL)
  grad_LR = align.gradients_
  np.savetxt(os.path.join(graddir, 'LR', i), grad_LR, delimiter = ',')
  
  # FC RL  
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  
  fc_RL = np.array(pd.read_csv(os.path.join(datadir, 'RL', i), header=None))
  align.fit(fc_RL,reference=group_grad_LL)
  grad_RL = align.gradients_
  np.savetxt(os.path.join(graddir, 'RL', i), grad_RL, delimiter = ',')
  print('finish   ' + i)
  
# quality check: get the correlations between individual gradients and 
# the group-level template gradient, swap the axis if negative
cons = ['LL', 'RR', 'LR', 'RL']

for dir in path_list:
  for con in cons:     
      print(dir, con)
      df = np.array(pd.read_csv(os.path.join(graddir, con, dir), header=None))
      r = [None] * 10
      corrected = [None]*10
      for i in range(10):
        r[i] = stats.pearsonr(group_grad_LL[:,i],df[:,i])    
        if r[i][0] > 0:
            corrected[i]=df[:,i]
        else:
            corrected[i]=-1*df[:,i]
            m = i +1
            correct = open(os.path.join(graddir, 'corrected.txt'), 'a')
            correct.write(str(con) + '  ' + dir + '  ' + str(m) +'\n')
            correct.close()
      corrected_x = np.concatenate((corrected[0].reshape(corrected[0].shape[0],1),
                                    corrected[1].reshape(corrected[1].shape[0],1),
                                    corrected[2].reshape(corrected[2].shape[0],1),
                                    corrected[3].reshape(corrected[3].shape[0],1),
                                    corrected[4].reshape(corrected[4].shape[0],1),
                                    corrected[5].reshape(corrected[5].shape[0],1),
                                    corrected[6].reshape(corrected[6].shape[0],1),
                                    corrected[7].reshape(corrected[7].shape[0],1),
                                    corrected[8].reshape(corrected[8].shape[0],1),
                                    corrected[9].reshape(corrected[9].shape[0],1))
                                    ,axis = 1)
      np.savetxt(os.path.join(graddir, con, dir), corrected_x, delimiter=',')

      if con == 'LL':
          corrected_ll = corrected_x
      elif con == 'RR':
          corrected_rr = corrected_x
      elif con == 'LR':
          corrected_lr = corrected_x
      elif con == 'RL':
          corrected_rl = corrected_x
          
  AI_llrr = corrected_ll - corrected_rr
  AI_lrrl = corrected_lr - corrected_rl
  np.savetxt(os.path.join(graddir, 'LL-RR', dir), AI_llrr, delimiter = ',')
  np.savetxt(os.path.join(graddir, 'LR-RL', dir), AI_lrrl, delimiter = ',')     
  # mean AI
  total_llrr = total_llrr + AI_llrr
  total_lrrl = total_lrrl + AI_lrrl
  
mean_llrr = total_llrr/len(path_list)
mean_lrrl = total_lrrl/len(path_list)

np.savetxt(os.path.join(graddir, 'mean_asym_LLRR.csv'), mean_llrr, 
           delimiter = ',')
np.savetxt(os.path.join(graddir, 'mean_asym_LRRL.csv'), mean_lrrl, 
           delimiter = ',')

# read Cole-Anticevic network parcels for left and right hemispheres,
# for each subject, we parcellate the gradients into 12 networks,
# for intra-hemispheric (LL, RR) and inter-hemispheric (LR, RL) differences, 
# we subtract the mean gradient scores in each network

cadir= os.path.join(os.path.dirname(codedir), 'data')

ca_l = np.array(pd.read_csv(os.path.join(cadir, 'ca_glasser_network.csv'),
                                         header=None))[:,0][:180]
ca_r = np.array(pd.read_csv(os.path.join(cadir, 'ca_glasser_network.csv'),
                                         header=None))[:,0][180:]

for n in range(len(path_list)):
  ll = np.array(pd.read_csv(os.path.join(graddir,'LL',path_list[n]),header=None))
  rr = np.array(pd.read_csv(os.path.join(graddir,'RR',path_list[n]),header=None))
  lr = np.array(pd.read_csv(os.path.join(graddir,'LR',path_list[n]),header=None))
  rl = np.array(pd.read_csv(os.path.join(graddir,'RL',path_list[n]),header=None))
  
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

  np.savetxt(os.path.join(graddir, 'intra_ca', path_list[n]), 
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

  np.savetxt(os.path.join(graddir, 'inter_ca', path_list[n]), 
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
  intra_g1[n] = np.array(pd.read_csv(os.path.join(graddir, 'LL-RR', path_list[n]), header=None))[:,0]
  intra_g2[n] = np.array(pd.read_csv(os.path.join(graddir, 'LL-RR', path_list[n]), header=None))[:,1]
  intra_g3[n] = np.array(pd.read_csv(os.path.join(graddir, 'LL-RR', path_list[n]), header=None))[:,2]
  inter_g1[n] = np.array(pd.read_csv(os.path.join(graddir, 'LR-RL', path_list[n]), header=None))[:,0]
  inter_g2[n] = np.array(pd.read_csv(os.path.join(graddir, 'LR-RL', path_list[n]), header=None))[:,1]
  inter_g3[n] = np.array(pd.read_csv(os.path.join(graddir, 'LR-RL', path_list[n]), header=None))[:,2]

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

