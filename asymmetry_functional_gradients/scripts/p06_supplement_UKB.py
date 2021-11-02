"""
calculates the rs-fMRI functional connectivity (FC) matrix for each UKB subject
saves the intra- and inter-hemispheric FC patterns separately
then calculate the gradients for intra- and inter-hemispheric FC patterns, 
correct them, and compute Asymmetry Index, finally doing Cole-Anticevic network
CA 12 networks comparisons with 0.
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
from numpy import inf
from nilearn.connectome import ConnectivityMeasure
from brainspace.gradient import GradientMaps

codedir = os.path.dirname(os.path.abspath(__file__))

fcdir = os.path.join(os.path.dirname(codedir), 
                     'data/data_UKB_glasser_fc')

outdir = os.path.join(os.path.dirname(codedir), 
                     'data/data_results/supplementary/ukb')

# read-in UK-Biobank resting-state fMRI time series for each subject (n=35143),
# and compute the functional connectivity for all 4 fashions:
# LL, RR, LR, RL
path = os.path.join(os.path.dirname(codedir), 
                    'data/BPFiltered')
path_list = os.listdir(path)
path_list.sort()

correlation_measure = ConnectivityMeasure(kind='correlation')
for i in range(len(path_list)):
  ts = np.array(pd.read_csv(path+path_list[i],header=None))[16:].T # first 16 ROIs are subcortical
  if ts.shape[1]==360:
    corr_matrix = correlation_measure.fit_transform([ts])[0]
    fc = np.arctanh(corr_matrix)
    fc[fc == inf] = 0
    np.savetxt((fcdir + '/' + '%.10s' % path_list[i] + '.csv'), fc, delimiter = ',')
    fc_LL = fc[0:180,0:180]
    fc_RR = fc[180:360,180:360]
    fc_LLRR = fc_LL - fc_RR
    fc_LR = fc[0:180,180:360]
    fc_RL = fc[180:360,0:180]
    fc_LRRL = fc_LR - fc_RL
    np.savetxt(outdir + '/FC/LL/'+'%.10s' % path_list[i]+'.csv', fc_LL, delimiter = ',')
    np.savetxt(outdir + '/FC/RR/'+'%.10s' % path_list[i]+'.csv', fc_RR, delimiter = ',')
    np.savetxt(outdir + '/FC/intra/'+'%.10s' % path_list[i]+'.csv', fc_LLRR, delimiter = ',')
    np.savetxt(outdir + '/FC/LR/'+'%.10s' % path_list[i]+'.csv', fc_LR, delimiter = ',')
    np.savetxt(outdir + '/FC/RL/'+'%.10s' % path_list[i]+'.csv', fc_RL, delimiter = ',')
    np.savetxt(outdir + '/FC/inter/'+'%.10s' % path_list[i]+'.csv', fc_LRRL, delimiter = ',')
    print('finish......'+'%.10s' % path_list[i], i)

## mean FC matrix
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

path_list = os.listdir(os.path.join(outdir, 'FC/LL/'))


for i in range(n):
  matrix_fc_LL[i] = np.array(pd.read_csv(os.path.join(outdir, 'FC/LL/', 
                                          path_list[i]),header=None))
  total_fc_LL += matrix_fc_LL[i]
  matrix_fc_RR[i] = np.array(pd.read_csv(os.path.join(outdir, 'FC/RR/', 
                                          path_list[i]),header=None))
  total_fc_RR += matrix_fc_RR[i]
  matrix_fc_LR[i] = np.array(pd.read_csv(os.path.join(outdir, 'FC/LR/', 
                                          path_list[i]),header=None))
  total_fc_LR += matrix_fc_LR[i]
  matrix_fc_RL[i] = np.array(pd.read_csv(os.path.join(outdir, 'FC/RL/', 
                                          path_list[i]),header=None))
  total_fc_RL += matrix_fc_RL[i]
  matrix_fc_LLRR[i] = np.array(pd.read_csv(os.path.join(outdir, 'FC/intra/', 
                                            path_list[i]),header=None))
  total_fc_LLRR += matrix_fc_LLRR[i]
  matrix_fc_LRRL[i] = np.array(pd.read_csv(os.path.join(outdir, 'FC/inter/', 
                                            path_list[i]),header=None))
  total_fc_LRRL += matrix_fc_LRRL[i]

mean_fc_LL = total_fc_LL/n
mean_fc_RR = total_fc_RR/n
mean_fc_RL = total_fc_RL/n
mean_fc_LR = total_fc_LR/n
mean_fc_LLRR = total_fc_LLRR/n
mean_fc_LRRL = total_fc_LRRL/n
np.savetxt(os.path.join(outdir, 'FC/LL_groupmean.csv'), mean_fc_LL, delimiter = ',')
np.savetxt(os.path.join(outdir, 'FC/RR_groupmean.csv'), mean_fc_RR, delimiter = ',')
np.savetxt(os.path.join(outdir, 'FC/LR_groupmean.csv'), mean_fc_LR, delimiter = ',')
np.savetxt(os.path.join(outdir, 'FC/RL_groupmean.csv'), mean_fc_RL, delimiter = ',')
np.savetxt(os.path.join(outdir, 'FC/intra_groupmean.csv'), mean_fc_LLRR, delimiter = ',')
np.savetxt(os.path.join(outdir, 'FC/inter_groupmean.csv'), mean_fc_LRRL, delimiter = ',')

# compute the gradients of functional connectivity for the mean LL fc
# -> this will be our reference gradients for the UK-Biobank data
np.random.seed(0)
gm = GradientMaps(approach='dm', kernel='normalized_angle',
                  n_components=10,random_state=0)
LL = np.array(pd.read_csv(os.path.join(outdir, 'FC/LL_groupmean.csv'), header=None))
gm.fit(LL)
group_grad_LL = gm.gradients_
print(group_grad_LL)

np.savetxt(os.path.join(outdir, 'gradient/group_grad_LL.csv'), group_grad_LL, delimiter = ',')
np.savetxt(os.path.join(outdir, 'gradient/group_grad_LL_lambdas.csv'), gm.lambdas_, delimiter = ',')

# get the gradients for each subject and each FC (LL, RR, LR, RL) and
# align them to the group-level LL gradients
path_list = os.listdir(os.path.join(outdir, 'FC/LL/'))
n = len(path_list)

for i in range(n):
  np.random.seed(0)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LL = np.array(pd.read_csv(os.path.join(outdir, 'FC/LL', path_list[i]),
                                            header=None))
  align.fit(fc_LL, reference=group_grad_LL)
  grad_LL = align.gradients_
  np.savetxt(os.path.join(outdir, 'gradient/LL', path_list[i]), 
                          grad_LL, delimiter = ',')
  np.random.seed(0)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  fc_RR = np.array(pd.read_csv(os.path.join(outdir, 'FC/RR', path_list[i]), 
                                            header=None))
  align.fit(fc_RR,reference=group_grad_LL)  
  grad_RR = align.gradients_
  np.savetxt(os.path.join(outdir, 'gradient/RR', path_list[i]), 
                          grad_RR, delimiter = ',')
  np.random.seed(0)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LR = np.array(pd.read_csv(os.path.join(outdir, 'FC/LR', path_list[i]), 
                                            header=None))
  align.fit(fc_LR,reference=group_grad_LL)
  grad_LR = align.gradients_
  np.savetxt(os.path.join(outdir, 'gradient/LR', path_list[i]), 
                         grad_LR, delimiter = ',')
  np.random.seed(0)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  fc_RL = np.array(pd.read_csv(os.path.join(outdir, 'FC/RL', path_list[i]), 
                                            header=None))
  align.fit(fc_RL,reference=group_grad_LL)
  grad_RL = align.gradients_
  np.savetxt(os.path.join(outdir, 'gradient/RL', path_list[i]),
                          grad_RL, delimiter = ',')
  print(i)

# quality check: get the correlations between individual gradients and 
# the group-level template gradient, swap the axis if negative,
# and compute the Asymmetry Index for intra- and inter-hemisphere gradients
group_grad_LL = np.array(pd.read_csv(os.path.join(outdir, 
                                   'gradient/group_grad_LL.csv'), header=None))
total_llrr = 0
total_lrrl = 0

cons = ['LL', 'RR', 'LR', 'RL']
for dir in path_list:
    for con in cons:
      gradFname = os.path.join(outdir, 'gradient', con, dir)    
      df = np.array(pd.read_csv(gradFname, header=None))
      r = [None] * 10
      corrected = [None]*10
      for i in range(10):
        r[i] = stats.pearsonr(group_grad_LL[:,i], df[:,i])    
        if r[i][0] > 0:
            corrected[i]=df[:,i]
        else:
            corrected[i]=-1*df[:,i]

      correctedX = np.array(corrected).T
      np.savetxt(gradFname, correctedX, delimiter = ',')
      
      if con == 'LL':
          corrected_ll = correctedX
      elif con == 'RR':
          corrected_rr = correctedX
      elif con == 'LR':
          corrected_lr = correctedX
      elif con == 'RL':
          corrected_rl = correctedX

    # RR-LL, RL-RL
    AI_llrr = corrected_ll - corrected_rr
    AI_lrrl = corrected_lr - corrected_rl    
    np.savetxt(os.path.join(outdir, 'gradient/intra', dir), AI_llrr, delimiter = ',')
    np.savetxt(os.path.join(outdir, 'gradient/inter', dir), AI_lrrl, delimiter = ',')
    print('finish   ' + dir)
    # mean AI
    total_llrr = total_llrr + AI_llrr
    total_lrrl = total_lrrl + AI_lrrl

mean_llrr = total_llrr/len(path_list)
mean_lrrl = total_lrrl/len(path_list)
np.savetxt(os.path.join(outdir, 'gradient/mean_asym_intra.csv'), mean_llrr, delimiter = ',')
np.savetxt(os.path.join(outdir, 'gradient/mean_asym_inter.csv'), mean_lrrl, delimiter = ',')

# Stats: implement one-sample t-test on Asymmetry Index scores, 
# for intra-hemispheric (LL, RR) and inter-hemispheric (LR, RL) fashions,
# for each parcel (1,...,180), and for each gradient (1,2,3), separately
path = os.path.join(outdir, 'gradient/intra')
path_list = os.listdir(path)
path_list.sort()

intra_g1 = [None] * len(path_list)
intra_g2 = [None] * len(path_list)
intra_g3 = [None] * len(path_list)
inter_g1 = [None] * len(path_list)
inter_g2 = [None] * len(path_list)
inter_g3 = [None] * len(path_list)
for n in range(len(path_list)):

    intra_g1[n] = np.array(pd.read_csv(os.path.join(outdir, 'gradient/intra', 
            path_list[n]),header=None))[:,0]
    intra_g2[n] = np.array(pd.read_csv(os.path.join(outdir, 'gradient/intra', 
            path_list[n]),header=None))[:,1]
    intra_g3[n] = np.array(pd.read_csv(os.path.join(outdir, 'gradient/intra', 
            path_list[n]),header=None))[:,2]
    inter_g1[n] = np.array(pd.read_csv(os.path.join(outdir, 'gradient/inter', 
            path_list[n]),header=None))[:,0]
    inter_g2[n] = np.array(pd.read_csv(os.path.join(outdir, 'gradient/inter', 
            path_list[n]),header=None))[:,1]
    inter_g3[n] = np.array(pd.read_csv(os.path.join(outdir, 'gradient/inter', 
            path_list[n]),header=None))[:,2]
    print(n)
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
    
pd.DataFrame(np.array(intra_g1_stats)).to_csv(os.path.join(outdir, 'gradient/intra_g1_stats.csv'))
pd.DataFrame(np.array(intra_g2_stats)).to_csv(os.path.join(outdir, 'gradient/intra_g2_stats.csv'))
pd.DataFrame(np.array(intra_g3_stats)).to_csv(os.path.join(outdir, 'gradient/intra_g3_stats.csv'))
pd.DataFrame(np.array(inter_g1_stats)).to_csv(os.path.join(outdir, 'gradient/inter_g1_stats.csv'))
pd.DataFrame(np.array(inter_g2_stats)).to_csv(os.path.join(outdir, 'gradient/inter_g2_stats.csv'))
pd.DataFrame(np.array(inter_g3_stats)).to_csv(os.path.join(outdir, 'gradient/inter_g3_stats.csv'))

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
                        to_csv(os.path.join(outdir, 
                                            'gradient/g_stats_fdr.csv'),
                        header = ['intra_g1_fdr','intra_g2_fdr','intra_g3_fdr',
                                  'inter_g1_fdr','inter_g2_fdr','inter_g3_fdr'],
                                  index=None)

# read Cole-Anticevic network parcels for left and right hemispheres,
# for each subject, we parcellate the gradients into 12 networks,
# for intra-hemispheric (LL, RR) and inter-hemispheric (LR, RL) differences, 
# we subtract the mean gradient scores in each network
ca_l = np.array(pd.read_csv('../data/ca_glasser_network.csv',header=None))[:,0][:180]
ca_r = np.array(pd.read_csv('../data/ca_glasser_network.csv',header=None))[:,0][180:]

for n in range(len(path_list)):
  ll = np.array(pd.read_csv(os.path.join(outdir, 'gradient/LL', path_list[n]), header=None))
  rr = np.array(pd.read_csv(os.path.join(outdir, 'gradient/RR', path_list[n]), header=None))
  lr = np.array(pd.read_csv(os.path.join(outdir, 'gradient/LR', path_list[n]), header=None))
  rl = np.array(pd.read_csv(os.path.join(outdir, 'gradient/RL', path_list[n]), header=None))
  intra = [None] * 3  
  for i in range(3):
    intra[i] = (np.mean(ll[:,i][ca_l==1])-np.mean(rr[:,i][ca_r==1]),
                np.mean(ll[:,i][ca_l==2])-np.mean(rr[:,i][ca_r==2]),
                np.mean(ll[:,i][ca_l==3])-np.mean(rr[:,i][ca_r==3]),
                np.mean(ll[:,i][ca_l==4])-np.mean(rr[:,i][ca_r==4]),
                np.mean(ll[:,i][ca_l==5])-np.mean(rr[:,i][ca_r==5]),
                np.mean(ll[:,i][ca_l==6])-np.mean(rr[:,i][ca_r==6]),
                np.mean(ll[:,i][ca_l==7])-np.mean(rr[:,i][ca_r==7]),
                np.mean(ll[:,i][ca_l==8])-np.mean(rr[:,i][ca_r==8]),
                np.mean(ll[:,i][ca_l==9])-np.mean(rr[:,i][ca_r==9]),
                np.mean(ll[:,i][ca_l==10])-np.mean(rr[:,i][ca_r==10]),
                np.mean(ll[:,i][ca_l==11])-np.mean(rr[:,i][ca_r==11]),
                np.mean(ll[:,i][ca_l==12])-np.mean(rr[:,i][ca_r==12]))
  
  np.savetxt(os.path.join(outdir, 'gradient/network/intra', path_list[n]),
             np.array(intra).T, delimiter = ',')
    
  inter = [None] * 3  
  for i in range(3):
    inter[i] = (np.mean(lr[:,i][ca_l==1])-np.mean(rl[:,i][ca_r==1]),
                np.mean(lr[:,i][ca_l==2])-np.mean(rl[:,i][ca_r==2]),
                np.mean(lr[:,i][ca_l==3])-np.mean(rl[:,i][ca_r==3]),
                np.mean(lr[:,i][ca_l==4])-np.mean(rl[:,i][ca_r==4]),
                np.mean(lr[:,i][ca_l==5])-np.mean(rl[:,i][ca_r==5]),
                np.mean(lr[:,i][ca_l==6])-np.mean(rl[:,i][ca_r==6]),
                np.mean(lr[:,i][ca_l==7])-np.mean(rl[:,i][ca_r==7]),
                np.mean(lr[:,i][ca_l==8])-np.mean(rl[:,i][ca_r==8]),
                np.mean(lr[:,i][ca_l==9])-np.mean(rl[:,i][ca_r==9]),
                np.mean(lr[:,i][ca_l==10])-np.mean(rl[:,i][ca_r==10]),
                np.mean(lr[:,i][ca_l==11])-np.mean(rl[:,i][ca_r==11]),
                np.mean(lr[:,i][ca_l==12])-np.mean(rl[:,i][ca_r==12]))
  
  np.savetxt(os.path.join(outdir, 'gradient/network/inter', path_list[n]),
             np.array(inter).T, delimiter = ',')  
  print('finish......'+path_list[n])

