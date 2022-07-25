"""
This analysis script includes two parts: 1). Normalize AI (asymmetry index) formula 
and 2). compute FC gradients using the DK atlas.

Normalize AI formula: (l-r)(l+r): to aviod l+r becomeing 0, we change the scale
using +1 for all gradients values to make them all positive. Finally to investigate the
difference of using the AI formula, used in the main manuscript, and using normalized AI formula.

DK atlas: computes the gradients of functional connectivity (LL & RR intra-hemispheric 
and LR & RL inter-hemispheric) using the DK atlas, after we compute the gradient scores 
we implement one sample t-test on intra-hemispheric and inter-hemispheric gradient 
score (LL-RR differences and LR-RL differences).
"""
import os
import pandas as pd
import numpy as np
from scipy import stats

codedir = os.path.dirname(os.path.abspath(__file__))

grdir = os.path.join(os.path.dirname(codedir), 
                     'data/data_results/gradient')

outdir = os.path.join(os.path.dirname(codedir), 
                      'data/data_results/supplementary')

path = os.path.join(grdir, 'LL') 
path_list = os.listdir(path)
path_list.sort()

n = len(path_list)

# compute the Asymmetry Index using another formula on the HCP gradients
# for all 4 fashions: LL, RR, LR, RL connectivity spaces
for i in range(n):
  ll = np.array(pd.read_csv(os.path.join(grdir, 'LL', path_list[i]), header=None)) + 1
  rr = np.array(pd.read_csv(os.path.join(grdir, 'RR', path_list[i]), header=None)) + 1
  lr = np.array(pd.read_csv(os.path.join(grdir, 'LR', path_list[i]), header=None)) + 1
  rl = np.array(pd.read_csv(os.path.join(grdir, 'RL', path_list[i]), header=None)) + 1
  intra = (ll-rr)/(ll+rr)
  inter = (lr-rl)/(lr+rl)
  np.savetxt(os.path.join(outdir, 'normalize/intra', path_list[i]), intra, 
             delimiter = ',')
  np.savetxt(os.path.join(outdir, 'normalize/inter', path_list[i]), inter, 
             delimiter = ',')

# read-in the Asymmetry Index scores along the first three gradients
intra_g1 = [None] * n
intra_g2 = [None] * n
intra_g3 = [None] * n
inter_g1 = [None] * n
inter_g2 = [None] * n
inter_g3 = [None] * n
for n in range(len(path_list)):
        intra_g1[n] = np.array(pd.read_csv(os.path.join(outdir, 
                'normalize/intra', path_list[n]),header=None))[:,0]
        intra_g2[n] = np.array(pd.read_csv(os.path.join(outdir, 
                'normalize/intra', path_list[n]),header=None))[:,1]
        intra_g3[n] = np.array(pd.read_csv(os.path.join(outdir, 
                'normalize/intra', path_list[n]),header=None))[:,2]
        inter_g1[n] = np.array(pd.read_csv(os.path.join(outdir, 
                'normalize/inter', path_list[n]),header=None))[:,0]
        inter_g2[n] = np.array(pd.read_csv(os.path.join(outdir, 
                'normalize/inter', path_list[n]),header=None))[:,1]
        inter_g3[n] = np.array(pd.read_csv(os.path.join(outdir, 
                'normalize/intra', path_list[n]),header=None))[:,2]
        
# Stats: implement one-sample t-test on Asymmetry Index scores, 
# for intra-hemispheric (LL, RR) and inter-hemispheric (LR, RL) fashions,
# for each parcel (1,...,180), and for each gradient (1,2,3), separately
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

print(np.shape(intra_g1_stats))

pd.DataFrame(np.array(intra_g1_stats)).to_csv(os.path.join(outdir, 
                'normalize', 'intra_g1_stats.csv'))
pd.DataFrame(np.array(intra_g2_stats)).to_csv(os.path.join(outdir, 
                'normalize', 'intra_g2_stats.csv'))
pd.DataFrame(np.array(intra_g3_stats)).to_csv(os.path.join(outdir, 
                'normalize', 'intra_g3_stats.csv'))
pd.DataFrame(np.array(inter_g1_stats)).to_csv(os.path.join(outdir, 
                'normalize', 'inter_g1_stats.csv')) 
pd.DataFrame(np.array(inter_g2_stats)).to_csv(os.path.join(outdir, 
                'normalize', 'inter_g2_stats.csv'))
pd.DataFrame(np.array(inter_g3_stats)).to_csv(os.path.join(outdir, 
                'normalize', 'inter_g3_stats.csv'))

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
                        'normalize', 'g_stats_fdr.csv'),header=['intra_g1_fdr',
                        'intra_g2_fdr','intra_g3_fdr','inter_g1_fdr',
                        'inter_g2_fdr','inter_g3_fdr'],index=None)

# Desikan-Killiany Atlas parcellation on HCP subjects (n=1014)
# read-in functional connectivity data and save for all 4 fashions:
# LL, RR, LR, RL, and compute inter- and intra-hem. differences
dkdir = os.path.join(os.path.dirname(codedir), 
                     'data/data_results/supplementary/dk')
path_list = os.listdir(os.path.join(dkdir, 'FC/raw'))
path_list.sort()

for dir in path_list:
    fc = np.array(pd.read_csv(os.path.join(dkdir, 'FC/raw', dir), header = None))
    fc_LL = fc[0:35,0:35]
    np.savetxt(os.path.join(dkdir, 'FC/LL', dir), fc_LL, delimiter = ',')
    fc_RR = fc[35:70,35:70]
    np.savetxt(os.path.join(dkdir, 'FC/RR', dir), fc_RR, delimiter = ',')
    fc_LR = fc[0:35,35:70]
    np.savetxt(os.path.join(dkdir, 'FC/LR', dir), fc_LR, delimiter = ',')
    fc_RL = fc[35:70,0:35]
    np.savetxt(os.path.join(dkdir, 'FC/RL', dir), fc_RL, delimiter = ',')
    fc_LLRR = fc_LL - fc_RR
    np.savetxt(os.path.join(dkdir, 'FC/intra', dir), fc_LLRR, delimiter = ',')    
    fc_LRRL = fc_LR - fc_RL   
    np.savetxt(os.path.join(dkdir, 'FC/inter', dir), fc_LRRL, delimiter = ',')

# get the mean FC matrix (34x34 cortical and 1x1 white matter) across subjects 
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
  matrix_fc_LL[i] = np.array(pd.read_csv(os.path.join(dkdir, 'FC/LL', 
                             path_list[i]), header=None))
  total_fc_LL += matrix_fc_LL[i]
  matrix_fc_RR[i] = np.array(pd.read_csv(os.path.join(dkdir, 'FC/RR', 
                             path_list[i]), header=None))
  total_fc_RR += matrix_fc_RR[i]
  matrix_fc_LR[i] = np.array(pd.read_csv(os.path.join(dkdir, 'FC/LR', 
                             path_list[i]), header=None))
  total_fc_LR += matrix_fc_LR[i]
  matrix_fc_RL[i] = np.array(pd.read_csv(os.path.join(dkdir, 'FC/RL', 
                             path_list[i]), header=None))
  total_fc_RL += matrix_fc_RL[i]
  matrix_fc_LLRR[i] = np.array(pd.read_csv(os.path.join(dkdir, 'FC/intra', 
                             path_list[i]), header=None))
  total_fc_LLRR += matrix_fc_LLRR[i]
  matrix_fc_LRRL[i] = np.array(pd.read_csv(os.path.join(dkdir, 'FC/inter', 
                             path_list[i]), header=None))
  total_fc_LRRL += matrix_fc_LRRL[i]

mean_fc_LL = total_fc_LL/n
mean_fc_RR = total_fc_RR/n
mean_fc_RL = total_fc_RL/n
mean_fc_LR = total_fc_LR/n
mean_fc_LLRR = total_fc_LLRR/n
mean_fc_LRRL = total_fc_LRRL/n
np.savetxt(os.path.join(dkdir, 'FC/LL_groupmean.csv'), mean_fc_LL, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'FC/RR_groupmean.csv'), mean_fc_RR, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'FC/LR_groupmean.csv'), mean_fc_LR, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'FC/RL_groupmean.csv'), mean_fc_RL, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'FC/intra_groupmean.csv'), mean_fc_LLRR, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'FC/inter_groupmean.csv'), mean_fc_LRRL, delimiter = ',')

# compute the gradients of functional connectivity for the mean LL fashion
# -> this will be our reference gradient space for the Desikan-Killiany parcel
from brainspace.gradient import GradientMaps
gm = GradientMaps(approach='dm', kernel='normalized_angle',
                  n_components=10,random_state=0)
LL = np.array(pd.read_csv(os.path.join(dkdir, 'FC/LL_groupmean.csv'), 
                          header=None))
gm.fit(LL)
group_grad_LL = gm.gradients_
np.savetxt(os.path.join(dkdir, 'gradient/group_grad_LL.csv'), group_grad_LL, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'gradient/group_grad_LL_lambdas.csv'), gm.lambdas_, delimiter = ',')

# get the individual gradients for each subject and each FC (LL, RR, LR, RL),
# and align them to the group-level LL gradients
for i in path_list:
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LL = np.array(pd.read_csv(os.path.join(dkdir, 'FC/LL', i), header=None))
  align.fit(fc_LL,reference=group_grad_LL)
  grad_LL = align.aligned_
  np.savetxt(os.path.join(dkdir, 'gradient/LL', i), grad_LL, delimiter = ',')
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  fc_RR = np.array(pd.read_csv(os.path.join(dkdir, 'FC/RR', i), header=None))
  align.fit(fc_RR,reference=group_grad_LL)
  grad_RR = align.aligned_
  np.savetxt(os.path.join(dkdir, 'gradient/RR', i), grad_RR, delimiter = ',')
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  fc_LR = np.array(pd.read_csv(os.path.join(dkdir, 'FC/LR', i), header=None))
  align.fit(fc_LR,reference=group_grad_LL)
  grad_LR = align.aligned_
  np.savetxt(os.path.join(dkdir, 'gradient/LR', i), grad_LR, delimiter = ',')
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  fc_RL = np.array(pd.read_csv(os.path.join(dkdir, 'FC/RL', i), header=None))
  align.fit(fc_RL,reference=group_grad_LL)
  grad_RL = align.aligned_
  np.savetxt(os.path.join(dkdir, 'gradient/RL', i), grad_RL, delimiter = ',')

# quality check: get the correlations between individual gradients and 
# the group-level template gradient, swap the axis if negative,
# and compute the Asymmetry Index for intra- and inter-hemisphere gradients
group_grad_LL = np.array(pd.read_csv(os.path.join(dkdir, 
                                'gradient/group_grad_LL.csv'), header=None))
cons = ['LL', 'RR', 'LR', 'RL']
for dir in path_list:
    for con in cons:
      gradFname = os.path.join(dkdir, 'gradient', con, dir)      
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
    np.savetxt(os.path.join(dkdir, 'gradient/intra', dir), AI_llrr, delimiter = ',')
    np.savetxt(os.path.join(dkdir, 'gradient/inter', dir), AI_lrrl, delimiter = ',')
    print('finish   ' + dir)
    # mean AI
    total_llrr = total_llrr + AI_llrr
    total_lrrl = total_lrrl + AI_lrrl

mean_llrr = total_llrr/len(path_list)
mean_lrrl = total_lrrl/len(path_list)
np.savetxt(os.path.join(dkdir, 'gradient/mean_asym_intra.csv'), mean_llrr, delimiter = ',')
np.savetxt(os.path.join(dkdir, 'gradient/mean_asym_inter.csv'), mean_lrrl, delimiter = ',')

# read-in AI scores for each subject, each gradient, intra- and inter-hemiphere 
path = os.path.join(dkdir, 'gradient/intra')
path_list = os.listdir(path)
path_list.sort()

intra_g1 = [None] * len(path_list)
intra_g2 = [None] * len(path_list)
intra_g3 = [None] * len(path_list)
inter_g1 = [None] * len(path_list)
inter_g2 = [None] * len(path_list)
inter_g3 = [None] * len(path_list)
for n in range(len(path_list)):
    intra_g1[n] = np.array(pd.read_csv(os.path.join(dkdir, 'gradient/intra', 
                        path_list[n]), header=None))[:,0]
    intra_g2[n] = np.array(pd.read_csv(os.path.join(dkdir, 'gradient/intra', 
                        path_list[n]), header=None))[:,1]
    intra_g3[n] = np.array(pd.read_csv(os.path.join(dkdir, 'gradient/intra', 
                        path_list[n]), header=None))[:,2]
    inter_g1[n] = np.array(pd.read_csv(os.path.join(dkdir, 'gradient/inter', 
                        path_list[n]), header=None))[:,0]
    inter_g2[n] = np.array(pd.read_csv(os.path.join(dkdir, 'gradient/inter', 
                        path_list[n]), header=None))[:,1]
    inter_g3[n] = np.array(pd.read_csv(os.path.join(dkdir, 'gradient/inter', 
                        path_list[n]), header=None))[:,2]

# Stats: implement one-sample t-test on Asymmetry Index scores, 
# for intra-hemispheric (LL, RR) and inter-hemispheric (LR, RL) fashions,
# for each parcel (1,...,35), and for each gradient (1,2,3), separately
intra_g1_stats = [None] * 35
intra_g2_stats = [None] * 35
intra_g3_stats = [None] * 35
inter_g1_stats = [None] * 35
inter_g2_stats = [None] * 35
inter_g3_stats = [None] * 35

for i in range(35):
    intra_g1_stats[i] = stats.ttest_1samp(np.array(intra_g1)[:,i],0)
    intra_g2_stats[i] = stats.ttest_1samp(np.array(intra_g2)[:,i],0)
    intra_g3_stats[i] = stats.ttest_1samp(np.array(intra_g3)[:,i],0)
    inter_g1_stats[i] = stats.ttest_1samp(np.array(inter_g1)[:,i],0)
    inter_g2_stats[i] = stats.ttest_1samp(np.array(inter_g2)[:,i],0)
    inter_g3_stats[i] = stats.ttest_1samp(np.array(inter_g3)[:,i],0)
    
pd.DataFrame(np.array(intra_g1_stats)).to_csv(os.path.join(dkdir, 
            'gradient/intra_g1_stats.csv'))
pd.DataFrame(np.array(intra_g2_stats)).to_csv(os.path.join(dkdir, 
            'gradient/intra_g2_stats.csv'))
pd.DataFrame(np.array(intra_g3_stats)).to_csv(os.path.join(dkdir, 
            'gradient/intra_g3_stats.csv'))
pd.DataFrame(np.array(inter_g1_stats)).to_csv(os.path.join(dkdir, 
            'gradient/inter_g1_stats.csv'))
pd.DataFrame(np.array(inter_g2_stats)).to_csv(os.path.join(dkdir, 
            'gradient/inter_g2_stats.csv'))
pd.DataFrame(np.array(inter_g3_stats)).to_csv(os.path.join(dkdir, 
            'gradient/inter_g3_stats.csv'))

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
                        fdr(np.array(inter_g3_stats)[:,1]))).T)\
                        .to_csv(os.path.join(dkdir, 
                                             'gradient/g_stats_fdr.csv'), 
                        header = ['intra_g1_fdr','intra_g2_fdr','intra_g3_fdr',
                                  'inter_g1_fdr','inter_g2_fdr','inter_g3_fdr'],
                                  index=None)
