"""
computes the gradients of functional connectivity for 4 different fashions:
LL & RR intra-hemispheric and LR & RL inter-hemispheric fashion for macaques.
we implement one sample t-test on intra-hemispheric and 
inter-hemispheric gradient score differences in each parcel (LL-RR differences
and LR-RL differences) and finally correct p-values using FDR.
"""
import os
import numpy as np
import pandas as pd
from numpy import inf
import scipy.io

codedir = os.path.dirname(os.path.abspath(__file__))

macdir = os.path.join(os.path.dirname(codedir), 
                      'data/macaque_data')

outdir = os.path.join(os.path.dirname(codedir), 
                      'data/data_results/macaque')

# load macaque functional connectivity data (n=19 macaque, 182 parcels)
data = scipy.io.loadmat(os.path.join(macdir, '19_macaques.mat'))
corr = data['mcc_fcconn'] # correlation, size (19, 182, 182)

# Fisher r to z transform for each macaque fc
macaque = [None]*19
for i in range(19):
  macaque[i] = np.arctanh(corr[i]) 
  macaque[i][macaque[i] == inf] = 0

# get fc for 4 fashions: intra- (LL, RR) and inter-hemispheric (LR, RL)
mq_ll = [None] * 19
mq_rr = [None] * 19
mq_lr = [None] * 19
mq_rl = [None] * 19
ll_total = 0
for i in range (19):
    mq_ll[i] = macaque[i][0:91,0:91]
    mq_rr[i] = macaque[i][91:182,91:182]
    mq_lr[i] = macaque[i][0:91,91:182]    
    mq_rl[i] = macaque[i][91:182,0:91]
    ll_total = ll_total + mq_ll[i]

# get the mean fc for LL
ll_mean = ll_total/19

# get the gradients of mean fc LL -> this will be the reference grads. set
from brainspace.gradient import GradientMaps
gm= GradientMaps(approach = 'dm', kernel='normalized_angle',
                 n_components=10, random_state=0)
gm.fit(ll_mean)
grad_ref = gm.gradients_
np.savetxt(os.path.join(outdir, 'group_grad_LL.csv'), grad_ref, 
           delimiter = ',')
np.savetxt(os.path.join(outdir,'group_grad_LL_lambdas.csv'), gm.lambdas_, 
           delimiter = ',')

# get the individual gradients by aligning them to the reference grads.
ll_grad = [None] * 19
rr_grad = [None] * 19
lr_grad = [None] * 19
rl_grad = [None] * 19

for i in range(19):
    # ll
    gm_ll = GradientMaps(approach = 'dm', kernel='normalized_angle', 
                        alignment='procrustes',n_components=10, random_state=0)
    gm_ll.fit(mq_ll[i],reference = grad_ref)
    ll_grad[i] = gm_ll.aligned_
    # rr
    gm_rr = GradientMaps(approach = 'dm', kernel='normalized_angle', 
                        alignment='procrustes',n_components=10, random_state=0)
    gm_rr.fit(mq_rr[i],reference = grad_ref)
    rr_grad[i] = gm_rr.aligned_
    # lr
    gm_lr = GradientMaps(approach = 'dm', kernel='normalized_angle', 
                        alignment='procrustes',n_components=10, random_state=0)
    gm_lr.fit(mq_lr[i],reference = grad_ref)
    lr_grad[i] = gm_lr.aligned_
    # rl
    gm_rl = GradientMaps(approach = 'dm', kernel='normalized_angle', 
                        alignment='procrustes',n_components=10, random_state=0)
    gm_rl.fit(mq_rl[i],reference = grad_ref)
    rl_grad[i] = gm_rl.aligned_

# quality check: get the correlations between individual gradients and 
# the group-level template gradient, swap the axis if negative    
# repeat it all for LL, RR, LR and RL individual gradients
from scipy import stats
corrected_ll = [None] * 19
for n in range(19):
  a = ll_grad[n]
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(grad_ref[:,i],a[:,i])    
    if r[i][0] > 0:
      corrected[i]=a[:,i]
    else:
      corrected[i]=-1*a[:,i]
  corrected_ll[n] = np.concatenate((corrected[0].reshape(corrected[0].shape[0],1),
                                    corrected[1].reshape(corrected[1].shape[0],1),
                                    corrected[2].reshape(corrected[2].shape[0],1),
                                    corrected[3].reshape(corrected[3].shape[0],1),
                                    corrected[4].reshape(corrected[4].shape[0],1),
                                    corrected[5].reshape(corrected[5].shape[0],1),
                                    corrected[6].reshape(corrected[6].shape[0],1),
                                    corrected[7].reshape(corrected[7].shape[0],1),
                                    corrected[8].reshape(corrected[8].shape[0],1),
                                    corrected[9].reshape(corrected[9].shape[0],1)),
                                   axis = 1)
  np.savetxt(os.path.join(outdir, 'LL', str(n)+'.csv'), 
                          corrected_ll[n], delimiter = ',')

corrected_rr = [None] * 19
for n in range(19):
  a = rr_grad[n]
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(grad_ref[:,i],a[:,i])    
    if r[i][0] > 0:
      corrected[i]=a[:,i]
    else:
      corrected[i]=-1*a[:,i]
  corrected_rr[n] = np.concatenate((corrected[0].reshape(corrected[0].shape[0],1),
                                    corrected[1].reshape(corrected[1].shape[0],1),
                                    corrected[2].reshape(corrected[2].shape[0],1),
                                    corrected[3].reshape(corrected[3].shape[0],1),
                                    corrected[4].reshape(corrected[4].shape[0],1),
                                    corrected[5].reshape(corrected[5].shape[0],1),
                                    corrected[6].reshape(corrected[6].shape[0],1),
                                    corrected[7].reshape(corrected[7].shape[0],1),
                                    corrected[8].reshape(corrected[8].shape[0],1),
                                    corrected[9].reshape(corrected[9].shape[0],1)),
                                    axis = 1)
  np.savetxt(os.path.join(outdir, 'RR', str(n)+'.csv'),
                          corrected_rr[n], delimiter = ',')

corrected_lr = [None] * 19
for n in range(19):
  a = lr_grad[n]
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(grad_ref[:,i],a[:,i])    
    if r[i][0] > 0:
      corrected[i]=a[:,i]
    else:
      corrected[i]=-1*a[:,i]
  corrected_lr[n] = np.concatenate((corrected[0].reshape(corrected[0].shape[0],1),
                                    corrected[1].reshape(corrected[1].shape[0],1),
                                    corrected[2].reshape(corrected[2].shape[0],1),
                                    corrected[3].reshape(corrected[3].shape[0],1),
                                    corrected[4].reshape(corrected[4].shape[0],1),
                                    corrected[5].reshape(corrected[5].shape[0],1),
                                    corrected[6].reshape(corrected[6].shape[0],1),
                                    corrected[7].reshape(corrected[7].shape[0],1),
                                    corrected[8].reshape(corrected[8].shape[0],1),
                                    corrected[9].reshape(corrected[9].shape[0],1)),axis = 1)
  np.savetxt(os.path.join(outdir, 'LR', str(n)+'.csv'),
                          corrected_lr[n], delimiter = ',')

corrected_rl = [None] * 19
for n in range(19):
  a = rl_grad[n]
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(grad_ref[:,i],a[:,i])    
    if r[i][0] > 0:
      corrected[i]=a[:,i]
    else:
      corrected[i]=-1*a[:,i]
  corrected_rl[n] = np.concatenate((corrected[0].reshape(corrected[0].shape[0],1),
                                    corrected[1].reshape(corrected[1].shape[0],1),
                                    corrected[2].reshape(corrected[2].shape[0],1),
                                    corrected[3].reshape(corrected[3].shape[0],1),
                                    corrected[4].reshape(corrected[4].shape[0],1),
                                    corrected[5].reshape(corrected[5].shape[0],1),
                                    corrected[6].reshape(corrected[6].shape[0],1),
                                    corrected[7].reshape(corrected[7].shape[0],1),
                                    corrected[8].reshape(corrected[8].shape[0],1),
                                    corrected[9].reshape(corrected[9].shape[0],1)),axis = 1)
  np.savetxt(os.path.join(outdir, 'RL', str(n)+'.csv'), 
                          corrected_rl[n], delimiter = ',')

# compute the Asymmetry Index (AI) for each macaque: 
# intra-hemispheric AI = LL (gradients) - RR (gradients)
# inter-hemispheric AI = LR (gradients) - RL (gradients)
AI_intra = [None] * 19
AI_inter = [None] * 19
for i in range(19):
  AI_intra[i] = corrected_ll[i] - corrected_rr[i]
  AI_inter[i] = corrected_lr[i] - corrected_rl[i]

# get the mean AI score for each parcel along individual gradients (intra-)
import statistics
AI_intra_mean = np.array([[None] * 91] * 10)
for n in range(10):
  for i in range(91):
      AI_intra_mean[:,i][n] = statistics.mean(np.array(AI_intra).T[:,i][n])
        
pd.DataFrame(AI_intra_mean).T.to_csv(os.path.join(outdir, 
                                     'macaque_91_gradient_intra_mean.csv'),
                                     header=None, index=None)

# get the mean AI score for each parcel along individual gradients (inter-)
AI_inter_mean = np.array([[None] * 91] * 10)
for n in range(10):
  for i in range(91):
      AI_inter_mean[:,i][n] = statistics.mean(np.array(AI_inter).T[:,i][n])

pd.DataFrame(AI_inter_mean).T.to_csv(os.path.join(outdir, 
                                     'macaque_91_gradient_inter_mean.csv'),
                                     header=None, index=None)

path = os.path.join(outdir, 'LL')
path_list = os.listdir(path)
path_list.sort()

for i in range(19):
    intra = np.array(pd.read_csv(os.path.join(outdir, 'LL', path_list[i]), header = None)) \
            - np.array(pd.read_csv(os.path.join(outdir, 'RR', path_list[i]), header = None))
    np.savetxt(os.path.join(outdir, 'intra', path_list[i]), 
               intra, delimiter = ',')
    
    inter = np.array(pd.read_csv(os.path.join(outdir, 'LR', path_list[i]), header = None)) \
            - np.array(pd.read_csv(os.path.join(outdir, 'RL', path_list[i]), header = None))
    np.savetxt(os.path.join(outdir, 'inter', path_list[i]), 
               inter, delimiter = ',')


# Stats: implement one-sample t-test on intra-hemispheric (LL, RR) and 
# inter-hemispheric (LR, RL) AI-scores,
# for each parcel (1,...,91) and for each gradient (1,2,3) separately
from scipy import stats
AI_intra_g1_stats = [None] * 91
AI_intra_g2_stats = [None] * 91
AI_intra_g3_stats = [None] * 91
AI_inter_g1_stats = [None] * 91
AI_inter_g2_stats = [None] * 91
AI_inter_g3_stats = [None] * 91

for i in range(91):
  AI_intra_g1_stats[i] = stats.ttest_1samp(np.array(AI_intra)[0:19,:,0][:,i],0)
  AI_intra_g2_stats[i] = stats.ttest_1samp(np.array(AI_intra)[0:19,:,1][:,i],0)
  AI_intra_g3_stats[i] = stats.ttest_1samp(np.array(AI_intra)[0:19,:,2][:,i],0)
  AI_inter_g1_stats[i] = stats.ttest_1samp(np.array(AI_inter)[0:19,:,0][:,i],0)
  AI_inter_g2_stats[i] = stats.ttest_1samp(np.array(AI_inter)[0:19,:,1][:,i],0)
  AI_inter_g3_stats[i] = stats.ttest_1samp(np.array(AI_inter)[0:19,:,2][:,i],0)

a = np.vstack((np.array(AI_intra_g1_stats).T,
               np.array(AI_intra_g2_stats).T,
               np.array(AI_intra_g3_stats).T,
               np.array(AI_inter_g1_stats).T,
               np.array(AI_inter_g2_stats).T,
               np.array(AI_inter_g3_stats).T))

pd.DataFrame(a.T, columns = ['intra_g1_t','intra_g1_p','intra_g2_t',
                             'intra_g2_p','intra_g3_t','intra_g3_p',
                             'inter_g1_t','inter_g1_p','inter_g2_t',
                             'inter_g2_p','inter_g3_t','inter_g3_p'],
                              index = None).to_csv(os.path.join(outdir, 
                         'macaque_asymmetric_gradients_stats.csv'), index=None)

# FDR correction on the p-values
def fdr(p_vals):
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

b = np.vstack((fdr(np.array(AI_intra_g1_stats)[:,1]).T,
               fdr(np.array(AI_intra_g2_stats)[:,1]).T,
               fdr(np.array(AI_intra_g3_stats)[:,1]).T,
               fdr(np.array(AI_inter_g1_stats)[:,1]).T,
               fdr(np.array(AI_inter_g2_stats)[:,1]).T,
               fdr(np.array(AI_inter_g3_stats)[:,1]).T))

pd.DataFrame(b.T, columns = ['intra_g1_fdr','intra_g2_fdr','intra_g3_fdr',
                             'inter_g1_fdr','inter_g2_fdr','inter_g3_fdr'],
             index = None).to_csv(os.path.join(outdir, 
                         'macaque_asymmetric_gradients_fdr.csv'), index=None)


