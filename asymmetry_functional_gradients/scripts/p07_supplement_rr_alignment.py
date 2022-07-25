'''
To address the potential issue that alignment to different hemispheres would cause bias results, here we
aligned the individual gradeints to the RR template. Following we correlate the results of 
LL with RR to see whether the patterns are similar (or not).
'''
import numpy as np
import nibabel as nib
import os
import shutil
import pandas as pd
from brainspace.gradient import GradientMaps
from scipy import stats

path = '../data/data_results/FC/'

try:
    os.remove(path+'LL/'+'.DS_Store')
    os.remove(path+'RR/'+'.DS_Store')
    os.remove(path+'LR/'+'.DS_Store')
    os.remove(path+'RL/'+'.DS_Store')
except:    
    pass

try:
  shutil.rmtree('../data/data_results/supplementary/rr_align/grad/', ignore_errors=False, onerror=None)
except:
  pass

try:
  os.mkdir('../data/data_results/supplementary/rr_align/grad/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/LL/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/RR/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/LR/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/RL/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/intra/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/inter/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/network/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/network/intra/')
  os.mkdir('../data/data_results/supplementary/rr_align/grad/network/inter/')
except:
  pass

path_list = os.listdir(path+'LL/')
path_list.sort()

# group gradients
RR = np.array(pd.read_csv('../data/data_results/FC/RR_groupmean.csv', header=None))
gm = GradientMaps(approach='dm', kernel='normalized_angle',n_components=10,random_state=0)
gm.fit(RR)
group_grad = gm.gradients_

path_add = '../data/data_results/supplementary/rr_align/grad/'

np.savetxt(path_add+'group_grad_RR.csv', group_grad)
np.savetxt(path_add+'group_grad_RR_lambdas.csv', gm.lambdas_)

# individual gradients
for i in path_list:
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')  
  FC_LL = np.array(pd.read_csv('../data/data_results/FC/LL/'+i, header=None))
  align.fit(FC_LL,reference=group_grad)
  grad_LL = align.aligned_
  np.savetxt(path_add+'LL/'+i, grad_LL)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  FC_RR = np.array(pd.read_csv('../data/data_results/FC/RR/'+i, header=None))
  align.fit(FC_RR,reference=group_grad)
  grad_RR = align.aligned_
  np.savetxt(path_add+'RR/'+i, grad_RR)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  FC_LR = np.array(pd.read_csv('../data/data_results/FC/LR/'+i, header=None))
  align.fit(FC_LR,reference=group_grad)
  grad_LR = align.aligned_
  np.savetxt(path_add+'LR/'+i, grad_LR)
  align = GradientMaps(n_components=10, random_state=0, approach='dm', 
                       kernel='normalized_angle', alignment='procrustes')
  FC_RL = np.array(pd.read_csv('../data/data_results/FC/RL/'+i, header=None))
  align.fit(FC_RL,reference=group_grad)
  grad_RL = align.aligned_
  np.savetxt(path_add+'RL/'+i, grad_RL)
  print('finish   ' + i)

total_llrr = 0
total_lrrl = 0
# correct individual gradient if the correlation is negative
for dir in path_list:
  #LL
  df = np.loadtxt(path_add+'LL/'+dir)
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(group_grad[:,i],df[:,i])    
    if r[i][0] > 0:
      corrected[i]=df[:,i]
    else:
      corrected[i]=-1*df[:,i]
  corrected_ll = np.array(corrected).T
  np.savetxt(path_add+'LL/'+dir, corrected_ll)
  
  # RR
  df = np.loadtxt(path_add+'RR/'+dir)
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(group_grad[:,i],df[:,i])    
    if r[i][0] > 0:
      corrected[i]=df[:,i]
    else:
      corrected[i]=-1*df[:,i]
  corrected_rr = np.array(corrected).T
  np.savetxt(path_add+'RR/'+dir, corrected_rr)
    
  # LR
  df = np.loadtxt(path_add+'LR/'+dir)
  r = [None] * 10   
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(group_grad[:,i],df[:,i])    
    if r[i][0] > 0:
      corrected[i]=df[:,i]
    else:
      corrected[i]=-1*df[:,i]
  corrected_lr = np.array(corrected).T
  np.savetxt(path_add+'LR/'+dir, corrected_lr)
    
  # RL
  df = np.loadtxt(path_add+'RL/'+dir)
  r = [None] * 10
  corrected = [None]*10
  for i in range(10):
    r[i] = stats.pearsonr(group_grad[:,i],df[:,i])    
    if r[i][0] > 0:
      corrected[i]=df[:,i]
    else:
      corrected[i]=-1*df[:,i]
  corrected_rl = np.array(corrected).T
  np.savetxt(path_add+'RL/'+dir, corrected_rl)

  # RR-LL, RL-RL
  AI_llrr = corrected_ll - corrected_rr
  AI_lrrl = corrected_lr - corrected_rl
  np.savetxt(path_add+'intra/'+dir, AI_llrr)
  np.savetxt(path_add+'inter/'+dir, AI_lrrl)
  print('finish   ' + dir)
  total_llrr = total_llrr + AI_llrr
  total_lrrl = total_lrrl + AI_lrrl
  
mean_llrr = total_llrr/len(path_list)
mean_lrrl = total_lrrl/len(path_list)

np.savetxt(path_add+'mean_asym_LLRR.csv', mean_llrr)
np.savetxt(path_add+'mean_asym_LRRL.csv', mean_lrrl)

# ca network
ca_l = np.array(pd.read_csv('../data/ca_glasser_network.csv',header=None))[:,0][:180]
ca_r = np.array(pd.read_csv('../data/ca_glasser_network.csv',header=None))[:,0][180:]

for n in range(len(path_list)):
  ll = np.loadtxt(path_add+'LL/'+path_list[n])
  rr = np.loadtxt(path_add+'RR/'+path_list[n])
  lr = np.loadtxt(path_add+'LR/'+path_list[n])
  rl = np.loadtxt(path_add+'RL/'+path_list[n])
  intra = [None] * 3  
  for i in range(3):
    intra[i] = [np.mean(ll[:,i][ca_l==1])-np.mean(rr[:,i][ca_r==1]),
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
                np.mean(ll[:,i][ca_l==12])-np.mean(rr[:,i][ca_r==12])]
  np.savetxt(path_add+'network/intra/'+path_list[n], np.array(intra).T)

  inter = [None] * 3  
  for i in range(3):
    inter[i] = [np.mean(lr[:,i][ca_l==1])-np.mean(rl[:,i][ca_r==1]),
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
                np.mean(lr[:,i][ca_l==12])-np.mean(rl[:,i][ca_r==12])]
  
  np.savetxt(path_add+'network/inter/'+path_list[n], np.array(inter).T)
  print('finish......'+path_list[n])


intra_g1 = [None] * len(path_list) # diff between LL and RR gradient 1
intra_g2 = [None] * len(path_list)
intra_g3 = [None] * len(path_list)
inter_g1 = [None] * len(path_list) # diff between LR and RL gradient 1
inter_g2 = [None] * len(path_list)
inter_g3 = [None] * len(path_list)    


for n in range(len(path_list)):
  intra_g1[n] = np.loadtxt(path_add+'intra/'+path_list[n])[:,0]
  intra_g2[n] = np.loadtxt(path_add+'intra/'+path_list[n])[:,1]
  intra_g3[n] = np.loadtxt(path_add+'intra/'+path_list[n])[:,2]
  inter_g1[n] = np.loadtxt(path_add+'inter/'+path_list[n])[:,0]
  inter_g2[n] = np.loadtxt(path_add+'inter/'+path_list[n])[:,1]
  inter_g3[n] = np.loadtxt(path_add+'inter/'+path_list[n])[:,2]

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

pd.DataFrame(np.array(intra_g1_stats)).to_csv(os.path.join(path_add, 'intra_g1_stats.csv' ))
pd.DataFrame(np.array(intra_g2_stats)).to_csv(os.path.join(path_add, 'intra_g2_stats.csv' ))
pd.DataFrame(np.array(intra_g3_stats)).to_csv(os.path.join(path_add, 'intra_g3_stats.csv' ))
pd.DataFrame(np.array(inter_g1_stats)).to_csv(os.path.join(path_add, 'inter_g1_stats.csv' ))
pd.DataFrame(np.array(inter_g2_stats)).to_csv(os.path.join(path_add, 'inter_g2_stats.csv' ))
pd.DataFrame(np.array(inter_g3_stats)).to_csv(os.path.join(path_add, 'inter_g3_stats.csv' ))

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
                        to_csv(os.path.join(path_add, 'g_stats_fdr.csv'),  
                               header = ['intra_g1_fdr','intra_g2_fdr',
                                         'intra_g3_fdr','inter_g1_fdr',
                                         'inter_g2_fdr','inter_g3_fdr'],
                                         index=None)
