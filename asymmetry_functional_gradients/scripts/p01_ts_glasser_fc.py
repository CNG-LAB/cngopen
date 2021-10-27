"""
calculates the rs-fMRI functional connectivity (FC) matrix for each HCP subject
and for each scan session, then gets the mean FC for each subject, 
saves the intra- and inter-hemispheric FC patterns separately 
"""
import numpy as np
from nilearn.connectome import ConnectivityMeasure
import nibabel as nib
import hcp_utils as hcp
import os
from numpy import inf
import pandas as pd

codedir = os.path.dirname(os.path.abspath(__file__))

datadir_day2 = os.path.join(os.path.dirname(codedir), 
                       'data/data_raw/hcp-functional-connectivity/')

datadir_day1 = os.path.join(os.path.dirname(codedir), 
                            'data/data_raw/HCP_S1200_rfMRI_DAY1/')

dataout = os.path.join(os.path.dirname(codedir), 
                       'data/data_HCP_glasser_fc/')

fcout = os.path.join(os.path.dirname(codedir), 
                       'data/data_results/FC/')

path_list = os.listdir(datadir_day2)
path_list.sort()

correlation_measure = ConnectivityMeasure(kind='correlation')

# read the time series of scan session day 2, LR, ** rfMRI_REST2_LR **,
# parcellate the time series to Glasser 360,
# compute the functional connectivity matrix and save

for dir in path_list:
    if os.path.exists(datadir_day2 + dir +
                      '/MNINonLinear/Results/rfMRI_REST2_LR/'\
                      'rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'):
        print('executing subject: '+ dir +'......')
        img = nib.load(datadir_day2 + dir+
                       '/MNINonLinear/Results/rfMRI_REST2_LR/'\
                       'rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii')
        data = img.get_fdata()
        data_parc = hcp.parcellate(data, hcp.mmp)
        data_cortex = data_parc[:,0:360]
        corr_matrix = correlation_measure.fit_transform([data_cortex])[0]
          
        fname_corr_matrix = dataout + 'DAY2_LR/corr/'+ dir +'.csv'
        np.savetxt(fname_corr_matrix, corr_matrix, delimiter = ',')
        print('finished subject: '+dir+'......')      
    else:
        fname_corr_matrix_failed = dataout + 'DAY2_LR/file_failed.txt'
        fail_file = open(fname_corr_matrix_failed, 'a')
        fail_file.write(dir+'\n')
        fail_file.close()
        print('subject '+ dir + ' file dose not exist.....')

# read the time series of scan session day 2, RL, ** rfMRI_REST2_RL **
# parcellate the time series to Glasser 360,
# compute the functional connectivity matrix and save

for dir in path_list:
    if os.path.exists(datadir_day2 + dir +
                      '/MNINonLinear/Results/rfMRI_REST2_RL/'\
                      'rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'):
        print('executing subject: '+ dir +'......')
        img = nib.load(datadir_day2 + dir +
                       '/MNINonLinear/Results/rfMRI_REST2_RL/'\
                       'rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii')
        data = img.get_fdata()
        data_parc = hcp.parcellate(data, hcp.mmp)
        data_cortex = data_parc[:,0:360]
        corr_matrix = correlation_measure.fit_transform([data_cortex])[0]
        fname_corr_matrix = dataout + 'DAY2_RL/corr/'+dir+'.csv'
        np.savetxt(fname_corr_matrix, corr_matrix, delimiter = ',')
        print('finished subject: '+dir+'......')
    else:
        fname_corr_matrix_failed = dataout + 'DAY2_RL/file_failed.txt'
        fail_file = open(fname_corr_matrix_failed, 'a')
        fail_file.write(dir+'\n')
        fail_file.close()
        print('subject '+ dir + ' file dose not exist.....')

# read the time series of scan session day 1, LR, ** rfMRI_REST1_LR **,
# parcellate the time series to Glasser 360,
# compute the functional connectivity matrix and save

for dir in path_list:
    if os.path.exists(datadir_day1 + dir + 
                      '/MNINonLinear/Results/rfMRI_REST1_LR/'\
                      'rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'):
        print('executing sub.'+ dir + '......')
        img = nib.load(datadir_day1 + dir +
                       '/MNINonLinear/Results/rfMRI_REST1_LR/'\
                       'rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii')
        data = img.get_fdata()
        data_parc = hcp.parcellate(data, hcp.mmp)
        data_cortex = data_parc[:,0:360]
        corr_matrix = correlation_measure.fit_transform([data_cortex])[0]
        fname_corr_matrix = dataout + 'DAY1_LR/corr/'+dir+'.csv'
        np.savetxt(fname_corr_matrix, corr_matrix, delimiter = ',')
        print('finish sub.'+dir+'......')
    else:
        fname_corr_matrix_failed = dataout + 'DAY1_LR/file_failed.txt'
        fail_file = open(fname_corr_matrix_failed, 'a')
        fail_file.write(dir+'\n')
        fail_file.close()
        print('sub.'+dir+'  file dose not exist.....')

# read the time series of scan session day 1, RL, ** rfMRI_REST1_RL **,
# parcellate the time series to Glasser 360,
# compute the functional connectivity matrix and save

for dir in path_list:
    if os.path.exists(datadir_day1 + dir + 
                      '/MNINonLinear/Results/rfMRI_REST1_RL/'\
                      'rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'):
        print('executing sub.'+dir+'......')
        img = nib.load(datadir_day1 + dir + 
                       '/MNINonLinear/Results/rfMRI_REST1_RL/'\
                       'rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii')
        data = img.get_fdata()
        data_parc = hcp.parcellate(data, hcp.mmp)
        data_cortex = data_parc[:,0:360]
        corr_matrix = correlation_measure.fit_transform([data_cortex])[0]
        fname_corr_matrix = dataout + 'DAY1_RL/corr/'+dir+'.csv'
        np.savetxt(fname_corr_matrix, corr_matrix, delimiter = ',')
        print('finish sub.'+dir+'......')
    else:
        fname_corr_matrix_failed = dataout + 'DAY1_RL/file_failed.txt'
        fail_file = open(fname_corr_matrix_failed, 'a')
        fail_file.write(dir+'\n')
        fail_file.close()
        print('sub.'+dir+'  file dose not exist.....')

# average functional connectivity matrices across four sessions,
# (across rfMRI_REST1_LR, rfMRI_REST1_RL, rfMRI_REST2_LR, rfMRI_REST2_RL)

for dir in path_list:
    if (os.path.exists(dataout + 'DAY1_LR/corr/' + dir + '.csv') and
        os.path.exists(dataout + 'DAY1_RL/corr/' + dir + '.csv') and
        os.path.exists(dataout + 'DAY2_LR/corr/' + dir + '.csv') and
        os.path.exists(dataout + 'DAY2_RL/corr/' + dir + '.csv')):
        
        print('executing sub.'+dir+'......')
        
        LR1 = np.array(pd.read_csv(dataout + 'DAY1_LR/corr/' + dir + '.csv',
                                 header=None))
        RL1 = np.array(pd.read_csv(dataout + 'DAY1_RL/corr/' + dir + '.csv',
                                 header=None))
        LR2 = np.array(pd.read_csv(dataout + 'DAY2_LR/corr/' + dir + '.csv',
                                 header=None))
        RL2 = np.array(pd.read_csv(dataout + 'DAY2_RL/corr/' + dir + '.csv',
                                 header=None))
        # average across four sessions
        mean = (LR1+RL1+LR2+RL2)/4
        # Fisher r to z transform
        fc = np.arctanh(mean) 
        fc[fc == inf] = 0
        np.savetxt(dataout + 'Mean/' + dir + '.csv', fc, delimiter = ',')

        # compute FC for left intra-hemisphere (LL), right intra- (RR),
        # for inter-hemispheres from left to right (LR) & right to left (RL)
        fc_LL = fc[0:180,0:180]
        fc_RR = fc[180:360,180:360]
        fc_LLRR = fc_LL - fc_RR
        fc_LR = fc[0:180,180:360]
        fc_RL = fc[180:360,0:180]
        fc_LRRL = fc_LR - fc_RL
        
        np.savetxt(fcout + 'LL/' + dir + '.csv', fc_LL, delimiter = ',')
        np.savetxt(fcout + 'RR/' + dir + '.csv', fc_RR, delimiter = ',')
        np.savetxt(fcout + 'LL-RR/' + dir + '.csv', fc_LLRR, delimiter = ',')
        np.savetxt(fcout + 'LR/' + dir + '.csv', fc_LR, delimiter = ',')
        np.savetxt(fcout + 'RL/' + dir + '.csv', fc_RL, delimiter = ',')
        np.savetxt(fcout + 'LR-RL/' + dir + '.csv', fc_LRRL, delimiter = ',')        
    else:
        unmatch = open(dataout + 'unmatch.txt', 'a')
        unmatch.write(dir+'\n')
        unmatch.close()
        print('sub.'+dir+'  file dose not exist.....')
