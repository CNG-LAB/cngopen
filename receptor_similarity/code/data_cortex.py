#!/usr/bin/env python3
"""
This script generates the data used in manuscript XXX
"""
from netneurotools.networks import struct_consensus
from netneurotools.freesurfer import find_parcel_centroids
from nilearn.input_data import NiftiLabelsMasker
from nilearn import datasets
from nilearn._utils import check_niimg
import os
import pandas as pd
from scipy.stats import zscore
import re
import numpy as np

parcels = 100
input_path = 'path/to/data/'

#pet data
#parcellate study maps
atl='schaefer'
parc = 100
scale = '{}Parcels7Networks'.format(parc)

dataset=datasets.fetch_atlas_schaefer_2018(n_rois=parc, yeo_networks=7)
outpath=input_path + "Parcellated\{}\{}".format(atl, scale)
if not os.path.exists(outpath):
    os.makedirs(outpath)

receptors_nii = [input_path+'5HT1a_way_hc36_savli.nii',
                 input_path+'5HT1a_cumi_hc8_beliveau.nii',
                 input_path+'5HT1b_az_hc36_beliveau.nii',
                 input_path+'5HT1b_p943_hc22_savli.nii',
                 input_path+'5HT1b_p943_hc65_gallezot.nii.gz',
                 input_path+'5HT2a_cimbi_hc29_beliveau.nii',
                 input_path+'5HT2a_alt_hc19_savli.nii',
                 input_path+'5HT2a_mdl_hc3_talbot.nii.gz',
                 input_path+'5HT4_sb20_hc59_beliveau.nii',
                 input_path+'5HT6_gsk_hc30_radnakrishnan.nii.gz',
                 input_path+'5HTT_dasb_hc100_beliveau.nii',
                 input_path+'5HTT_dasb_hc30_savli.nii',
                 input_path+'A4B2_flubatine_hc30_hillmer.nii.gz',
                 input_path+'CB1_omar_hc77_normandin.nii.gz',
                 input_path+'CB1_FMPEPd2_hc22_laurikainen.nii',
                 input_path+'D1_SCH23390_hc13_kaller.nii',
                 input_path+'D2_fallypride_hc49_jaworska.nii',
                 input_path+'D2_flb457_hc37_smith.nii.gz',
                 input_path+'D2_flb457_hc55_sandiego.nii.gz',
                 input_path+'DAT_fpcit_hc174_dukart_spect.nii',
                 input_path+'DAT_fepe2i_hc6_sasaki.nii.gz',
                 input_path+'GABAa-bz_flumazenil_hc16_norgaard.nii',
                 input_path+'GABAa_flumazenil_hc6_dukart.nii',
                 input_path+'H3_cban_hc8_gallezot.nii.gz',
                 input_path+'M1_lsn_hc24_naganawa.nii.gz',
                 input_path+'mGluR5_abp_hc22_rosaneto.nii',
                 input_path+'mGluR5_abp_hc28_dubois.nii',
                 input_path+'mGluR5_abp_hc73_smart.nii',
                 input_path+'MU_carfentanil_hc204_kantonen.nii',
                 input_path+'MU_carfentanil_hc39_turtonen.nii',
                 input_path+'NAT_MRB_hc77_ding.nii.gz',
                 input_path+'NAT_MRB_hc10_hesse.nii',
                 input_path+'NMDA_ge179_hc29_galovic.nii.gz',
                 input_path+'VAChT_feobv_hc4_tuominen.nii',
                 input_path+'VAChT_feobv_hc5_bedard_sum.nii',
                 input_path+'VAChT_feobv_hc18_aghourian_sum.nii']

parcellated = {}
mask = NiftiLabelsMasker(dataset['maps'], resampling_target='data', strategy='mean').fit()

for receptor in receptors_nii:
    img = check_niimg(receptor, atleast_4d=True)
    parcellated[receptor] = mask.transform(img).squeeze()
    name = receptor.split('\\')[-1]  # get nifti file name
    name = name.split('.')[0]  # remove .nii
    np.savetxt(outpath+'\\'+ name+'.csv', parcellated[receptor], delimiter=',')

#generate comprehensive df
all_files = os.listdir(outpath)
with open(input_path + 'NTRM_interest.txt', 'r') as f:
    s = f.read()
    s = s.split(',\n')
    studies = [re.search('.+?(?=\.)', study).group() for study in s]
receptors = {}
for study in studies:
    for i in all_files:
        if study in i:
            receptors[study] = i

# make dataframe from dict and returns output
vals={}
for key, val in receptors.items():
    vals[key]=np.genfromtxt(outpath + val, delimiter=',')
df=pd.DataFrame(vals)
df.columns = [x + '.csv' for x in df.columns]
# generate weighted averages
_5HT1b = (zscore(df['5HT1b_p943_hc22_savli.csv']) * 22 +
          zscore(df['5HT1b_p943_hc65_gallezot.csv']) * 65) / (22 + 65)
_D2 = (zscore(df['D2_flb457_hc37_smith.csv']) * 37 +
       zscore(df['D2_flb457_hc55_sandiego.csv']) * 55) / (37 + 55)
_mGluR5 = (zscore(df['mGluR5_abp_hc22_rosaneto.csv']) * 22 +
           zscore(df['mGluR5_abp_hc28_dubois.csv']) * 28 +
           zscore(df['mGluR5_abp_hc73_smart.csv']) * 73) / (22 + 28 + 73)
_VAChT = (zscore(df['VAChT_feobv_hc18_aghourian_sum.csv']) * 18 +
          zscore(df['VAChT_feobv_hc4_tuominen.csv']) * 4 +
          zscore(df['VAChT_feobv_hc5_bedard_sum.csv'])* 5) / (18 + 4 + 5)

d = {'5HT1a': df['5HT1a_way_hc36_savli.csv'], '5HT1b': _5HT1b, '5HT2a': df['5HT2a_cimbi_hc29_beliveau.csv'],
     '5HT4': df['5HT4_sb20_hc59_beliveau.csv'], '5HT6': df['5HT6_gsk_hc30_radnakrishnan.csv'],
     '5HTT': df['5HTT_dasb_hc100_beliveau.csv'],
     'A4B2': df['A4B2_flubatine_hc30_hillmer.csv'], 'CB1': df['CB1_omar_hc77_normandin.csv'],
     'D1': df['D1_SCH23390_hc13_kaller.csv'], 'D2': _D2, 'DAT': df['DAT_fpcit_hc174_dukart_spect.csv'],
     'GABAa': df['GABAa-bz_flumazenil_hc16_norgaard.csv'], 'H3': df['H3_cban_hc8_gallezot.csv'],
     'M1': df['M1_lsn_hc24_naganawa.csv'], 'mGluR5': _mGluR5, 'MU': df['MU_carfentanil_hc204_kantonen.csv'],
     'NAT': df['NAT_MRB_hc77_ding.csv'],'NMDA' : df['NMDA_ge179_hc29_galovic.csv'], 'VAChT': _VAChT}

df_wm_reg = pd.DataFrame(d)
df_wm_reg.to_csv(input_path +'{}_receptorprofiles.csv'.format(scale))

#SC, FC, MPC matrices
path_mics='path/to/mics/'
subjects=os.listdir(path)
subjects.pop(0)

sc_f=path_mics+'{}\ses-01\dwi\{}_ses-01_space-dwinative_atlas-schaefer{}_desc-sc.txt'
fc_f=path_mics+'{}\ses-01\func\{}_ses-01_space-fsnative_atlas-schaefer{}_desc-fc.txt'
mpc_f=path_mics+'{}\ses-01\anat\micro_profiles\{}_ses-01_space-fsnative_atlas-schaefer{}_desc-mpc.txt'

#fc
fc=np.zeros((parc,parc))
for i, subject in enumerate(subjects):
    s_sc=np.genfromtxt(fc_f.format(subject, subject, parc), delimiter=',')
    to_drop=list(range(15))
    to_drop.append(int(parc / 2 ) + 15)
    s_sc=np.delete(s_sc, to_drop, axis=0)
    s_sc=np.delete(s_sc, to_drop, axis=1)
    FCz = np.arctanh(s_sc)

    # replace inf with 0
    FCz[~np.isfinite(FCz)] = 0

    # Mirror the matrix
    FCz = np.triu(FCz,1)+FCz.T
    fc += FCz
fc = fc / 50

#mpc
mpc=np.zeros((parc,parc))
for file in subjects:
    s_mpc=np.genfromtxt(mpc_f.format(file, file, parc))
    to_drop=[0,int(parc / 2) + 1]
    s_mpc=np.delete(s_mpc, to_drop, axis=0)
    s_mpc=np.delete(s_mpc, to_drop, axis=1)
    s_mpc=np.triu(s_mpc,1)+s_mpc.T
    mpc+= s_mpc
mpc = mpc / 50

#sc
scs=[]
for i, subject in enumerate(subjects):
    s_sc=np.genfromtxt(sc_f.format(subject, subject, parc), delimiter=',')
    to_drop=list(range(15))
    to_drop.append(int(parc / 2 ) + 15)
    s_sc=np.delete(s_sc, to_drop, axis=0)
    s_sc=np.delete(s_sc, to_drop, axis=1)
    s_sc = np.log(np.triu(s_sc,1)+s_sc.T)
    s_sc[np.isneginf(s_sc)] = 0
    scs.append(s_sc)
sc=np.dstack(scs)

lh=input_path + 'lh.Schaefer2018_{}Parcels_7Networks_order.annot'.format(parc)
rh=input_path + 'rh.Schaefer2018_{}Parcels_7Networks_order.annot'.format(parc)

centroids=find_parcel_centroids(lhannot=lh, rhannot=rh, version='fsaverage5')
points=centroids[0]
dist=np.empty((parc,parc))
for i in range(len(points)):
    for j in range(len(points)):
        dist[i,j]=np.linalg.norm(points[i]-points[j])


hemiid=np.concatenate((np.ones(int(parc / 2)),np.zeros(int(parc / 2))))
hemiid=hemiid.reshape(-1,1)

SC=struct_consensus(sc, dist, hemiid, True)


np.save(input_path + 'mics_sc_{}.npy'.format(parc), SC)
np.save(input_path + 'mics_fc_{}.npy'.format(parc), fc)
np.save(input_path + 'mics_mpc_{}.npy'.format(parc), mpc)





