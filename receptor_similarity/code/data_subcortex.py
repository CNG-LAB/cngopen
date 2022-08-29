#!/usr/bin/env python3
"""
This script contains the subcortical data used in manuscript 
NEUROTRANSMITTER TRANSPORTER/RECEPTOR CO-EXPRESSION SHARES ORGANIZATIONAL TRAITS WITH BRAIN STRUCTURE AND FUNCTION
https://doi.org/10.1101/2022.08.26.505274
"""
from nilearn.input_data import NiftiMasker
from nilearn import image as nimg
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import zscore

input_path='path/to/data/'
mask_p='path/to/masks/'
#get subcortical data
def masking():
    path_out=input_path +'subcortex_masked/'
    images=os.listdir(input_path)
    mask = mask_p+'tian_binary_mask_total.nii.gz'
    masker=NiftiMasker(mask_img=mask).fit()
    for image in images:
        pet=nimg.load_img(input_path + image)
        masked=masker.transform(pet)
        name = re.search('.+?(?=\.)', image).group()
        np.savetxt(path_out + name + '_' + 'tian'+ '.csv', masked, delimiter=',')
    return

masking()

#make voxel-wise frame
path_in=input_path + 'subcortex_masked/'
all_files = os.listdir(path_in)
with open(input_path + 'receptor_list.txt', 'r') as f:
    s = f.read()
    s = s.split(',\n')
    studies = [re.search('.+?(?=\.)', study).group() for study in s]
receptors = {}
for study in studies:
    for i in all_files:
        if study in i and 'tian' in i:
            receptors[study] = i



rec_d = {}
for key, value in receptors.items():
    voxel_vals = np.genfromtxt(path_in + value, delimiter=',')
    # get average value for region
    rec_d[key] = voxel_vals
# returns the intensity per voxel

df = pd.DataFrame(rec_d)
df.columns = [x + '.csv' for x in df.columns]
_5HT1b = (zscore(df['5HT1b_p943_hc22_savli.csv']) * 22 +
          zscore(df['5HT1b_p943_hc65_gallezot.csv']) * 65) / (22 + 65)
_D2 = (zscore(df['D2_flb457_hc37_smith.csv']) * 37 +
       zscore(df['D2_flb457_hc55_sandiego.csv']) * 55) / (37 + 55)
_mGluR5 = (zscore(df['mGluR5_abp_hc22_rosaneto.csv']) * 22 +
           zscore(df['mGluR5_abp_hc28_dubois.csv']) * 28 +
           zscore(df['mGluR5_abp_hc73_smart.csv']) * 73) / (22 + 28 + 73)
_VAChT = (zscore(df['VAChT_feobv_hc18_aghourian_sum.csv']) * 18 +
          zscore(df['VAChT_feobv_hc4_tuominen.csv']) * 4 +
          zscore(df['VAChT_feobv_hc5_bedard_sum.csv']) * 5) / (18 + 4 + 5)

d = {'5HT1a': df['5HT1a_way_hc36_savli.csv'], '5HT1b': _5HT1b, '5HT2a': df['5HT2a_cimbi_hc29_beliveau.csv'],
     '5HT4': df['5HT4_sb20_hc59_beliveau.csv'], '5HT6': df['5HT6_gsk_hc30_radnakrishnan.csv'],
     '5HTT': df['5HTT_dasb_hc100_beliveau.csv'],
     'A4B2': df['A4B2_flubatine_hc30_hillmer.csv'], 'CB1': df['CB1_omar_hc77_normandin.csv'],
     'D1': df['D1_SCH23390_hc13_kaller.csv'], 'D2': _D2, 'DAT': df['DAT_fpcit_hc174_dukart_spect.csv'],
     'GABAa': df['GABAa-bz_flumazenil_hc16_norgaard.csv'], 'H3': df['H3_cban_hc8_gallezot.csv'],
     'M1': df['M1_lsn_hc24_naganawa.csv'], 'mGluR5': _mGluR5, 'MU': df['MU_carfentanil_hc204_kantonen.csv'],
     'NAT': df['NAT_MRB_hc77_ding.csv'], 'NMDA': df['NMDA_ge179_hc29_galovic.csv'], 'VAChT': _VAChT}

df_wm_reg = pd.DataFrame(d)
df_wm_reg.to_csv(path_in +'subcortex_tian.csv')

#make regional frame
whole_mask=mask_p+ 'tian_binary_mask_total.nii.gz'
whole_mask=NiftiMasker(mask_img=whole_mask).fit()
in_file=df_wm_reg
in_file=in_file.apply(zscore)
path_out=input_path + 'subcortex_tian_regions.csv'
all_files=os.listdir(mask_path)
del all_files[-1] #delete total mask
del all_files[-7] #delete amygdala lh
del all_files[-7] #delete amygdala rh
all_files.sort()
results={}
for subregion in range(len(all_files)):
    mask=NiftiMasker(mask_path + all_files[subregion]).fit()
    for receptor in in_file:
        demasked=whole_mask.inverse_transform(in_file[receptor])
        raw=mask.transform(demasked)
        avg=np.mean(raw)
        if receptor in results.keys():
            results[receptor].append(avg)
        else:
            results[receptor]=[avg]

df=pd.DataFrame(results)
df.to_csv(path_in +'subcortex_tian_regions.csv')
