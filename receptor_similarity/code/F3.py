#!/usr/bin/env python3
"""
This script contains the code behind the results in F3 in manuscript 
Cerebral chemoarchitecture shares organizational traits with brain structure and function
Elife. 2023 Jul 13;12:e83843. doi: 10.7554/eLife.83843.
"""
from enigmatoolbox import datasets
from enigmatoolbox.utils.parcellation import parcel_to_surface
import nilearn
from brainspace.null_models import SpinPermutations
import matplotlib.pyplot as plt
import plotly.express as px
from mn_funcs import spin, schaefer_to_surf
from plotly.graph_objs import *
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import nibabel as nib
from nilearn.input_data import NiftiLabelsMasker
import re
import os

parcels = 100
input_path = 'path/to/data/'
res_p = 'path/to/results/F3/'

rc_g1 = np.load(input_path + 'rc_g1_{}.npy')
rc_g2 = np.load(input_path + 'rc_g2_{}.npy')
rc_g3 = np.load(input_path + 'rc_g3_{}.npy')

# cortical thinning
disorders = ['22q', 'adhd', 'asd', 'bipolar', 'depression', 'epilepsy', 'ocd', 'schizophrenia']
d = {}

# 22q
sum_stats = datasets.load_summary_stats('22q')
CT = sum_stats['CortThick_case_vs_controls']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['22q'] = ct_surf

# adhd
sum_stats = datasets.load_summary_stats('adhd')
CT = sum_stats['CortThick_case_vs_controls_adult']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['adhd'] = ct_surf

# asd
sum_stats = datasets.load_summary_stats('asd')
CT = sum_stats['CortThick_case_vs_controls_meta_analysis']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['asd'] = ct_surf

# bipolar
sum_stats = datasets.load_summary_stats('bipolar')
CT = sum_stats['CortThick_case_vs_controls_adult']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['bipolar'] = ct_surf

# epilepsy
sum_stats = datasets.load_summary_stats('epilepsy')
CT = sum_stats['CortThick_case_vs_controls_allepilepsy']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['allepilepsy'] = ct_surf

# depression
sum_stats = datasets.load_summary_stats('depression')
CT = sum_stats['CortThick_case_vs_controls_adult']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['depression_all'] = ct_surf

# ocd
sum_stats = datasets.load_summary_stats('ocd')
CT = sum_stats['CortThick_case_vs_controls_adult']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['ocd'] = ct_surf

# schizophrenia
sum_stats = datasets.load_summary_stats('schizophrenia')
CT = sum_stats['CortThick_case_vs_controls']
CT = CT['d_icv']
ct_surf = parcel_to_surface(CT, target_lab='aparc_fsa5')
d['schizophrenia'] = ct_surf

fsavg = nilearn.datasets.fetch_surf_fsaverage('fsaverage5')
sphere_lh_fs = surface.load_surf_mesh(fsavg['sphere_left'])[0]
sphere_rh_fs = surface.load_surf_mesh(fsavg['sphere_right'])[0]
sp_fs = SpinPermutations(n_rep=1000, random_state=0)
sp_fs.fit(sphere_lh_fs, points_rh=sphere_rh_fs)


def permute(inp):
    left = inp[:10242]
    right = inp[10242:]
    return np.hstack(sp_fs.randomize(left, right))


g1 = schaefer_to_surf(parcels, rc_g1)
g2 = schaefer_to_surf(parcels, rc_g2)
g3 = schaefer_to_surf(parcels, rc_g3)

enig_g1_spin = {}
enig_g2_spin = {}
enig_g3_spin = {}
select = ['asd', 'ocd', 'schizophrenia', 'allepilepsy', 'depression_all', '22q', 'adhd', 'bipolar']
for key, val in d.items():
    if key in select:
        print(key)
        permuted = permute(val)
        enig_g1_spin[key] = spin(g1, val, permuted)
        enig_g2_spin[key] = spin(g2, val, permuted)
        enig_g3_spin[key] = spin(g3, val, permuted)


def plot_corr(inp, fname):
    enigma_keys = {'DGS': inp['22q'],
                   'ADHD': inp['adhd'],
                   'ASD': inp['asd'],
                   'BPD': inp['bipolar'],
                   'EPS': inp['allepilepsy'],
                   'MDD': inp['depression_all'],
                   'OCD': inp['ocd'],
                   'SCZ': inp['schizophrenia']}
    plot = pd.DataFrame.from_dict(enigma_keys, orient='index')
    plot.columns = ['r', 'p_vgm']
    plot.sort_values(by='r', inplace=True)
    fig, ax = plt.subplots(1, figsize=[11, 12])
    color = ['lightskyblue' if x > 0.05 else 'dodgerblue' for x in plot['p_vgm']]
    ax.barh(range(len(plot)), plot['r'], color=color)
    ax.set_yticks(range(0, len(plot)))
    ax.set_yticklabels(plot.index, fontsize=32)
    ax.set_xlabel("spearman's r", fontsize=32)
    ax.tick_params(labelsize=32)
    plt.xlim([-0.50, 0.50])
    plt.tight_layout()
    fig.savefig(res_p + fname)
    return


plot_corr(enig_g1_spin, 'rc_g1_enigma.png')
plot_corr(enig_g2_spin, 'rc_g2_enigma.png')
plot_corr(enig_g3_spin, 'rc_g3_enigma.png')

# functional activation

'''
This script correlates meta-analytical terms with significance testing.
'''
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import nibabel as nib
from nilearn.input_data import NiftiLabelsMasker
from nilearn import datasets
import re
import os
import argparse

path=input_path
rc_g1=np.load(path + 'rc_g1_100.npy')
rc_g2=np.load(path + 'rc_g2_100.npy')
rc_g3=np.load(path + 'rc_g3_100.npy')
rc_g1_vgm=np.load(path + 'rc_g1_100_vgm.npy')
rc_g2_vgm=np.load(path + 'rc_g2_100_vgm.npy')
rc_g3_vgm=np.load(path + 'rc_g3_100_vgm.npy')

grads=[rc_g1, rc_g2, rc_g3]
grads_vgm=[rc_g1_vgm, rc_g2_vgm, rc_g3_vgm]

dataset=datasets.fetch_atlas_schaefer_2018(n_rois=100, yeo_networks=7)
mask = NiftiLabelsMasker(dataset['maps'], resampling_target='data', strategy='mean').fit()

names=[]
path=path + r"lda50_images\\"
select=os.listdir(path)
def vgm_matcher(inp:list, inp_vgm:list):
    n_rand = 1000
    correlations=[]
    p=[]
    for i in range(len(inp)):
        correlations.append(np.zeros(50))
        p.append(np.zeros(50))
    for idx, val in enumerate(select):
        names.append(val[1])
        feature_data = nib.load(path+val[0])
        trafo=mask.transform(feature_data).squeeze()
        for i in range(len(inp)):
            grad = inp[i]
            grad_vgm = inp_vgm[i]
            r_obs=spearmanr(grad, trafo)[0]
            correlations[i][idx]=r_obs
            r_spin=np.empty(n_rand)
            for j, perm in enumerate(grad_vgm):
                r_spin[j]= spearmanr(perm, trafo)[0]
            p[i][idx]=np.mean(np.abs(r_spin) >= np.abs(r_obs))

    for i in range(len(inp)):
        df=pd.DataFrame({"Spearman's r": correlations[i], 'p_vgm': p[i]}, index=names)
        df.to_csv(res_p + r'g{}_{}_nsynth_vgm.csv'.format(i+1,100))
        
vgm_matcher(grads, grads_vgm)

