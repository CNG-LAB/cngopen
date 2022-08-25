#!/usr/bin/env python3
"""
This script contains the code behind the results in F3 in manuscript xxx
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

g1 = np.load(input_path + 'rc_g1_{}.npy'.format(parcels))
g2 = np.load(input_path + 'rc_g2_{}.npy'.format(parcels))
g3 = np.load(input_path + 'rc_g3_{}.npy'.format(parcels))

g1_vgm = np.load(input_path + 'rc_g1_{}_vgm.npy'.format(parcels))
g2_vgm = np.load(input_path + 'rc_g2_{}_vgm.npy'.format(parcels))
g3_vgm = np.load(input_path + 'rc_g3_{}_vgm.npy'.format(parcels))

grads = [g1, g2, g3]
grads_vgm = [g1_vgm, g2_vgm, g3_vgm]

path_neurosynth = '/path/to/brainstat/data/'
nii_files = os.listdir(path_neurosynth)
dataset = nilearn.datasets.fetch_atlas_schaefer_2018(n_rois=parcels, yeo_networks=7)
mask = NiftiLabelsMasker(dataset['maps'], resampling_target='data', strategy='mean').fit()
names = []
pat = re.compile('__([A-Za-z0-9 ]+).+z$')


def vgm_matcher(inp: list, inp_vgm: list):
    n_rand = 1000
    correlations = []
    p = []
    for i in range(len(inp)):
        correlations.append(np.zeros(len(nii_files)))
        p.append(np.zeros(len(nii_files)))
    for idx, val in enumerate(nii_files):
        names.append(re.search(pat, val)[1])
        feature_data = nib.load(path_neurosynth + val)
        trafo = mask.transform(feature_data).squeeze()
        for i in range(len(inp)):
            grad = inp[i]
            grad_vgm = inp_vgm[i]
            r_obs = spearmanr(grad, trafo)[0]
            correlations[i][idx] = r_obs
            r_spin = np.empty(n_rand)

            for j, perm in enumerate(grad_vgm):
                r_spin[j] = spearmanr(perm, trafo)[0]
            p[i][idx] = np.mean(np.abs(r_spin) >= np.abs(r_obs))

    for i in range(len(inp)):
        df = pd.DataFrame({"Spearman's r": correlations[i], 'p_vgm': p[i]}, index=names)
        df.sort_values(by="Spearman's r", inplace=True)
        df.to_csv(input_path + 'g{}_{}_nsynth_vgm.csv'.format(i + 1, parcels))
    return


vgm_matcher(grads, grads_vgm)

n_g1 = pd.read_csv(input_path + 'g1_{}_nsynth_vgm.csv'.format(parcels), index_col=0)
n_g2 = pd.read_csv(input_path + 'g2_{}_nsynth_vgm.csv'.format(parcels), index_col=0)
n_g3 = pd.read_csv(input_path + 'g3_{}_nsynth_vgm.csv'.format(parcels), index_col=0)

with open(input_path + 'supp_l1.txt', 'r') as l:
    interest = l.readlines()

g1_interest = n_g1.transpose()[interest].transpose()
g2_interest = n_g2.transpose()[interest].transpose()
g3_interest = n_g3.transpose()[interest].transpose()

# g1
g1_interest = g1_interest[g1_interest['p_vgm'] < 0.05]
keep_g1 = ['face recognition', 'autobiographical memory', 'mind tom', 'secondary somatosensory', 'parkinson',
           'primary somatosensory','coordination', 'motor imagery', 'painful', 'response selection', 'anticipation',
           'inhibitory control','executive functions', 'control network', 'adhd', 'goal directed', 'insight',
           'focusing', 'illusion','manipulation', 'interference']

total = g1_interest.transpose()[keep_g1].transpose()
for_treemap = total.reset_index()
for_treemap['abs corr'] = np.abs(for_treemap["Spearman's r"])
test = [float("{:.2f}".format(x)) for x in for_treemap['abs corr']]
for_treemap['abs corr'] = test
breaks = [x.replace(' ', '<br>') for x in list(for_treemap['index'])]
for_treemap['breaks'] = breaks
for_treemap.iloc[1, 4] = 'autobio-<br>graphical<br>memory'
for_treemap.iloc[17, 4] = 'focu-<br>sing'

fig = px.treemap(for_treemap, names='breaks', path=['breaks'], values="abs corr", color="Spearman's r",
                 range_color=[-0.7, 0.7], color_continuous_scale='rdbu_r')
fig.update_layout(
    uniformtext=dict(minsize=28, mode='show'), margin=dict(t=50, l=25, r=25, b=25), width=800, height=1000,
    font=dict(family='arial'),
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)')

fig.update_layout({'plot_bgcolor': 'rgba(255, 255, 255, 255)', 'paper_bgcolor': 'rgba(255, 255, 255, 255)', })
fig.update_coloraxes(colorbar_orientation='h', colorbar_thickness=10,
                     colorbar_title_font_size=10, colorbar_tickfont_size=10)
fig.write_image(res_p + 'g1_nsynth.png')

# g2
g2_interest = g2_interest[g2_interest['p_vgm'] < 0.05]
keep_g2 = ['primary visual', 'executive functions', 'control network', 'response selection', 'belief', 'motor pre',
           'intention', 'painful', 'intelligence', 'working memory', 'disorder ocd', 'parkinson', 'motor sma',
           'uncertainty', 'goal', 'subtraction', 'judge', 'rules', 'secondary somatosensory', 'primary somatosensory',
           'efficiency', 'risky', 'decision making', 'attention', 'cognition', 'thoughts', 'bipolar disorder', 'theory',
           'self', 'mental state', 'beliefs', 'mind','planning']

total = g2_interest.transpose()[keep_g2].transpose()
for_treemap = total.reset_index()
for_treemap['abs corr'] = np.abs(for_treemap["Spearman's r"])
test = [float("{:.2f}".format(x)) for x in for_treemap['abs corr']]
for_treemap['abs corr'] = test
breaks = [x.replace(' ', '<br>') for x in list(for_treemap['index'])]
for_treemap['breaks'] = breaks
for_treemap.iloc[19, 4] = 'secondary<br>somato-<br>sensory'
for_treemap.iloc[20, 4] = 'primary<br>somato-<br>sensory'
for_treemap.iloc[25, 4] = 'cog-<br>nition'
for_treemap.iloc[24, 4] = 'atten-<br>tion'
fig = px.treemap(for_treemap, names='breaks', path=['breaks'], values="abs corr", color="Spearman's r",
                 range_color=[-0.7, 0.7], color_continuous_scale='rdbu_r')
fig.update_layout(
    uniformtext=dict(minsize=28, mode='show'), margin=dict(t=50, l=25, r=25, b=25), width=800, height=1000,
    font=dict(family='arial'),
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)')

fig.update_layout({'plot_bgcolor': 'rgba(255, 255, 255, 255)', 'paper_bgcolor': 'rgba(255, 255, 255, 255)', })
fig.update_coloraxes(colorbar_orientation='h', colorbar_thickness=10,
                     colorbar_title_font_size=10, colorbar_tickfont_size=10)
fig.write_image(res_p + 'g2_nsynth.png')

# g3
g3_interest = g3_interest[g3_interest['p_vgm'] < 0.05]
keep_g3 = ['early visual', 'navigation', 'visual attention', 'empathy', 'social cognitive', 'primary auditory',
           'listening', 'consolidation', 'dementia']
total = g3_interest.transpose()[keep_g3].transpose()
for_treemap = total.reset_index()
for_treemap['abs corr'] = np.abs(for_treemap["Spearman's r"])
test = [float("{:.2f}".format(x)) for x in for_treemap['abs corr']]
for_treemap['abs corr'] = test
breaks = [x.replace(' ', '<br>') for x in list(for_treemap['index'])]
for_treemap['breaks'] = breaks
for_treemap.iloc[19, 4] = 'secondary<br>somato-<br>sensory'
for_treemap.iloc[20, 4] = 'primary<br>somato-<br>sensory'
for_treemap.iloc[25, 4] = 'cog-<br>nition'
for_treemap.iloc[24, 4] = 'atten-<br>tion'
fig = px.treemap(for_treemap, names='breaks', path=['breaks'], values="abs corr", color="Spearman's r",
                 range_color=[-0.7, 0.7], color_continuous_scale='rdbu_r')
fig.update_layout(
    uniformtext=dict(minsize=28, mode='show'), margin=dict(t=50, l=25, r=25, b=25), width=800, height=1000,
    font=dict(family='arial'),
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)')

fig.update_layout({'plot_bgcolor': 'rgba(255, 255, 255, 255)', 'paper_bgcolor': 'rgba(255, 255, 255, 255)', })
fig.update_coloraxes(colorbar_orientation='h', colorbar_thickness=10,
                     colorbar_title_font_size=10, colorbar_tickfont_size=10)
fig.write_image(res_p + 'g3_nsynth.png')
