#!/usr/bin/env python3
"""
This script contains the code behind the results in F1 in manuscript xxx
"""

from scipy.stats import zscore
import numpy as np
from nilearn import datasets, plotting
from nilearn.input_data import NiftiMasker
from palettable.colorbrewer.qualitative import Set2_4, Set2_6
import os
import pandas as pd
from brainspace.datasets import load_parcellation, load_conte69
from brainspace.gradient import GradientMaps
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels
import seaborn as sns
import matplotlib.pyplot as plt
import ptitprince as pt
from mn_funcs import spearman
from palettable.scientific.diverging import Roma_4

input_path = 'path/to/data/'
res_p = 'path/to/results/F1/'
# load data


mask = NiftiMasker(
    r'C:/Users/benja/Documents/md_thesis/data/tian_subcortex/masks/S1_3T/tian_binary_mask_total.nii.gz').fit()

# hierarchical clustering
subcortex_regions = pd.read_csv(input_path + 'subcortex_tian_regions.csv', index_col=0)
sns.set(font_scale=2.3)
sns.clustermap(subcortex_regions, method='weighted', cmap='coolwarm', figsize=(15, 10), vmin=-2.4, vmax=2.4,
               tree_kws={'linewidths': 2})
plt.tight_layout()
plt.savefig(res_p + 'subcortex_hierarchical.png')

# subcortical gradients
# make and subcortical display
subcortex = pd.read_csv(input_path + 'tian_subcortex.csv', index_col=0)
subcortex = subcortex.apply(zscore)
sc_corr = subcortex.transpose().corr('spearman')
gm = GradientMaps(n_components=10, approach='dm', kernel='normalized_angle')
gm.fit(sc_corr.values)
s_g1 = gm.gradients_[:, 0]

sg1 = mask.inverse_transform(s_g1)
bg = datasets.load_mni152_template(1)

fig = plt.figure(figsize=(14, 7))
plotting.plot_stat_map(sg1, bg_img=bg, draw_cross=False, display_mode='ortho',
                       black_bg=False, cmap='viridis_r', annotate=False, figure=fig, colorbar=False,
                       output_file=res_p + 'subcort_grad_1_proj.png', dim=2, cut_coords=(12, -4, 2))
# scree
scat = [(x / sum(gm.lambdas_)) * 100 for x in gm.lambdas_]
fig, ax = plt.subplots(figsize=(8, 8))
sns.set_style('ticks')
ax.scatter(range(1, len(scat) + 1), scat, s=140)
ax.plot(range(1, len(scat) + 1), scat, '--')
ax.set_xlabel('# component', fontsize=28)
ax.set_xticks(range(1, 11))
ax.set_ylabel('% variance explained', fontsize=28)
ax.tick_params(labelsize=24)
plt.tight_layout()
sns.despine()
fig.savefig(res_p + 'sc_grad_scree.png')

# raincloud plot
masks = os.listdir(input_path + 'masks/S1_3T/')
del masks[-1]
del masks[-7]
del masks[-7]
masks.sort()
with open(input_path + 'Tian_Subcortex_S1_3T_label.txt') as f:
    s = f.readlines()
    s = [x.strip() for x in s]
s.sort()
del s[-9]
del s[-9]
path = input_path + 'masks/S1_3T/'
ct = NiftiMasker(path + 'tian_binary_mask_total.nii.gz').fit()
s_grad_total = sg1
res_sgrad = {}
for idx, mask in enumerate(masks):
    roi_mask = NiftiMasker(path + mask).fit()
    arr = roi_mask.transform(s_grad_total)
    res_sgrad[s[idx]] = arr.squeeze()
labels_s = []
values = []
for key, val in res_sgrad.items():
    n_key = [key for x in val]
    labels_s.extend(n_key)
    values.extend(list(val))
s_grad_df = pd.DataFrame({'Region': labels_s, 'Value': values})
alt = list(s_grad_df['Region'])
alt = [x.strip('-lh') for x in alt]
alt = [x.strip('-rh') for x in alt]
alt = ['THA' if x in ('pTHA', 'aTHA') else x for x in alt]
s_grad_df['alternate lable'] = alt
s_grad_df.sort_values(by='alternate lable', inplace=True)
order = list(set(alt))
order.sort()
r = []
for i in order:
    sub = s_grad_df[s_grad_df['alternate lable'] == i]
    mean = float(sub.median())
    s = [mean] * len(sub)
    r.extend(s)
s_grad_df['mean alt'] = r
s_grad_df.sort_values(by='mean alt', inplace=True)
f, ax = plt.subplots(figsize=(8, 8))
pal = Set2_6.mpl_colors
dy = "Value"
dx = "alternate lable"
ort = "h"
sigma = 0.2
pt.RainCloud(x=dx, y=dy, data=s_grad_df, palette=pal, bw=sigma,
             width_viol=1.35, ax=ax, orient=ort)
ax.set_ylabel('', fontsize=14)
ax.set_xlabel('Principal Gradient 1', fontsize=26)
ax.tick_params(labelsize=24)
plt.tight_layout()
f.savefig(res_p + 'sg1_raincloud.png')

# subcortico-cortical
cortex = pd.read_csv(r'C:\Users\benja\Documents\md_thesis\data\stable\100Parcels7Networks_receptorprofiles.csv',
                     index_col=0)
cortex = cortex.apply(zscore)
sxc = spearman(cortex.transpose().values, subcortex.transpose().values)
gm_sxc = GradientMaps(approach='dm', kernel='normalized_angle')
gm_sxc.fit(sxc)
surf_lh, surf_rh = load_conte69()
labeling = load_parcellation('schaefer', scale=100, join=True)
grad = map_to_labels(gm_sxc.gradients_[:, 0], labeling, mask=labeling != 0, fill=np.nan)
plot_hemispheres(surf_lh, surf_rh, array_name=grad, size=(750, 600), color_bar=False, embed_nb=True,
                 layout_style='grid', zoom=1.2,
                 cmap='viridis_r', screenshot=True, filename=res_p + 'sxc_g1.png')
grad = map_to_labels(gm_sxc.gradients_[:, 1], labeling, mask=labeling != 0, fill=np.nan)
plot_hemispheres(surf_lh, surf_rh, array_name=grad, size=(750, 600), color_bar=False, embed_nb=True,
                 layout_style='grid', zoom=1.2,
                 cmap='viridis_r', screenshot=True, filename=res_p + 'sxc_g2.png')
grad = map_to_labels(gm_sxc.gradients_[:, 2], labeling, mask=labeling != 0, fill=np.nan)
plot_hemispheres(surf_lh, surf_rh, array_name=grad, size=(750, 600), color_bar=False, embed_nb=True,
                 layout_style='grid', zoom=1.2,
                 cmap='viridis_r', screenshot=True, filename=res_p + 'sxc_g3.png')

