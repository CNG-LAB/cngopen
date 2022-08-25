#!/usr/bin/env python3
"""
This script contains the code behind the results in F1 in manuscript xxx
"""
from scipy.stats import zscore
from brainspace.gradient import GradientMaps
from brainspace.datasets import load_parcellation, load_conte69
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mn_funcs import vgm
from brainspace.null_models import SurrogateMaps
import seaborn as sns

parcels=100
labeling = load_parcellation('schaefer', scale=parcels, join=True)
surf_lh, surf_rh = load_conte69()

input_path='path/to/data/'
res_p='path/to/results/F1/'
#NTRM gradients and scree
ntrm=pd.read_csv(input_path + '{}Parcels7Networks_receptorprofiles.csv'.format(parcels), index_col=0)
ntrm=ntrm.apply(zscore)
surf_corr=ntrm.transpose().corr('spearman')
ntrm_grad=GradientMaps(approach='dm', kernel='normalized_angle',random_state=1)
ntrm_grad.fit(surf_corr.values)


def plot_surf(arr, fname):
    plot_hemispheres(surf_lh, surf_rh, array_name=arr, size=(1600, 450), color_bar=True,
                     cmap='viridis_r', screenshot=True, filename=res_p + fname)
    return

rc_g1=ntrm_grad.gradients_[:,0]
rc_g2=ntrm_grad.gradients_[:,1]
rc_g3=ntrm_grad.gradients_[:,2]

np.save(input_path + 'rc_g1_{}.npy'.format(parcels), rc_g1)
np.save(input_path + 'rc_g2_{}.npy'.format(parcels), rc_g2)
np.save(input_path + 'rc_g3_{}.npy'.format(parcels), rc_g3)

grad=map_to_labels(rc_g1, labeling, mask=labeling != 0, fill=np.nan)
plot_surf (grad, 'G1_on_surf.png')

grad=map_to_labels(rc_g2, labeling, mask=labeling != 0, fill=np.nan)
plot_surf (grad, 'G2_on_surf.png')

grad=map_to_labels(rc_g3, labeling, mask=labeling != 0, fill=np.nan)
plot_surf (grad, 'G3_on_surf.png')


sns.set_style('white')
fig, ax=plt.subplots(figsize=(7,7))
var=[(s / sum(ntrm_grad.lambdas_)) * 100 for s in ntrm_grad.lambdas_]
ax.scatter(range(1, len(var) + 1), var, s=140)
ax.plot(range(1,len(var) + 1), var, '--')
ax.set_xlabel('# component', fontsize=28)
ax.set_xticks(range(1,11))
ax.set_ylabel('% variance explained', fontsize=28)
ax.tick_params(labelsize=24)
plt.tight_layout()
sns.despine()
fig.savefig(res_p +'RC_scree.png')

#gradient-receptor-correlations
r_1={}
r_2={}
r_3={}

#generate permuted brain maps
dist=np.load(input_path + 'tian_euclidean_distance.npy')
ssm=SurrogateMaps(kernel='invdist')
ssm.fit(dist)
def gen_vgm(grad):
    n_surrogate_datasets = 1000
    g_vgm=ssm.randomize(grad, n_rep=n_surrogate_datasets)
    return g_vgm

g1_vgm=gen_vgm(rc_g1)
g2_vgm=gen_vgm(rc_g2)
g3_vgm=gen_vgm(rc_g3)

np.save(input_path + 'rc_g1_{}_vgm.npy'.format(parcels), g1_vgm)
np.save(input_path + 'rc_g2_{}_vgm.npy'.format(parcels), g2_vgm)
np.save(input_path + 'rc_g3_{}_vgm.npy'.format(parcels), g3_vgm)

for i in ntrm.columns:
    sub=ntrm[i]
    r_1[i]=vgm(sub, rc_g1, g1_vgm)
    r_2[i]=vgm(sub, rc_g2, g2_vgm)
    r_3[i]=vgm(sub, rc_g3, g3_vgm)

def plot_dens_corr(inp, fname):
    df1=pd.DataFrame.from_dict(inp, orient='index')
    df1.columns=["Spearman's r", 'p']
    df1.sort_values(by="Spearman's r", inplace=True)
    fig, ax=plt.subplots(figsize=(15,5))
    color=['lightskyblue' if x > 0.05 else 'dodgerblue' for x in df1['p']]
    ax.bar(range(len(df1)), df1["Spearman's r"], color=color)
    ax.set_xticks(range(len(df1)), labels=list(df1.index))
    ax.set_ylabel("Spearman's r", fontsize=28)
    ax.tick_params(labelsize=26)
    plt.xticks(rotation=30+270)
    plt.tight_layout()
    fig.savefig(res_p + fname)
    return

plot_dens_corr(r_1, 'G1_receptors.png')
plot_dens_corr(r_2, 'G2_receptors.png')
plot_dens_corr(r_3, 'G3_receptors.png')




