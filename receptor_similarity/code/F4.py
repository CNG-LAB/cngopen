#!/usr/bin/env python3
"""
This script contains the code behind the results in F4 in manuscript xxx
"""
from scipy.stats import zscore, spearmanr
import numpy as np
import pandas as pd
from brainspace.gradient import GradientMaps
from mn_funcs import vgm
import seaborn as sns
import matplotlib.pyplot as plt
import os
from brainspace.datasets import load_parcellation
from palettable.cartocolors.qualitative import Vivid_4

parcels = 100
input_path = 'path/to/data/'
res_p = 'path/to/results/F4/'


# generate data
def data_generator(size):
    path = input_path
    RC = pd.read_csv(path + '{}Parcels7Networks_receptorprofiles.csv'.format(size), index_col=0)
    RC = RC.apply(zscore)
    RC = RC.transpose().corr('spearman')
    FC = np.load(path + 'mics_fc_{}.npy'.format(size))
    SC = np.load(path + 'mics_sc_{}.npy'.format(size))
    MPC = np.load(path + 'mics_mpc_{}.npy'.format(size))

    labeling = load_parcellation('schaefer', scale=size, join=True)

    ntrm_grad = GradientMaps(approach='dm', kernel='normalized_angle', random_state=1)
    ntrm_grad.fit(RC.values)
    gm_SC_L = GradientMaps(approach='dm', kernel='normalized_angle', random_state=1)
    gm_SC_L.fit(SC[0:int(size / 2), 0:int(size / 2)], sparsity=0.9)
    # SC Right hemi
    gm_SC_R = GradientMaps(alignment='procrustes', kernel='normalized_angle',
                           random_state=1)
    gm_SC_R.fit(SC[int(size / 2):size, int(size / 2):size], sparsity=0.9, reference=gm_SC_L.gradients_)
    sc_grad = np.concatenate((gm_SC_L.gradients_, gm_SC_R.gradients_), axis=0)

    fc_grad = GradientMaps(approach='dm', kernel='normalized_angle', random_state=1)
    fc_grad.fit(FC)

    mpc_grad = GradientMaps(approach='dm', kernel='normalized_angle', random_state=1)
    mpc_grad.fit(MPC)

    bigbrain_g = np.genfromtxt(
        r'C:\Users\benja\Anaconda3\envs\sklearn-env\Lib\site-packages\enigmatoolbox-1.1.3-py3.8.egg\enigmatoolbox\histology\bb_gradient_schaefer_{}.csv'.format(
            size), delimiter=',')

    return {'labeling': labeling, 'NTRM': ntrm_grad.gradients_, 'FC': fc_grad.gradients_, 'MPC': mpc_grad.gradients_,
            'SC': sc_grad, 'BB': bigbrain_g}


g100 = data_generator(parcels)

g1_vgm = np.load(input_path + 'rc_g1_{}_vgm.npy'.format(parcels))
g2_vgm = np.load(input_path + 'rc_g2_{}_vgm.npy'.format(parcels))
g3_vgm = np.load(input_path + 'rc_g3_{}_vgm.npy'.format(parcels))


# gradient correlations

def correlate(dic):
    rc_g1 = dic['NTRM'][:, 0]
    rc_g2 = dic['NTRM'][:, 1]
    rc_g3 = dic['NTRM'][:, 2]

    fc_g1 = dic['FC'][:, 0]
    fc_g2 = dic['FC'][:, 1]
    sc_g1 = dic['SC'][:, 0]
    sc_g2 = dic['SC'][:, 1]
    mpc_g1 = dic['MPC'][:, 0]
    bigbrain_g = dic['BB']

    g1_res = []
    g2_res = []
    g3_res = []

    l = [sc_g1, fc_g1, mpc_g1, bigbrain_g, sc_g2, fc_g2]
    for i in l:
        g1_res.append(vgm(i, rc_g1, g1_vgm))
        g2_res.append(vgm(i, rc_g2, g2_vgm))
        g3_res.append(vgm(i, rc_g3, g3_vgm))
    arr = np.array([g1_res, g2_res, g3_res])
    return arr


corr = correlate(g100)

corrs = corr[:, :, 0]
sigs = corr[:, :, 1]
corrs_r = np.around(corrs, 2)
corrs_s = np.where(sigs < 0.05, corrs_r, '')
stars = np.where(sigs < 0.05, '*\n', '')
ann = np.core.defchararray.add(stars, corrs_s)

fig, ax = plt.subplots(figsize=(30, 14))
sns.heatmap(np.abs(corrs), vmin=0, vmax=1, annot=ann, fmt='', annot_kws={'fontsize': 48},
            cbar_kws={}, cmap='summer_r')
ax.set_yticklabels(['RC G1', 'RC G2', 'RC G3'], fontsize=42)
ax.set_xticklabels(['SC G1', 'FC G1', 'MPC G1', 'BB G1', 'SC G2', 'FC G2'], fontsize=42)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=38)
plt.xticks(rotation=0)
plt.yticks(rotation=90)
plt.tight_layout()
fig.savefig(res_p + 'corr_heatmap.png')


# coupling

def coupler(size):
    path = r'C:\Users\benja\Documents\md_thesis\data\stable\\'
    RC = pd.read_csv(path + '{}Parcels7Networks_receptorprofiles.csv'.format(size), index_col=0)
    RC = RC.apply(zscore)
    RC = RC.transpose().corr('spearman')
    RC_val = RC.values
    FC = np.load(path + 'mics_fc_{}.npy'.format(size))
    SC = np.load(path + 'mics_sc_{}.npy'.format(size))
    MPC = np.load(path + 'mics_mpc_{}.npy'.format(size))

    labeling = load_parcellation('schaefer', scale=size, join=True)

    fc_cpl = np.empty((size))
    sc_cpl = np.empty((size))
    mpc_cpl = np.empty((size))
    for i in range(size):
        mask = RC_val[i] != 0
        full_m = (mask) & (FC[i] != 0)
        fc_cpl[i] = spearmanr(RC_val[i][full_m], FC[i][full_m])[0]
        full_m = (mask) & (SC[i] != 0)
        sc_cpl[i] = spearmanr(RC_val[i][full_m], SC[i][full_m])[0]
        full_m = (mask) & (MPC[i] != 0)
        mpc_cpl[i] = spearmanr(RC_val[i][full_m], MPC[i][full_m])[0]

    return {'labeling': labeling, 'FC': fc_cpl, 'MPC': mpc_cpl, 'SC': sc_cpl}


cpl_100 = coupler(parcels)

types_100 = np.load(r'C:\Users\benja\Documents\md_thesis\data\stable\mesulam_100parc.npy')

col_types = {'paralimbic': Vivid_4.mpl_colors[3], 'unimodal': Vivid_4.mpl_colors[1],
             'idiotypic': Vivid_4.mpl_colors[0], 'heteromodal': Vivid_4.mpl_colors[2]}


def cpl_types(modality, fname):
    plot = pd.DataFrame(cpl_100[modality], types_100)
    plot.reset_index(inplace=True)
    plot.columns = ['Laminar Differentiation', modality]
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(y=modality, x='Laminar Differentiation', data=plot, order=['idiotypic', 'unimodal',
                                                                           'heteromodal', 'paralimbic'],
                palette=col_types
                , width=0.65)
    ax.set_ylabel('RC-{} Coupling'.format(modality), fontsize=28)
    ax.set_xlabel('')
    ax.tick_params(labelsize=26)
    ax.set_xticklabels(['idio-\ntypic', 'uni-\nmodal', 'hetero-\nmodal', 'para-\nlimbic'])
    sns.despine(trim=True, offset=10)
    plt.tight_layout()
    fig.savefig(res_p + fname)
    return


cpl_types('FC', 'RC_FC_cpl.png')
cpl_types('SC', 'RC_SC_cpl.png')
cpl_types('MPC', 'RC_MPC_cpl.png')

# modular stability
networks = np.load(input_path + 'mesulam_7_100.npy')
path = input_path + 'c_partition/igraph/schaefer100/'
files = os.listdir(path)
d = []
for i in range(len(files)):
    temp = np.load(path + files[i], allow_pickle=True).item()
    d.append(temp)
parts = pd.DataFrame(d)
parts.sort_values(by='gamma', inplace=True)
label_arr = np.array(networks)
score = {}
sizes = {}
for i in set(networks):
    sizes[i] = sum(label_arr == i) / 100
resolutions = []
for i in d:
    part = i['partition']
    resolutions.append(i['gamma'])
    num_parts = len(set(part))
    for label in set(networks):
        ind = label_arr == label
        subset = part[ind]
        uniques = np.unique(subset, return_counts=True)
        biggest_part = max(uniques[1]) / sum(uniques[1])
        parts_in_network = len(uniques[0])
        parts_outside_network = num_parts - parts_in_network
        network_size = sizes[label]
        mds = (biggest_part / (parts_in_network / num_parts)) * network_size
        if label not in score.keys():
            score[label] = [mds]
        else:
            score[label].append(mds)

scoreframe = pd.DataFrame(score)
means = scoreframe
means['gamma'] = parts['gamma']
means.sort_values(by='gamma', inplace=True)
order = 2
size = 15

fig,ax=plt.subplots(figsize=(15,10))
idiotypic=sns.regplot(x='gamma', y='idiotypic', color=Vivid_4.mpl_colors[0], label='idiotypic', data=means, scatter=True, order=order,
                  scatter_kws={'s':size})

unimodal=sns.regplot(x='gamma', y='unimodal', color=Vivid_4.mpl_colors[1], label='unimodal', data=means, scatter=True, order=order,
                  scatter_kws={'s':size})

heteromodal=sns.regplot(x='gamma', y='heteromodal', color=Vivid_4.mpl_colors[2], label='heteromodal', data=means, scatter=True,
                          order=order, scatter_kws={'s':size})

paralimbic=sns.regplot(x='gamma', y='paralimbic', color=Vivid_4.mpl_colors[3], label='paralimbic', data=means, scatter=True,
                             order=order,
                            scatter_kws={'s':size})
lgnd = ax.legend(fontsize=34, frameon=True)
for handle in lgnd.legendHandles:
    handle.set_sizes([150])
ax.set_xlabel('gamma', fontsize=34)
ax.set_ylabel('Modular stability score', fontsize=34)
ax.tick_params(labelsize=30)
plt.tight_layout()
fig.savefig(res_p + 'modular_stability_ctypes.png')
