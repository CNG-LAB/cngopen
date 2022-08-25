#!/usr/bin/env python3
"""
This script contains the code behind the results in F1 in manuscript xxx
"""
import numpy as np
from scipy.stats import spearmanr

def vgm(inp, g1, g1_vgm):
    n_rand=1000
    r_vgm = np.empty(n_rand)
    r_obs = spearmanr(g1, inp)[0]
    # Compute perm pval
    for i, perm in enumerate(g1_vgm):
        r_vgm[i] = spearmanr(perm, inp)[0]
    pv_vgm = np.mean(np.abs(r_vgm) >= np.abs(r_obs))
    return [r_obs, pv_vgm]

def spin(inp, grad, spin_grad):
    r_spin = np.empty(1000)
    mask = ~np.isnan(grad) & ~ np.isnan(inp)
    r_obs = spearmanr(grad[mask], inp[mask])[0]
    for i, perm in enumerate(spin_grad):
        mask_rot = mask & ~np.isnan(perm)  # Remove midline
        r_spin[i] = spearmanr(perm[mask_rot], inp[mask_rot])[0]
    pv_spin = np.mean(np.abs(r_spin) >= np.abs(r_obs))
    return [r_obs, pv_spin]

def schaefer_to_surf(parcels, inp):
    table=np.genfromtxt(r'C:\Users\benja\Documents\md_thesis\data\mappings\schaefer_{}_fsa5.csv'.format(parcels), delimiter=',')
    lookup={}
    for ind, val in enumerate(inp):
        lookup[ind]=val
    return np.array([lookup[int(x)] if x != -1 else np.NaN for x in table])

def spearman(x, y):
    """
    :param x: One matrix
    :param y: Second matrix
    :return: matrix of spearman correlation coefficients

    This function calculates spearman rank correlations between two matrices in a vectorized manner.

    Adapted from https://cancerdatascience.org/blog/posts/pearson-correlation/
    """
    import numpy as np
    x = np.argsort(np.argsort(x, axis=0), axis=0)
    y = np.argsort(np.argsort(y, axis=0), axis=0)
    xv = x - x.mean(axis=0)
    yv = y - y.mean(axis=0)
    xvss = (xv * xv).sum(axis=0)
    yvss = (yv * yv).sum(axis=0)
    result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
    # bound the values to -1 to 1 in the event of precision issues
    return np.maximum(np.minimum(result, 1.0), -1.0)