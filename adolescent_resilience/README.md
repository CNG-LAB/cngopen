# Dynamic change in adolescent resilience and vulnerability to environmental adversity links to cortical myelination trajectories
## This repository contains scripts and data to reproduce analyses preprinted here: https://osf.io/preprints/psyarxiv/2dv68

## Instructions

If you'd like to reproduce the presented analyses, you can run these 3 scripts in the following order:

- **1. distress_prediction.ypnb:** This Jupyter notebook generates resilience scores by predicting well-being / distress from adversity measures. As features, sumscores from the following questionnaires are required: The Life Events Questionnaire (LEQ), Child Trauma Questionnaire (CTQ), Alabama Parenting Questionnaire (APQ), Measure of Parenting Style (MOPS), and socioeconomic status (as approximated by zip codes/IMD). 
    _Expected outputs_: Resilience scores  resilience_scores_prediction_repeated_measures.csv
      _Expected runtime_: ~1 day
  
  <img width="533" alt="Bildschirmfoto 2024-02-21 um 09 46 12" src="https://github.com/CNG-LAB/cngopen/assets/93781179/634dafd1-0021-4139-9761-1a22d4a7c556">

- **MT_delta_Figures_1_2.m:** This script computes change in MT and FC and their relation to changes in stressor resilience scores. _Expected outputs_:  Figures 1&2.
  _Expected runtime_: ~40 minutes
  
  <img width="1253" alt="Bildschirmfoto 2024-02-21 um 09 48 54" src="https://github.com/CNG-LAB/cngopen/assets/93781179/15ae6902-6b6d-4711-b47c-0d519671c06d">

- **Maturational_Index_Gradients_Figures_3_4.m:** This script computes MPC and FC maturational indices, as well as MPC gradients. It compares groups of individuals showing an increase vs. decrease in resilience capacities, and compares this group difference to nullmodels of random/permuted groups.
   _Expected outputs_: Figures 3 and 4.
  _Expected runtime_: On a cluster: ~1-2 days, on a desktop computer: ~3 days.
  
  <img width="1349" alt="Bildschirmfoto 2024-02-21 um 09 50 15" src="https://github.com/CNG-LAB/cngopen/assets/93781179/c8186bb8-ad50-452b-8f3a-fb84d1b3b1e4">


## Data:
This study is based on data from the Neuroscience in Psychiatry Network. You can find more information and apply for access here: https://portal.ide-cam.org.uk/overview/6
- **NSPN_mt_subsample.mat:** Preprocessed MT profiles for all participants. Functional connectivity data can be found here: https://zenodo.org/records/6390852. 
- **resilience_scores_prediction_repeatedmeasures_totalscores.csv:** Resilience scores derived from the prediction.
- **MI_FC.mat & MI_MPC.mat:** Maturational index based on the full imaging sample.
- **MI_delta_groups_MPC/FC_res_groups.mat:** Group difference map for each maturational index.
- **mpc_diff_ROIs.csv:** Group difference results table.
- Other uploaded data includes colormaps, cortical types and yeo atlases, permutation indices.

## Requirements 
Matlab scripts were run in Matlab 2022b, Python scripts were run in Python 3.10. 

_Matlab_: Our analysis code makes use of open software: Gradient mapping analyses were carried out using BrainSpace (v. 0.1.2;  https://brainspace.readthedocs.io/en/latest/) and surface visualizations were based on code from the ENIGMA Toolbox (v.1.1.3; https://enigma-toolbox.readthedocs.io/en/latest/; 104) in combination with ColorBrewer (v. 1.0.0; https://github.com/scottclowe/cbrewer2). Statistical analyses were carried out using SurfStat (https://www.math.mcgill.ca/keith/surfstat/). Equivolumetric surfaces were computed using code from: https://github.com/MICA-MNI/micaopen/tree/master/a_moment_of_change. Z-tests were performed using the compare correlation coefficients function (Sisi Ma (2024), compare_correlation_coefficients (https://www.mathworks.com/matlabcentral/fileexchange/44658-compare_correlation_coefficients). 

_Python_:
We made use of the following packages: scipy 1.10.1, sklearn 0.0.post1, matplotlib 3.7.1, numpy 1.24.2, pandas 1.5.3, seaborn 0.11.0

The scripts can be run on a standard desktop computer, however, the computation of Maturational Indices and the respective permutations will likely take 3-4 days. We recommend running these computations on a cluster to parallelize permutations.

## Install
Follow these instructions to install Matlab (license required): https://de.mathworks.com/help/install/, then clone this repository:
git clone https://github.com/CNG-LAB/cngopen/tree/main/adolescent_resilience/
(This should take a few seconds only.)
