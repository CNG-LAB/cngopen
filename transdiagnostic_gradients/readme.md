# Welcome to the Transdiagnostic Gradients Project!

### Here you can find the code and data that was used and generated for the manuscript „COORDINATED CORTICAL THICKNESS ALTERATIONS ACROSS PSYCHIATRIC CONDITIONS: A TRANSDIAGNOSTIC ENIGMA STUDY“ (https://www.medrxiv.org/content/10.1101/2022.02.03.22270326v1.full)

## Scripts 

### Hettwer2022_Figure1_Transdiagnostic_Hubs_Epicenters.m

This script runs analyses presented in Figure 1 of the manuscript. It loads ENIGMA summary statistics (Cohen's d maps) and HCP connectivity data, computes cross-disorder co-alteration hubs as well as functional and structural disease epicenters. You can exchange selected disorders or use your own effect size maps (in DK or fsa5 space). 

Fig1.png![grafik](https://user-images.githubusercontent.com/93781179/155957059-68268733-d97a-40c8-b069-b603254e5801.png)


### Hettwer2022_Figure2_Transdiagnostic_Gradients.m

This script runs analyses presented in Figure 2 of the manuscript. It loads ENIGMA summary statistics (Cohen's d maps), computes cross-disorder similarity, derives transdiagnostic gradients and contextualizes derived gradients with cytoarchitectonic and functional data.

Fig2.png![grafik](https://user-images.githubusercontent.com/93781179/155957301-e9c621c6-73f4-4b13-bafd-25ee16875911.png)


### Dependencies
For the current work, we rely on the ENIGMA Toolbox (https://enigma-toolbox.readthedocs.io/en/latest/; we used ENIGMA-1.1.3), brainspace (download here: https://brainspace.readthedocs.io/en/latest/pages/install.html; used here: BrainSpace-0.1.2), surfstat: https://mica-mni.github.io/surfstat/, cbrewer https://de.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab, and NeuroSynth https://www.neurosynth.org.

We ran the code on Mac OS Monterey 12.2.1 in Matlab 2021a. No non-standard hardware is required.

Expected run time: 10 minutes
### Required input
If you'd like to run the code on your own data, you can either make use of the ENIGMA Toolbox to load/include other disorders or upload your own data. If you aim to use your own data, it should be case-control differences (Cohen's d maps) on a template implemented in the ENIGMA Toolbox (such as Desikan-Killiany Parcels, fsa5, or Conte69).

### Expected output
After running both scripts, you should have retrieved a co-alteration hub map, functional and structural disease epicenters (r and p values), 7 transdiagnostic gradients, their stratification according to cytoarchitectonic classes (i.e. mean gradient loading per class), correlating genes, and weighted z-values indicating the ordering of functional topic terms along the first 2 transdiagnostic gradients.

### Data that was generated for this study:

* _Disorder_covariance.mat_ = Cross-disorder inter-parcel correlation of COhen’s d values
* _Normative_Connectivity_Hubs.mat_ = Hubs computed based on HCP young adult sample (rs-fMRI and DTI)
* _Transdiagnostic_Covariance_Hubs.mat_ = Cross-disorder hubs of covariance of CT alterations
* _Epicenters.mat_ = All cortical and sub-cortical functional and structural disease epicenters
* _CT_psych_gradients.mat_ = the transdiagnostic gradients presented in this manuscript
* _Cyto_gradients.mat_ = The 2 trans diagnostic gradients stratified according to Von Economo-Koskinas cytoarchitectonic classes
* _Genes_G1.mat_ = Genes whose cortical expression pattern significantly correlates with G1

### Data that is required to run the analyses
* _Cohen’s_d_CT_6_disorders.mat_ = ENIGMA case-control difference Cohen’s d maps for cortical thickness.
* _Valk2020_CT_structural_covariance_gradient.mat_ = normative Gradient of cortical thickness covariance, published by Sofie Valk in Science Advances 2020
* _neurosynth_z_values_ = Functional topic terms and their fMRI correlates (maps)
* _DK_midbrain_parcels.mat_ / midbrain_vertices_fsa5.mat = vertices which you can set to 0 during visualization to mask the midbrain
