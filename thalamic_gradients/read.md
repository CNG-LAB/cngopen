# Low-Dimensional Thalamocortical Structural Connectivity Patterns

### Here you can find the costum code that was used for the manuscript 
### *"A Multimodal Characterization of Low-Dimensional Thalamocortical Structural Connectivity Patterns" (2024)*
preprinted: https://doi.org/10.1101/2024.02.01.578366

![alt text](https://github.com/CNG-LAB/cngopen/blob/main/thalamic_gradients/thala_gradients.png) 

## **Scripts** 
The repo includes code for analysis and plotting. 

### 1.  Thalamocortical Structural Connectivity Gradients.

- `00_refine_thalamus_mask.py` - to manually refine the Harvard-Oxford thalamus mask (integrated in FSL)
- `01_build_struc_connectivity_matrix.py` - loads probabilistic tractography data from each subject (output of probtrackx) and build thalamocortical structural connectivity matrix
- `02_thalamic_sc_gradient.py` - uses TC structural connectivity matrix as input -> computes TC structural connectivity gradients
- `fig1/01_thalamic_gradients.ipynb` - plotting gradients
- `fig1/02_Decoding_with_Thomas_atlas_lh.ipynb` - Decoding THOMAS

### 2.  Contextualization of Gradients with Microstructure and Functional Connectivity.

- `04_create_grouplevel_thalamic_T1q_map.py` - creates grouplevel qT1 map for left and right hemisphere
- `05_warp_rsfMRI_to_MNI` - script to warp rsfMRI to MNI using ANTS and the micapipe transformation files
- `06_create_functional_connectivity_matrix.py` - uses thalamus and cortex timeseries and creates TC functional connectivity matrix
- `07_compute_fc_gradients.py` - uses TC functional connectivity matrix as input -> computes TC functional connectivity gradients
- `brainsmash_variograms` - script for spatial autocorrelation using surrogates
- `fig2/Decoding_T1q.ipynb` - qT1
- `fig2/Decoding_corematrix.ipynb` - core-matrix
- `fig2/Decoding_func_conn_gradients.ipynb` - func gradients

### 3.  Cortical Projections of Structural Connectivity Gradients and their Associations to Functional Connectivity and Structural Covariance. 

- `03_load_intensity_profiles.py` - (for structural covariance computation): loads in the cortical qT1 intensity profiles (=output of micapipes MPC modul)
- `08_structural_covariance_voxelwise_mean.py` - script to compute qT1 structural covariance
- `fig3/on_cortex.ipynb` - cortex projections
- `fig3/Yeo_violin.ipynb` - Decoding with Yeo Networks

### Supp1.  Robustness of Structural Connectivity Gradients.
- `11_Supp_compute_gradients_with_diff_threshhold.py` - Supplementary, computes TC structural connectivity gradients with matrices tresholded at different percentiles
- `fig1/thalamic_gradients-diff_thresholds.ipynb` - plot gradients

### Supp2.   Thalamocortical Structural Connectivity Gradients (RH).
- `fig1/01_thalamic_gradients.ipynb` - plotting gradients
- `fig1/02_Decoding_with_Thomas_atlas_lh.ipynb` - Decoding with THOMAS

# work in progress

### Supp3.  

### Supp4.  

### Supp5.  

### Supp6.  




- `09_Supp_groundtruth.py` - Supplementary, mean projection of sc, fc, and scov of thomas nucleus (AV, VLP, MD) on surface 
- `10_Supp_QC_scov.py` - Supplementary, crosscorrelation: correlating between left voxelwise thalamic T1q and wholebrain cortex parcel T1q ->links this to left structural connectivity gradients by correlating columns (parcels) of structural covariance with G1 and G2 gradient loadings (same for right)


- `fig2` -  Supp3
- `fig3` -  Supp4
- `Supp_QC_scov` - Supp5
- `Supp_thomas_nuclei_projection` - Supp6


