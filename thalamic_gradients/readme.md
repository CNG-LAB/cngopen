# Low-Dimensional Thalamocortical Structural Connectivity Patterns

### Here you can find the costum code that was used for the manuscript 
### *"A Multimodal Characterization of Low-Dimensional Thalamocortical Structural Connectivity Patterns" (2024)*
preprinted: https://doi.org/10.1101/2024.02.01.578366

![img](https://github.com/CNG-LAB/cngopen/blob/main/thalamic_gradients/thala_gradients.png) 

## **Data**
Our Project is based on the MICA-MICs dataset (Royer et al., 2022), which can be downloaded from: https://portal.conp.ca/dataset?id=projects/mica-mics

## **Scripts** 
The repo includes code for preprocessing, analysis and plotting. 

### Preprocessing and Tractography



### Fig 1.  Thalamocortical Structural Connectivity Gradients.

- `00_refine_thalamus_mask.py` - to manually refine the Harvard-Oxford thalamus mask (integrated in FSL)
- `01_build_struc_connectivity_matrix.py` - loads probabilistic tractography data from each subject (output of probtrackx) and build thalamocortical structural connectivity matrix
- `02_thalamic_sc_gradient.py` - uses TC structural connectivity matrix as input -> computes TC structural connectivity gradients
- `fig1/01_thalamic_gradients.ipynb` - plotting gradients
- `fig1/02_Decoding_with_Thomas_atlas_lh.ipynb` - Decoding THOMAS atlas

### Fig 2.  Contextualization of Gradients with Microstructure and Functional Connectivity.

- `04_create_grouplevel_thalamic_T1q_map.py` - creates thalamic grouplevel qT1 map for left and right hemisphere
- `05_warp_rsfMRI_to_MNI` - script to warp rsfMRI to MNI using ANTS and the micapipe transformation files
- `06_create_functional_connectivity_matrix.py` - uses thalamus and cortex timeseries and creates TC functional connectivity matrix
- `07_compute_fc_gradients.py` - uses TC functional connectivity matrix as input -> computes TC functional connectivity gradients
- `brainsmash_variograms` - script for spatial autocorrelation using surrogates
- `fig2/Decoding_T1q.ipynb` - qT1
- `fig2/Decoding_corematrix.ipynb` - core-matrix
- `fig2/Decoding_func_conn_gradients.ipynb` - func gradients

### Fig 3.  Cortical Projections of Structural Connectivity Gradients and their Associations to Functional Connectivity and Structural Covariance. 

- `03_load_intensity_profiles.py` - (for structural covariance computation): loads in the cortical qT1 intensity profiles (=output of micapipes MPC modul)
- `08_structural_covariance_voxelwise_mean.py` - script to compute qT1 structural covariance
- `fig3/on_cortex.ipynb` - cortex projections
- `fig3/Yeo_violin.ipynb` - Decoding with Yeo Networks

### Supp Fig 1.  Robustness of Structural Connectivity Gradients.
- `11_Supp_compute_gradients_with_diff_threshhold.py` - Supplementary, computes TC structural connectivity gradients with matrices tresholded at different percentiles
- `fig1/thalamic_gradients-diff_thresholds.ipynb` - plot gradients


### Supp Fig 2.  Additional Thalamocortical Structural Connectivity Gradients.
- `12_Supp_gradient_3_4.py` - save gradient 3 and 4 as nifti
- `Supp_gradient3_4/plot_gradient_3_4.ipynb` - plot gradients
  
### Supp Fig 3.  Thalamocortical Structural Connectivity Gradients - Right Hemisphere (RH).
- `fig1/01_thalamic_gradients.ipynb` - plotting gradients
- `fig1/02_Decoding_with_Thomas_atlas_right.ipynb` - Decoding with THOMAS

### Supp Fig 4. Thalamocortical Functional Connectivity Gradients - Right Hemisphere (RH). 
- `fig2/func_gradients_right.ipynb` - plotting right hem functional gradient

### Supp Fig 5. Correlation between the Individual- and Corresponding Group-Level Maps
- `13_Supp_individual_gradients_fc.py` - calculate fc gradients at subject level and correlate with group-level
- `14_Supp_individual_gradients_sc.py` - calculate sc gradients at subject level and correlate with group-level
- `15_Supp_individual_T1q_maps.py` - get individual qT1 maps and correlate with group-level
- `Supp_individual_diff/individual_grouplevel_corr.ipynb` -  create violinplot

### Supp Fig 6.  Cortical Projections of Structural Connectivity Gradients and their Association to Functional Connectivity and Structural Covariance - Right Hemisphere (RH)
- `fig3/on_cortex-right_hem.ipynb` - plotting cortex projections right
- `fig3/Yeo_violin-right_hem.ipynb` - Decoding with Yeo Networks right

### Supp Fig 7.   Cross-Check of Structural Covariance Results. 
- `10_Supp_QC_scov.py` - Supplementary, crosscorrelation: correlating between left voxelwise thalamic T1q and wholebrain cortex parcel T1q ->links this to left structural connectivity gradients by correlating columns (parcels) of structural covariance with G1 and G2 gradient loadings (same for right)
- `Supp_QC_scov/QC_scov.ipynb` - plotting

### Supp Fig 8. Projections Based on Thomas Nuclei  
- `09_Supp_groundtruth.py` - Supplementary, mean projection of sc, fc, and scov of thomas nucleus (AV, VLP, MD) on surface
- `Supp_thomas_nuclei_projection/sc_thomas nuclei.ipynb` - plotting

### Supp Fig 9. Thalamic SNR and tSNR
- `Supp_snr/thalamic_snr.ipynb` - plotting


## **Dependencies**

`conda_env` contains the environment.yml files to recreate the conda environments


Key dependencies are 
- [Brainspace v. 0.1.3](https://brainspace.readthedocs.io/en/latest/index.html)
- [Brainsmash v. 0.11.0](https://brainsmash.readthedocs.io/en/latest/)







