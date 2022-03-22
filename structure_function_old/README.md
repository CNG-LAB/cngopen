
# GENETIC UNCOUPLING OF STRUCTURE AND FUNCTION IN HUMAN TRANSMODAL CORTEX

> Data and code to think about structure-function associations <br />
> From "Genetic and phylogenetic uncoupling of structure and function in human transmodal cortex", Valk et al., preprint  <br />
> https://www.biorxiv.org/content/10.1101/2021.06.08.447522v1 <br />

---

## Table of Contents

- [Dependencies](#step-by-step)
- [Step by step guide to analysis approach](#step-by-step)
- [Manuscript data](#manuscript-data)
- [Support](#support)

---

### Dependencies

- For the current work, we rely on brainspace (download here: https://brainspace.readthedocs.io/en/latest/pages/install.html), surfstat: https://mica-mni.github.io/surfstat/, scientific colourmaps: https://www.fabiocrameri.ch/colourmaps/, raincloudplots: https://github.com/RainCloudPlots/RainCloudPlots; abagen https://abagen.readthedocs.io/en/stable/ and NeuroSynth https://www.neurosynth.org, MPC code: https://github.com/MICA-MNI/micaopen/tree/master/MPC

- For the analyses we use surfaces constructed from freesurfer (fsaverage 5) and Schaefer 400 (7 networks) parcel solution.

### Step by step

Check out structure and function using the following steps. <br />

1. First we correlated the rows of the structural and functional matrices.

2. We did the same for the mean functional/structural matrix with their heritability maps - as well as heritability of structure-function coupling, and for macaques.

3. To compare humans and macaque data on the surface, we used Ting Xu's code and approach: https://github.com/CNG-LAB/PRIME-DE

4. To evaluate the topology we used  brainspace to construct gradients.

5. For decoding we combined maps from step 1 and 4 in a 2D framework and mapped cytoarchitecture/functional communities, as well as
phylogenetic models (dual origin and functional reorganisation between macaques and humans), transcriptomic data (AHBA, using abagen) and
functional terms (NeuroSynth).

### Manuscript data

- Main data is from HCP and PRIME-DE which is openly available. 
- Auxiliary data is in a .mat file. Do let me know if you think something is missing.

---

### Support

Feel free to get in touch if you have any questions (valk@cbs.mpg.de)
