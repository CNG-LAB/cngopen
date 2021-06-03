
# Changing the social brain
> Data and code to use the structural manifold as a standard space or create your own  <br />
> From "Changing the social brain: plasticity along macro-scale axes of functional connectivity following social mental training", Valk et al., 2021  <br />
> https://www.biorxiv.org/content/10.1101/2020.11.11.377895v2 <br />

---

## Table of Contents

- [Dependencies](#step-by-step)
- [Step by step guide to analysis approach](#step-by-step)
- [Manuscript data](#manuscript-data)
- [Support](#support)

---

### Dependencies

- For the current work, we rely on brainspace (download here: https://brainspace.readthedocs.io/en/latest/pages/install.html), surfstat: https://mica-mni.github.io/surfstat/, scientific colourmaps: https://www.fabiocrameri.ch/colourmaps/, raincloudplots: https://github.com/RainCloudPlots/RainCloudPlots and canlab-core: https://github.com/canlab/CanlabCore for prediction models.
- For the analyses we use surfaces constructed from freesurfer (fsaverage 5) and in case of parcellations, we used the Schaefer 400 (7 networks) parcel solution.

### Step by step

Compute your own functional eccentricity change the following steps. <br />
This tutorial can be in theory followed using any longitudinal dataset with the following structure: [connectome x participant x timepoint]

1. First we performed gradient analysis and aligned the individual gradients to the template of the human connectome young adult sample (Glasser et al, 2013).

2. Eccentricity; sum of squared gradient 1 - 3.

```matlab
  % vertex level eccentricity
  for i = 1:992 %number of datapoints
    GGG(i,:) = sqrt((G1(i,:).^2)+(G2(i,:).^2)+(G3(i,:).^2));
  end
```
3. Compute the difference scores between timepoints - see the matlab script: load_data for the first three gradients as well as the eccentricity measure.


### Manuscript data

- Data corresponding to the Figures presented in the manuscript can be found in the F1.mat, F2.mat and F3.mat files.
- These files contain the linear models associated with the main analysis comparing the three modules with each other (F1 and F2) as well as the outcome of the prediction analysis (F3)

---

### Support

Feel free to get in touch if you have any questions (valk@cbs.mpg.de)
