#!/bin/bash

echo "BEGIN host info"
hostname
uptime
free
ulimit -a
echo "END host info"
set -x

SUB=$(echo "$1")
echo $SUB

mkdir /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"
mkdir /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01
mkdir /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi
mkdir /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx

mrconvert /p_02666/mica-mics/singularity_out/SUB_$SUB/micapipe/sub-$SUB/ses-01/dwi/sub-${SUB}_ses-01_space-dwi_desc-dwi_preproc.mif /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/sub-"$SUB"_ses-01_space-dwi_desc-dwi_preproc.nii.gz -export_grad_fsl /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/bvecs /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/bvals

mv /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/sub-"$SUB"_ses-01_space-dwi_desc-dwi_preproc.nii.gz /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/data.nii.gz

cp /p_02666/mica-mics/singularity_out/SUB_"$SUB"/micapipe/sub-"$SUB"/ses-01/dwi/sub-"$SUB"_ses-01_space-dwi_desc-brain_mask.nii.gz /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/nodif_brain_mask.nii.gz

