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

. ./env_fsl
unset FSLPARALLEL
env

mkdir /p_02666/mica-mics/probtrackx_out/sub_"$SUB"
mkdir /p_02666/mica-mics/probtrackx_out/sub_"$SUB"/ses-01
mkdir /p_02666/mica-mics/probtrackx_out/sub_"$SUB"/ses-01/dwi
mkdir /p_02666/mica-mics/probtrackx_out/sub_"$SUB"/ses-01/dwi/left_hem
mkdir /p_02666/mica-mics/probtrackx_out/sub_"$SUB"/ses-01/dwi/right_hem

#left hemisphere
probtrackx2 -x /p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_lh.nii.gz -s /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/merged -m /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/nodif_brain_mask --xfm=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2diff_warp.nii.gz --invxfm=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2standard_warp.nii.gz --targetmasks=/p_02666/mica-mics/schaefer200_space_mni/targets.txt -l --pd -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --omatrix2 --target2=/fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_2mm_brain_mask.nii.gz --forcedir --opd --dir=/p_02666/mica-mics/probtrackx_out/sub_"$SUB"/ses-01/dwi/left_hem/ --os2t

#right hemisphere
probtrackx2 -x /p_02666/mica-mics/thalamus_space_mni/space-MNI125_atlas-thalamus_rh.nii.gz -s /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/merged -m /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/nodif_brain_mask --xfm=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2diff_warp.nii.gz --invxfm=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2standard_warp.nii.gz --targetmasks=/p_02666/mica-mics/schaefer200_space_mni/targets.txt -l --pd -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --omatrix2 --target2=/fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_2mm_brain_mask.nii.gz --forcedir --opd --dir=/p_02666/mica-mics/probtrackx_out/sub_"$SUB"/ses-01/dwi/right_hem/ --os2t


