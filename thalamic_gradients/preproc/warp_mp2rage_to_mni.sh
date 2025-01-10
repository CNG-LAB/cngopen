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

FSL
unset FSLPARALLEL
env

mkdir /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/

#denoise: remove noisy Hintergrund von uni_T1map
mask=/p_02666/mica-mics/rawdata/sub-$SUB/ses-01/anat/sub-"$SUB"_ses-01_acq-inv2_T1map.nii.gz  #use inv2 to make a mask
file=/p_02666/mica-mics/rawdata/sub-$SUB/ses-01/anat/sub-"$SUB"_ses-01_acq-uni_T1map.nii.gz
out=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/"$SUB"_ses-01_acq-uni_T1map_dn.nii.gz
fslmaths $mask -mul $file $out

#flirt from uni_T1 native space to MNI 1mm (just to use mat file as a prior for fnirt) [linear]
flirt -ref /fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_1mm.nii.gz -in /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/"$SUB"_ses-01_acq-uni_T1map_dn.nii.gz -out /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/flirt_uni_T1map_to_MNI1mm.nii.gz -omat /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/flirt_uni_T1map_to_MNI1mm.mat

#fnirt from uni_T1 native space to MNI 1mm [non-linear]
fnirt --ref=/fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_1mm.nii.gz --in=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/"$SUB"_ses-01_acq-uni_T1map_dn.nii.gz --aff=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/flirt_uni_T1map_to_MNI1mm.mat --cout=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_uni_T1map_to_MNI1mm_cout --iout=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_uni_T1map_to_MNI1mm_iout

#convert unimap into MNI 2mm 
flirt -ref /fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_2mm.nii.gz -in /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_uni_T1map_to_MNI1mm_iout -out /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_uni_T1map_to_MNI2mm_iout -nosearch -applyisoxfm 2

#apply transformationmatrix to mprage (T1q) image 
applywarp --ref=/fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_1mm.nii.gz --in=/p_02666/mica-mics/rawdata/sub-$SUB/ses-01/anat/sub-"$SUB"_ses-01_acq-mp2rage_T1map.nii.gz --warp=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_uni_T1map_to_MNI1mm_cout --rel --out=/p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_mp2rage_T1map_to_MNI1mm.nii.gz

#convert mp2rage into MNI 2mm 
flirt -ref /fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_2mm.nii.gz -in /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_mp2rage_T1map_to_MNI1mm -out /p_02666/mica-mics/space_mp2rage_T1map_to_MNI/SUB_$SUB/fnirt_mp2rage_T1map_to_MNI2mm -nosearch -applyisoxfm 2


