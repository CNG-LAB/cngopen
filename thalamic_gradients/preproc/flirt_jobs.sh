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

unset FSLPARALLEL
env

#apply BET to b0 image (nodif_brain)
bet /p_02666/mica-mics/singularity_out/SUB_"$SUB"/micapipe/sub-"$SUB"/ses-01/dwi/sub-"$SUB"_ses-01_space-dwi_desc-b0 /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/sub-"$SUB"_ses-01_space-dwi_desc-b0_nodif_brain  -f 0.5 -g 0

#create transformation matrix
diff_brain=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/sub-"$SUB"_ses-01_space-dwi_desc-b0_nodif_brain
structural_brain=/p_02666/mica-mics/singularity_out/SUB_"$SUB"/micapipe/sub-"$SUB"/ses-01/anat/sub-"$SUB"_ses-01_space-nativepro_t1w_brain.nii.gz
structural=/p_02666/mica-mics/singularity_out/SUB_"$SUB"/micapipe/sub-"$SUB"/ses-01/anat/sub-"$SUB"_ses-01_space-nativepro_t1w.nii.gz
standard_brain=/fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_2mm_brain
standard=/fsl/5.0.11/ubuntu-bionic-amd64/data/standard/MNI152_T1_2mm

#linear registration
#diff2str
flirt -in "$diff_brain" -out /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/qc_flirt_diff2str -ref "$structural_brain" -omat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2str.mat -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -cost corratio

#str2diff
convert_xfm -omat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2diff.mat -inverse /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2str.mat

#str2standard
flirt -in "$structural_brain" -ref "$standard_brain" -omat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard.mat -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -cost corratio

#standard2str
convert_xfm -omat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2str.mat -inverse /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard.mat

#concatenate: diff2standard
convert_xfm -omat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2standard.mat -concat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard.mat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2str.mat

#concatenate: standard2diff
convert_xfm -omat /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2diff.mat -inverse /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2standard.mat


#nonlinear registration
#str2standard_warp
fnirt --in="$structural" --aff=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard.mat --cout=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard_warp --iout=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard_warp_iout --config=T1_2_MNI152_2mm

#standard2str_warp
invwarp -w /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard_warp -o /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2str_warp -r "$structural_brain"

#diff2standard_warp (combination of flirt and fnirt results)
convertwarp -o /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2standard_warp -r "$standard" -m /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2str.mat -w /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2standard_warp

#standard2diff_warp (combination of fnirt and flirt results)
convertwarp -o /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2diff_warp -r /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx/nodif_brain_mask -w /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/standard2str_warp --postmat=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/str2diff.mat

#apply warp for qualitycheck
applywarp --ref="$standard" --in="$diff_brain" --out=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/qc_applywarp_diff2standard --warp=/p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx.bedpostX/xfms/diff2standard_warp


