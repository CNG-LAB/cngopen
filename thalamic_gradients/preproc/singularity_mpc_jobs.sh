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

/usr/bin/time -v singularity run --no-home --cleanenv  -B /p_02666/mica-mics/rawdata -B /p_02666/mica-mics/singularity_out/SUB_$SUB -B /hu_ajohn -B /freesurfer/licensekeys:/opt/freesurfer-6.0.0/license.txt     /micapipe/0.1.2/1  -bids /p_02666/mica-mics/rawdata -out /p_02666/mica-mics/singularity_out/SUB_$SUB  -sub $SUB -ses 01 -proc_freesurfer -post_structural -MPC
