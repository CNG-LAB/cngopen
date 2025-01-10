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

/usr/bin/time -v singularity run --no-home --cleanenv  -B /p_02666/mica-mics/rawdata -B /p_02666/mica-mics/singularity_out/SUB_$SUB -B /hu_ajohn -B /afs/cbs/software/freesurfer/licensekeys:/opt/freesurfer-6.0.0/license.txt      /t_SoftwareServiceLinux_sc/micapipe/0.1.2/1  -bids /p_02666/mica-mics/rawdata -out /p_02666/mica-mics/singularity_out/SUB_$SUB  -sub $SUB -ses 01  -keep_tck -proc_structural -proc_dwi  -nocleanup
