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


/usr/bin/time -v bedpostx /p_02666/mica-mics/bedpostx_out/SUB_"$SUB"/ses-01/dwi/bedpostx --nf=3 --fudge=1 --bi=1000 --model=2


