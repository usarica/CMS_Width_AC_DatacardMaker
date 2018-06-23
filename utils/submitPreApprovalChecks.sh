#!/bin/bash

fname=$1
wname="../../"$2
poi=$3
rangel=$4
rangeh=$5
extarg=$6

scr="doPreApprovalChecks.slurm.sh"

mkdir -p $fname"/Logs"
cp $scr $fname"/"
pushd $fname

sbatch --output="./Logs/lsflog_PAC_"$fname".txt" --error="./Logs/lsferr_PAC_"$fname".err" $scr $wname $fname $poi $rangel $rangeh $extarg

popd
