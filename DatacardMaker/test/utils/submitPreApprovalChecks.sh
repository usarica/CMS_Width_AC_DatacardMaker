#!/bin/bash

fname=$1
wname="../../"$2
poi=$3
rangel=$4
rangeh=$5
extarg=$6

scr=""
cmd=""
hname=$(hostname)
if [[ "$hname" == *"lxplus"* ]];then
  echo "Host is on LXPLUS, so need to use LXBATCH"
  scr="doPreApprovalChecks.lsf.sh"
  cmd="bsub -q 2nd -C 0 -oo ./Logs/lsflog_"$fname".txt -eo ./Logs/lsferr_"$fname".err"
elif [[ "$hname" == *"login-node"* ]] || [[ "$hname" == *"bc-login"* ]]; then
  echo "Host is on MARCC, so need to use SLURM batch"
  scr="doPreApprovalChecks.slurm.sh"
  cmd="sbatch --output=./Logs/lsflog_"$fname".txt --error=./Logs/lsferr_"$fname".err"
fi

mkdir -p $fname"/Logs"
scrcore=${scr%%"."*}
cp "$scrcore"* $fname"/"
pushd $fname

$cmd $scr $wname $fname $poi $rangel $rangeh $extarg

popd
