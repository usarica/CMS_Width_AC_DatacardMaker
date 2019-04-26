#!/bin/bash

fname=$1
wname=$2
sqrts=$3
opt=$4

mkdir -p Logs

scr=""
cmd=""
hname=$(hostname)
if [[ "$hname" == *"lxplus"* ]];then
  echo "Host is on LXPLUS, so need to use LXBATCH"
  scr="buildCards.lsf.sh"
  cmd="bsub -q 2nd -C 0 -oo ./Logs/lsflog_buildlog_"$wname".txt -eo ./Logs/lsferr_buildlog_"$wname".err"
elif [[ "$hname" == *"login-node"* ]] || [[ "$hname" == *"bc-login"* ]]; then
  echo "Host is on MARCC, so need to use SLURM batch"
  scr="buildCards.slurm.sh"
  cmd="sbatch --output=./Logs/lsflog_buildlog_"$wname".txt --error=./Logs/lsferr_buildlog_"$wname".err"
else
  echo "Could not identify type of host $hname"
fi

$cmd $scr $fname $wname $sqrts $opt
