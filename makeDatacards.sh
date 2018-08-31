#!/bin/bash

todaysdate=$1
massregion=$2
hypo=$3
data=$4

mkdir -p Logs
compoundname=$todaysdate"_"$massregion"_"$data"_"$hypo

scr=""
cmd=""
hname=$(hostname)
if [[ "$hname" == *"lxplus"* ]];then
  echo "Host is on LXPLUS, so need to use LXBATCH"
  scr="makeDatacards.lsf.sh"
  cmd="bsub -q 8nh -C 0 -o ./Logs/lsflog_"$compoundname".txt -e ./Logs/lsferr_"$compoundname".err"
elif [[ "$hname" == *"login-node"* ]] || [[ "$hname" == *"bc-login"* ]]; then
  echo "Host is on MARCC, so need to use SLURM batch"
  scr="makeDatacards.slurm.sh"
  cmd="sbatch --output=./Logs/lsflog_"$compoundname".txt --error=./Logs/lsferr_"$compoundname".err"
fi

$cmd $scr $todaysdate $massregion $hypo $data
