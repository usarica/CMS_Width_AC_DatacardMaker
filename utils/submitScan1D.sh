#!/bin/bash

let INCREMENT=1
let COUNTER=0
let CTP=0
let minvar=0
let maxVar=0

fname=$1
wname="../../"$2
poi=$3
let npoints=$4
rangel=$5
rangeh=$6
let minP=$7
let maxP=$8
extarg=$9

if [ $maxP -gt $npoints ]; then
  let maxP=$npoints
fi

scr=""
cmd=""
if [[ "$hname" == *"lxplus"* ]];then
  echo "Host is on LXPLUS, so need to use LXBATCH"
  scr="scan1D.lsf.sh"
  cmd="bsub -q 2nd -C 0 -oo ./Logs/lsflog_"$fname"_"$minVar"_"$maxVar".txt -eo ./Logs/lsferr_"$fname"_"$minVar"_"$maxVar".err"
elif [[ "$hname" == *"login-node"* ]] || [[ "$hname" == *"bc-login"* ]]; then
  echo "Host is on MARCC, so need to use SLURM batch"
  scr="scan1D.slurm.sh"
  cmd="sbatch --output=./Logs/lsflog_"$fname"_"$minVar"_"$maxVar".txt --error=./Logs/lsferr_"$fname"_"$minVar"_"$maxVar".err"
fi


mkdir -p $fname"/Logs"
cp $scr $fname"/"
pushd $fname

while [  $maxVar -lt $maxP ];
do
  let CTP=$COUNTER+1
  let minVar=$COUNTER*$INCREMENT
  let maxVar=$CTP*$INCREMENT
  if [ $maxVar -gt $maxP ]; then
    let maxVar=$maxP
  fi
  $cmd $scr $wname $fname $poi $rangel $rangeh $npoints $minVar $maxVar "$extarg"
  let COUNTER=COUNTER+1
done

popd
