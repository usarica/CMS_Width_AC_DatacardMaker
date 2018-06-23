#!/bin/bash

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

scr="doPreApprovalChecks.slurm.sh"

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
	sbatch --output="./Logs/lsflog_PAC_"$fname"_"$minVar"_"$maxVar".txt" --error="./Logs/lsferr_PAC_"$fname"_"$minVar"_"$maxVar".err" $scr $wname $fname $poi $rangel $rangeh $npoints $minVar $maxVar $extarg
	let COUNTER=COUNTER+1
done

popd
