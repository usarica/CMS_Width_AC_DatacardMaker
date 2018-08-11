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

scr="scan1D.slurm.sh"

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
	sbatch --output="./Logs/lsflog_ScanExp_"$fname"_"$minVar"_"$maxVar".txt" --error="./Logs/lsferr_ScanExp_"$fname"_"$minVar"_"$maxVar".err" $scr $wname $fname $poi $rangel $rangeh $npoints $minVar $maxVar $extarg
	let COUNTER=COUNTER+1
done

popd