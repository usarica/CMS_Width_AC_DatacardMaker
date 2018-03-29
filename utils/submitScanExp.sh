#!/bin/bash

let INCREMENT=5
let COUNTER=0
let CTP=0
let minvar=0
let maxVar=0

fname=$1
wname="../../"$2
poi=$3
let npoints=$4

scr="scanExp.slurm.sh"

mkdir -p $fname"/Logs"
cp $scr $fname"/"
pushd $fname

while [  $maxVar -lt $npoints ];
do
	let CTP=$COUNTER+1
	let minVar=$COUNTER*$INCREMENT
	let maxVar=$CTP*$INCREMENT
	if [ $maxVar -gt $npoints ]; then
		let maxVar=$npoints
	fi
	sbatch --output="./Logs/lsflog_ScanExp_"$fname"_"$minVar"_"$maxVar".txt" --error="./Logs/lsferr_ScanExp_"$fname"_"$minVar"_"$maxVar".err" $scr $wname $fname $poi $npoints $minVar $maxVar
	let COUNTER=COUNTER+1
done

popd
