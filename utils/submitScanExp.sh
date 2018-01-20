#!/bin/bash
NMIN=0
NMAX=40
INCREMENT=5
COUNTER=$NMIN
CTP=0
fname=$1
wname="../../"$2
poi=$3

scr="scanExp.slurm.sh"

mkdir -p $fname
cp $scr $fname"/"
pushd $fname

mkdir -p Logs
while [  $COUNTER -lt $NMAX ];
do
	let CTP=$COUNTER+1
	let minVar=$COUNTER*$INCREMENT
	let maxVar=$CTP*$INCREMENT
	sbatch --output="./Logs/lsflog_ScanExp_"$fname"_"$minVar"_"$maxVar".txt" --error="./Logs/lsferr_ScanExp_"$fname"_"$minVar"_"$maxVar".err" $scr $wname $fname $minVar $maxVar $poi
	let COUNTER=COUNTER+1
done
popd
