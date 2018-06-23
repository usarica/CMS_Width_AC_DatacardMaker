#!/bin/bash

fname=$1
wname=$2
sqrts=$3
opt=$4

scr="buildCards.slurm.sh"

mkdir -p Logs
sbatch --output="./Logs/lsflog_buildlog_"$wname".txt" --error="./Logs/lsferr_buildlog_"$wname".err" $scr $fname $wname $sqrts $opt
