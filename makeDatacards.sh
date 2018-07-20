#!/bin/bash

todaysdate=$1
massregion=$2
hypo=$3
data=$4

compoundname=$todaysdate"_"$massregion"_"$data"_"$hypo

scr="makeDatacards.slurm.sh"

mkdir -p Logs
sbatch --output="./Logs/lsflog_"$compoundname".txt" --error="./Logs/lsferr_"$compoundname".err" $scr $todaysdate $massregion $hypo $data
