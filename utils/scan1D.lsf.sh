#!/bin/bash

RUNDIR=${LS_SUBCWD}
cd $RUNDIR
echo "LSF job running in: "$(pwd)

eval `scram runtime -sh`

echo $CMSSW_VERSION

source scan1D.run.sh $@
