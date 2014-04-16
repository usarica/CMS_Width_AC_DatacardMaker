#!/bin/bash

First=$1
Last=$2

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd` with options $First $Last

eval `scram runtime -sh`

combine -M MultiDimFit hzz4l_all.root --algo=grid --points 200 -m 125.6 -n Exp_nLL_scan_$3_${First} -t -1 --expectSignal=1 --firstPoint $First --lastPoint $Last -v 1 --setPhysicsModelParameters CMS_zz4l_GGsm=$3 --setPhysicsModelParameterRanges CMS_zz4l_GGsm=0.000001,$4