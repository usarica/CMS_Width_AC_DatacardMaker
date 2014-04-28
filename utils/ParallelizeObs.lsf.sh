#!/bin/bash

First=$1
Last=$2
Kind=$3
Name=$3
Range=30
Points=200
if [ $Kind == "Combined" ]
then 
Kind=""
fi
if  [ "$4" != "" ]
then
Range=$4
Name="$3lowR"
Points=500
fi

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd` with options $First $Last

eval `scram runtime -sh`

combine -M MultiDimFit hzz${Kind}_all.root --algo=grid --points $Points -m 125.6 -n Obs_${Name}_nLL_scan_${First} --firstPoint $First --lastPoint $Last -v 1 --setPhysicsModelParameterRanges CMS_zz4l_GGsm=0.000001,$Range 