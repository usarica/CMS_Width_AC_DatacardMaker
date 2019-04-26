#!/bin/bash

wname=$1
outname=$2
poi=$3
rangel=$4
rangeh=$5
extarg=$6

cmdadd=""

if [[ "$poi" == "GGsm" ]];then
  echo "POI is GGsm"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "GGsm_fixMH" ]];then
  echo "POI is GGsm, but fixing MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF,MH --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "GGsm_floatMH" ]];then
  echo "POI is GGsm, but floating MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "fai1" ]];then
  echo "POI is fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=GGsm,kbkg_VBF --setParameterRanges CMS_zz4l_fai1="$rangel","$rangeh
fi

cmdcore=" -M Impacts -d "$wname" --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 -v 3 -t -1 "$cmdadd
if [[ "$extarg" == *"obs"* ]];then
  cmdcore=${cmdcore/"-t -1 "/" "}
elif [[ "$extarg" == *"exp"* ]];then
  cmdcore=${cmdcore/"-t -1 "/"-t 1 "}
fi

# Add protection for fluctuations
if [[ "$extarg" == *"ALTFIT"* ]];then
  cmdcore=$cmdcore" --startFromPreFit 1 --cminPreScan --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerPrecision 0.000001"
fi

cmd=$cmdcore" --doInitialFit "
echo "Command: combineTool.py "$cmd
combineTool.py $cmd
rename mH125.123456 mH125 *.root

cmd=$cmdcore" --doFits "
echo "Command: combineTool.py "$cmd
combineTool.py $cmd
rename mH125.123456 mH125 *.root

cmd=$cmdcore" -o impacts.json "
#cmd=${cmd/" --robustFit 1"/" "}
echo "Command: combineTool.py "$cmd
combineTool.py $cmd

plotImpacts.py -i impacts.json -o impacts
