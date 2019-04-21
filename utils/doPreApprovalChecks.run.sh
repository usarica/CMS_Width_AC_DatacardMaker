#!/bin/bash

wname=$1
outname=$2
poi=$3
rangel=$4
rangeh=$5
extarg=$6

cmdadd=""
cmdadd_bkgonly=""

if [[ "$poi" == "GGsm" ]];then
  echo "POI is GGsm"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
  cmdadd_bkgonly=$cmdadd_bkgonly" -m 125 --redefineSignalPOIs=GGsm --setParameters RF=0,RV=0 --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "GGsm_fixMH" ]];then
  echo "POI is GGsm, but fixing MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF,MH --setParameterRanges GGsm="$rangel","$rangeh
  cmdadd_bkgonly=$cmdadd_bkgonly" -m 125 --redefineSignalPOIs=GGsm --setParameters RF=0,RV=0 --freezeParameters=CMS_zz4l_fai1,kbkg_VBF,MH --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "GGsm_floatMH" ]];then
  echo "POI is GGsm, but floating MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
  cmdadd_bkgonly=$cmdadd_bkgonly" -m 125 --redefineSignalPOIs=GGsm --setParameters RF=0,RV=0 --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "fai1" ]];then
  echo "POI is fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=GGsm,kbkg_VBF --setParameterRanges CMS_zz4l_fai1="$rangel","$rangeh
  cmdadd_bkgonly=$cmdadd_bkgonly" -m 125 --redefineSignalPOIs=CMS_zz4l_fai1 --setParameters RF=0,RV=0 --freezeParameters=GGsm,kbkg_VBF --setParameterRanges CMS_zz4l_fai1="$rangel","$rangeh
elif [[ "$poi" == "GGsm_fai1" ]];then
  echo "POIs are GGsm and fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm,CMS_zz4l_fai1 --freezeParameters=kbkg_VBF  --setParameterRanges GGsm="$rangel":CMS_zz4l_fai1="$rangeh
  cmdadd_bkgonly=$cmdadd_bkgonly" -m 125 --redefineSignalPOIs=GGsm,CMS_zz4l_fai1 --setParameters RF=0,RV=0 --freezeParameters=kbkg_VBF  --setParameterRanges GGsm="$rangel":CMS_zz4l_fai1="$rangeh
fi

if [[ "$extarg" == *"obs"* ]];then
  cmdadd=${cmdadd/"-t -1 "/" "}
  cmdadd_bkgonly=${cmdadd_bkgonly/"-t -1 "/" "}
fi

# Add protection for fluctuations
if [[ "$extarg" == *"ALTFIT"* ]];then
  cmdadd=$cmdadd" --startFromPreFit 1 --cminPreScan --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerPrecision 0.000001"
  cmdadd_bkgonly=$cmdadd_bkgonly" --startFromPreFit 1 --cminPreScan --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerPrecision 0.000001"
fi


cmd="-M FitDiagnostics "$wname" -n fitdiagnositcs_bkgonly --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 -v 3 -t -1 "$cmdadd_bkgonly $extarg
echo "Command: combine "$cmd
combine $cmd

cmd="-M FitDiagnostics "$wname" -n fitdiagnositcs_sigbkg --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 -v 3 -t -1 "$cmdadd $extarg
echo "Command: combine "$cmd
combine $cmd
