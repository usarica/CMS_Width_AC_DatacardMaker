#!/bin/bash

wname=$1
outname=$2
poi=$3
rangel=$4
rangeh=$5
let npoints=$6
let firstpoint=$7
let lastpoint=$8
if [ $lastpoint -lt $npoints ]; then
  let lastpoint=$lastpoint-1
fi
if [ $lastpoint -lt $firstpoint ]; then
  let lastpoint=$firstpoint
fi
extarg=$9

cmdadd=""
if [ $firstpoint -ge 0 ];then
  cmdadd=$cmdadd" --firstPoint "$firstpoint
  outname=$outname"_"$firstpoint
fi

if [ $lastpoint -ge 0 ];then
  cmdadd=$cmdadd" --lastPoint "$lastpoint
  outname=$outname"_"$lastpoint
fi

if [[ "$poi" == "GGsm" ]];then
  echo "POI is GGsm, fixing fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "GGsm_floatfai1" ]];then
  echo "POIs are GGsm, floating fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "GGsm_floatfai1Narrow" ]];then
  echo "POIs are GGsm, floating fai1 but with narrow range"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh":CMS_zz4l_fai1=-0.005:0.005"
elif [[ "$poi" == "GGsm_fixMH" ]];then
  echo "POI is GGsm, but fixing MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF,MH --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "MH" ]];then
  echo "POI is MH, fixing other POIs"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=MH --freezeParameters=CMS_zz4l_fai1,GGsm,kbkg_VBF --setParameterRanges MH="$rangel","$rangeh
elif [[ "$poi" == "GGsm_floatMH" ]];then
  echo "POI is GGsm, but floating MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "fai1" ]];then
  echo "POI is fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=GGsm,kbkg_VBF --setParameterRanges CMS_zz4l_fai1="$rangel","$rangeh
elif [[ "$poi" == "GGsm_fai1" ]];then
  echo "POIs are GGsm and fai1"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=GGsm,CMS_zz4l_fai1 --freezeParameters=kbkg_VBF  --setParameterRanges GGsm="$rangel":CMS_zz4l_fai1="$rangeh
elif [[ "$poi" == "GGsm_MH" ]];then
  echo "POIs are GGsm and MH"
  cmdadd=$cmdadd" -m 125 --redefineSignalPOIs=MH,GGsm --freezeParameters=kbkg_VBF,CMS_zz4l_fai1  --setParameterRanges GGsm="$rangel":MH="$rangeh
fi
if [[ "$extarg" == *"fastScan"* ]];then
  cmdadd=$cmdadd" --fastScan"
fi

cmd="-M MultiDimFit "$wname" --algo=grid --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd MINIMIZER_no_analytic --alignEdges=1 --saveNLL --saveSpecifiedNuis=all --saveSpecifiedFunc=R,RV,RF,R_13TeV,RV_13TeV,RF_13TeV -v 3 -t -1 --points "$npoints" -n "$outname" "$cmdadd
if [[ "$extarg" == *"noSyst"* ]];then
  cmd=${cmd/"--freezeParameters="/"--freezeParameters=allConstrainedNuisances,"}
fi
if [[ "$extarg" == *"obs"* ]];then
  cmd=${cmd/"-t -1 "/" "}
fi

# Add protection for fluctuations
if [[ "$extarg" == *"ALTFIT"* ]];then
  cmd=$cmd" --startFromPreFit 1 --cminPreScan --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerPrecision 0.000001"
fi

echo "Command: combine "$cmd
combine $cmd
