#!/bin/bash

inname=$1
outname=$2
sqrts=$3
opt=$4

echo "Building $outname from $inname at sqrts=$sqrts and option $opt"


runCmd="text2workspace.py -m 125 $inname -o $outname -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs --PO=sqrts=${sqrts}"
if [[ "$opt" == *"SM"* ]]; then
  echo "Configuring for the SM case..."
  runCmd="${runCmd} --PO=fai1fixed"
else
  echo "Configuring for the BSM case..."
  runCmd="${runCmd} --PO=allowPMF"
fi


if [[ "$opt" == *"UseRVoverRF"* ]]; then
  echo "Using (R, RV) instead of (RF, RV) for the parametrization..."
  runCmd="${runCmd} --PO=uservoverrf"
fi

if [[ "$opt" == *"DecoupleOnshellMu"* ]]; then
  echo "Decouplings on-shell mu parameters..."
  runCmd="${runCmd} --PO=splitonshellmuvf"
  if [[ "$opt" == *"Uservoverrfonshell"* ]]; then
    echo "Using (r, rv) instead of (rf, rv) for the parametrization..."
    runCmd="${runCmd} --PO=usesplitonshellrvoverrf"
  fi
fi

if [[ "$opt" == *"DecoupleOffshellMu"* ]]; then
  echo "Decouplings off-shell mu parameters..."
  runCmd="${runCmd} --PO=splitoffshellmuvf"
  if [[ "$opt" == *"Uservoverrfoffshell"* ]]; then
    echo "Using (r_offshell, rv_offshell) instead of (rf_offshell, rv_offshell) for the parametrization..."
    runCmd="${runCmd} --PO=usesplitoffshellrvoverrf"
  fi
fi

if [[ "$opt" == *"OnshellOnly"* ]]; then
  echo "Only on-shell, so not including the off-shell option..."
else
  echo "Including the off-shell option..."
  runCmd="${runCmd} --PO=offshell"
fi

if [[ "$opt" == *"AddChannelMasks"* ]]; then
  echo "Adding channel masks..."
  runCmd="${runCmd} --channel-masks"
fi

echo "Command: ${runCmd}"
${runCmd}
