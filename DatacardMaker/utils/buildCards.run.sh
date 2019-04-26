#!/bin/bash

inname=$1
outname=$2
sqrts=$3
opt=$4

echo "Building $outname from $inname at sqrts=$sqrts and option $opt"

if [[ "$opt" == "SM" ]];then
   echo "SM case:"
   echo "Command:" text2workspace.py \
      -m 125 $inname -o $outname \
      -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
      --PO="offshell" --PO="fai1fixed" --PO="sqrts="$sqrts
   text2workspace.py \
      -m 125 $inname -o $outname \
      -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
      --PO="offshell" --PO="fai1fixed" --PO="sqrts="$sqrts
else
   echo "BSM case:"
   echo "Command:" text2workspace.py \
      -m 125 $inname -o $outname \
      -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
      --PO="offshell" --PO="allowPMF" --PO="sqrts="$sqrts
   text2workspace.py \
      -m 125 $inname -o $outname \
      -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
      --PO="offshell" --PO="allowPMF" --PO="sqrts="$sqrts
fi

