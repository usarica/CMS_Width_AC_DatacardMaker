#!/bin/bash

version=$1
wsname=$2
wdname=${wsname/".root"/""}
coupling=$3

if [[ "$coupling" == "SM" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_SM/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_SM/HCG/Scans/
pushd cards_"$version"_Combination_SM/HCG/Scans/
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_13TeV $wsname GGsm_fixMH 51 0 7 0 51
popd
fi

if [[ "$coupling" == "a3" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_a3/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_a3/HCG/Scans/
pushd cards_"$version"_Combination_a3/HCG/Scans/
#. submitScan1D.sh "$wdname"_GGsmfai1Floated_NoSyst_13TeV $wsname GGsm_fai1 441 0,10 -0.03,0.03 0 441 noSyst
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_13TeV $wsname GGsm 51 0 7 0 51
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_13TeV $wsname GGsm_floatfai1 51 0 7 0 51
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_13TeV $wsname fai1 31 -0.02 0.02 0 31
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_13TeV_extra1 $wsname fai1 21 -0.15 0.15 0 21
popd
fi

if [[ "$coupling" == "a2" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_a2/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_a2/HCG/Scans/
pushd cards_"$version"_Combination_a2/HCG/Scans/
#. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_NoSyst_13TeV $wsname GGsm_fai1 441 0,10 -0.03,0.04 0 441 noSyst
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_13TeV $wsname GGsm 51 0 7 0 51
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_13TeV $wsname GGsm_floatfai1 51 0 7 0 51
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_13TeV $wsname fai1 31 -0.015 0.015 0 31
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_13TeV_extra1 $wsname fai1 31 -0.15 0.25 0 31
popd
fi

if [[ "$coupling" == "L1" ]];then
#cp utils/submitScan1D.sh cards_"$version"_Combination_L1/HCG/Scans/
#cp utils/scan1D.slurm.sh cards_"$version"_Combination_L1/HCG/Scans/
pushd cards_"$version"_Combination_L1/HCG/Scans/
#. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_NoSyst_13TeV $wsname GGsm_fai1 441 0,10 -0.03,0.04 0 441 noSyst
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_13TeV $wsname GGsm 51 0 7 0 51
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_13TeV $wsname GGsm_floatfai1 51 0 7 0 51
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_13TeV $wsname fai1 31 -0.015 0.015 0 31
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_13TeV_extra1 $wsname fai1 41 -0.7 0.3 0 41
popd
fi

