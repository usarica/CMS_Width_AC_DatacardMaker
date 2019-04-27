#!/bin/bash

version=$1
wsname=$2
wdname=${wsname/".root"/""}
coupling=$3

strSqrts=""
if [[ "$wdname" == *"7813TeV"* ]]; then
  strSqrts="7813TeV"
elif [[ "$wdname" == *"7TeV"* ]]; then
  strSqrts="7TeV"
elif [[ "$wdname" == *"8TeV"* ]]; then
  strSqrts="8TeV"
elif [[ "$wdname" == *"13TeV_2015"* ]]; then
  strSqrts="13TeV_2015"
elif [[ "$wdname" == *"13TeV_2016"* ]]; then
  strSqrts="13TeV_2016"
elif [[ "$wdname" == *"13TeV_2017"* ]]; then
  strSqrts="13TeV_2017"
elif [[ "$wdname" == *"13TeV_2018"* ]]; then
  strSqrts="13TeV_2018"
elif [[ "$wdname" == *"13TeV"* ]]; then
  strSqrts="13TeV"
fi


if [[ "$coupling" == "SM" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_SM/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_SM/HCG/Scans/
pushd cards_"$version"_Combination_SM/HCG/Scans/
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts" $wsname GGsm_fixMH 51 0 7 0 51

. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts"_ALTFIT $wsname GGsm_fixMH 51 0 7 0 51 ALTFIT
popd
fi

if [[ "$coupling" == "a3" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_a3/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_a3/HCG/Scans/
pushd cards_"$version"_Combination_a3/HCG/Scans/
#. submitScan1D.sh "$wdname"_GGsmfai1Floated_NoSyst_"$strSqrts" $wsname GGsm_fai1 441 0,10 -0.03,0.03 0 441 noSyst
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts" $wsname GGsm 51 0 7 0 51
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_"$strSqrts" $wsname GGsm_floatfai1 51 0 7 0 51
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts" $wsname fai1 31 -0.02 0.02 0 31
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_extra1 $wsname fai1 21 -0.15 0.15 0 21

. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts"_ALTFIT $wsname GGsm 51 0 7 0 51 ALTFIT
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_"$strSqrts"_ALTFIT $wsname GGsm_floatfai1 51 0 7 0 51 ALTFIT
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_ALTFIT $wsname fai1 31 -0.02 0.02 0 31 ALTFIT
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_ALTFIT_extra1 $wsname fai1 21 -0.15 0.15 0 21 ALTFIT
popd
fi

if [[ "$coupling" == "a2" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_a2/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_a2/HCG/Scans/
pushd cards_"$version"_Combination_a2/HCG/Scans/
#. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_NoSyst_"$strSqrts" $wsname GGsm_fai1 441 0,10 -0.03,0.04 0 441 noSyst
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts" $wsname GGsm 51 0 7 0 51
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_"$strSqrts" $wsname GGsm_floatfai1 51 0 7 0 51
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts" $wsname fai1 31 -0.015 0.015 0 31
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_extra1 $wsname fai1 31 -0.15 0.25 0 31

. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts"_ALTFIT $wsname GGsm 51 0 7 0 51 ALTFIT
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_"$strSqrts"_ALTFIT $wsname GGsm_floatfai1 51 0 7 0 51 ALTFIT
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_ALTFIT $wsname fai1 31 -0.015 0.015 0 31 ALTFIT
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_ALTFIT_extra1 $wsname fai1 31 -0.15 0.25 0 31 ALTFIT
popd
fi

if [[ "$coupling" == "L1" ]];then
cp utils/submitScan1D.sh cards_"$version"_Combination_L1/HCG/Scans/
cp utils/scan1D.slurm.sh cards_"$version"_Combination_L1/HCG/Scans/
pushd cards_"$version"_Combination_L1/HCG/Scans/
#. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_NoSyst_"$strSqrts" $wsname GGsm_fai1 441 0,10 -0.03,0.04 0 441 noSyst
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts" $wsname GGsm 51 0 7 0 51
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_"$strSqrts" $wsname GGsm_floatfai1 51 0 7 0 51
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts" $wsname fai1 31 -0.015 0.015 0 31
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_extra1 $wsname fai1 41 -0.7 0.3 0 41

. submitScan1D.sh "$wdname"_GGsmFloated_fai1Fixed_"$strSqrts"_ALTFIT $wsname GGsm 51 0 7 0 51 ALTFIT
. submitScan1D.sh "$wdname"_GGsmFloated_fai1Floated_"$strSqrts"_ALTFIT $wsname GGsm_floatfai1 51 0 7 0 51 ALTFIT
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_ALTFIT $wsname fai1 31 -0.015 0.015 0 31 ALTFIT
. submitScan1D.sh "$wdname"_fai1Floated_GGsmFixed_"$strSqrts"_ALTFIT_extra1 $wsname fai1 41 -0.7 0.3 0 41 ALTFIT
popd
fi

