#!/bin/bash

version=$1
pushd cards_"$version"_Combination_a3/HCG/Scans/
#. submitScanExp.sh GGsmFloated_fai1Fixed_Old2016_13TeV hzz4l_Prop_All_13TeV_old.root GGsm 51 0 10 0 51
#. submitScanExp.sh fai1Floated_GGsmFixed_Old2016_13TeV hzz4l_Prop_All_13TeV_old.root fai1 31 -0.02 0.02 0 31
#. submitScanExp.sh fai1Floated_GGsmFixed_Old2016_13TeV_extra1 hzz4l_Prop_All_13TeV_old.root fai1 21 -0.15 0.15 0 21
popd

pushd cards_"$version"_Combination_a2/HCG/Scans/
#. submitScanExp.sh GGsmFloated_fai1Fixed_Old2016_13TeV hzz4l_Prop_All_13TeV_old.root GGsm 51 0 10 0 51
. submitScanExp.sh fai1Floated_GGsmFixed_Old2016_13TeV hzz4l_Prop_All_13TeV_old.root fai1 31 -0.015 0.015 0 31
. submitScanExp.sh fai1Floated_GGsmFixed_Old2016_13TeV_extra1 hzz4l_Prop_All_13TeV_old.root fai1 31 -0.15 0.25 0 31
popd

pushd cards_"$version"_Combination_L1/HCG/Scans/
#. submitScanExp.sh GGsmFloated_fai1Fixed_Old2016_13TeV hzz4l_Prop_All_13TeV_old.root GGsm 51 0 10 0 51
. submitScanExp.sh fai1Floated_GGsmFixed_Old2016_13TeV hzz4l_Prop_All_13TeV_old.root fai1 31 -0.015 0.015 0 31
. submitScanExp.sh fai1Floated_GGsmFixed_Old2016_13TeV_extra1 hzz4l_Prop_All_13TeV_old.root fai1 41 -0.7 0.3 0 41
popd

