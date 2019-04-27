#!/bin/bash


todaysdate=$1
option=$2
curdir=$(pwd)

if [[ "$option" == *"offshell"* ]];then
echo "Doing off-shell cards"
# Make off-shell datacards
for p in SM; do
  for data in 7TeV_2011 8TeV_2012; do
    . makeDatacards.sh $todaysdate "Offshell" $p $data
  done
done
fi

if [[ "$option" == *"onshell"* ]];then
echo "Doing on-shell cards"
# Make on-shell datacards
for p in SM L1 L1ZGs a2 a3; do
  for data in 7TeV_2011 8TeV_2012 13TeV_2015; do
    . makeDatacards.sh $todaysdate "Onshell" $p $data
  done
done
fi
