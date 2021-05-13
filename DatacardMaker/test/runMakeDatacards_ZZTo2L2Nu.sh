#!/bin/bash

date=$1
tpltag=$2
dcinputstag=$3


for hasBoostedVH in 0; do
  for hypo in SM a2 a3 L1; do

    for year in 2016 2017 2018; do
      ./makeDatacards_ZZTo2L2Nu.run.sh $date $hypo $year $tpltag $dcinputstag $hasBoostedVH
    done

    theCatScheme="nj012_2l2nu"
    if [[ $hasBoostedVH -eq 1 ]]; then
      theCatScheme="nj012_boostedhadvh_2l2nu"
    fi
    dirname=cards_${date}_ZZ2L2Nu_InputTag_${dcinputstag}_TplTag_${tpltag}_${hypo}_${theCatScheme}/HCG

    mkdir -p ${dirname}/13TeV
    for f in loadLib.C buildCards.run.sh; do
      cp utils/$f ${dirname}/13TeV/
    done

    cd ${dirname}/13TeV
    combineCards.py ../13TeV_201?/*.txt > hto2l2nu_allcats_13TeV.txt
    ./buildCards.run.sh hto2l2nu_allcats_13TeV.txt hto2l2nu_allcats_13TeV.root 13 $hypo
    for year in 2016 2017 2018; do
      combineCards.py ../13TeV_${year}/*.txt > hto2l2nu_allcats_13TeV_${year}.txt
      ./buildCards.run.sh hto2l2nu_allcats_13TeV_${year}.txt hto2l2nu_allcats_13TeV_${year}.root 13 $hypo
    done
    cd -

  done
done
