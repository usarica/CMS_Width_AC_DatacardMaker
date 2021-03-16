#!/bin/bash

date=$1
hypo=$2
year=$3
tpltag=$4
dcinputstag=$5
hasBoostedVH=$6

theCatScheme="nj012cat"
theCatDir="CatScheme_Nj"
theCats=( Nj_eq_0  Nj_eq_1 Nj_geq_2 )
if [[ $hasBoostedVH -eq 1 ]]; then
  theCatScheme="nj012boostedhadvhcat"
  theCatDir="CatScheme_Nj_BoostedHadVH"
  theCats+=( BoostedHadVH)
fi
inputcards=/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DatacardSpecs/${dcinputstag}/Offshell_inputs_13TeV_${year}/${theCatDir}/${hypo}
tpldir=/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/Templates/${tpltag}/${theCatDir}/${hypo}/${year}
outdirname=${date}_ZZ2L2Nu_InputTag_${dcinputstag}_TplTag_${tpltag}_${hypo}_${theCatScheme}

if [[ "$(hostname)" == *"lxplus"* ]]; then
  inputcards=${inputcards//'/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/'}
  tpldir=${tpldir//'/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/'}
fi

for channel in 2e2nu 2mu2nu; do
  echo "--------------"
  echo "Channel: ${channel}"
  echo "--------------"
  for cat in "${theCats[@]}"; do
    echo "Category: $cat"
    theCoords="mass:KD1"
    if [[ "$hypo" != "SM" ]] && [[ "$cat" == "Nj_geq_2" ]]; then
      theCoords="mass:KD1:KD2"
    elif [[ "$cat" == "BoostedHadVH" ]]; then
      theCoords="mass"
    fi
    echo "Coordinates = ${theCoords}"

    if [[ "$hypo" == "SM" ]]; then
      python makeWidthDatacards.py -b --writeoutput -i $inputcards -t $tpldir -a $outdirname --coord $theCoords --GHmodel 1 --CatScheme $theCatScheme --ac 0 --mLow 200 --mHigh 13000 --category $cat --channel $channel
    else
      python makeWidthDatacards.py -b --writeoutput -i $inputcards -t $tpldir -a $outdirname --coord $theCoords --GHmodel 1 --CatScheme $theCatScheme --ac 2 --mLow 200 --mHigh 13000 --category $cat --channel $channel
    fi

    RUN_STATUS=$?
    if [[ ${RUN_STATUS} -ne 0 ]]; then
      echo "Datacard compilation crashed  with error code ${RUN_STATUS}. Exiting..."
      exit ${RUN_STATUS}
    fi

    echo ""
  done
done
