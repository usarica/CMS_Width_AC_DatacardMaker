#!/bin/bash

date=$1
hypo=$2
year=$3
tpltag=$4
dcinputstag=$5
datainputtag=$6
hasBoostedVH=$7

getBestDirectory(){
  local dirmain=$1
  local chkdirs=(
  $(pwd) \
    /home/users/usarica/work/Width_AC_Run2/2L2Nu_MiniAOD/production/slc6_photonextra/CMSSW_10_2_22/src/CMS3/AnalysisTree/test \
    /hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker \
  )
  for chkdir in "${chkdirs[@]}"; do
    if [[ -e ${chkdir}/${dirmain} ]]; then
      echo ${chkdir}/${dirmain}
      break
    fi
  done
}

theCatScheme="nj012_3l1nu"
theCatDir="CatScheme_Nj"
theCats=( Nj_eq_0 Nj_eq_1 Nj_geq_2 )
if [[ $hasBoostedVH -eq 1 ]]; then
  echo "Boosted VH is not implemented in ZW final states."
  exit 1
fi

inputcards="$(getBestDirectory output/DatacardSpecs/ZWTo3L1Nu/${dcinputstag}/Offshell_inputs_13TeV_${year}/${theCatDir}/${hypo})"
tpldir="$(getBestDirectory output/Templates/${tpltag}/${theCatDir}/${hypo}/${year})"
#datadir="$(getBestDirectory output/DCDataTrees_ZWTo3L1Nu/${datainputtag}/${theCatDir}/${hypo})"
datadir="$(getBestDirectory output/DCDataTrees_ZWTo3L1Nu/${datainputtag}/${theCatDir})"

echo "Inputs: $inputcards"
echo "Templates: $tpldir"
echo "Data: $datadir"

outdirname=${date}_ZW3L1Nu_InputTag_${dcinputstag}_TplTag_${tpltag}_${hypo}_${theCatScheme}
mLowVal=150

for channel in 2l1e 2l1mu; do
  echo "--------------"
  echo "Channel: ${channel}"
  echo "--------------"
  for cat in "${theCats[@]}"; do
    echo "Category: $cat"
    theCoords="mass"
    echo "Coordinates = ${theCoords}"

    if [[ "$hypo" == "SM" ]]; then
      python makeWidthDatacards.py -b --writeoutput -i $inputcards -t $tpldir -r $datadir -a $outdirname --coord $theCoords --GHmodel 1 --CatScheme $theCatScheme --ac 0 --mLow ${mLowVal} --mHigh 13000 --category $cat --channel $channel
    else
      python makeWidthDatacards.py -b --writeoutput -i $inputcards -t $tpldir -r $datadir -a $outdirname --coord $theCoords --GHmodel 1 --CatScheme $theCatScheme --ac 2 --mLow ${mLowVal} --mHigh 13000 --category $cat --channel $channel
    fi

    RUN_STATUS=$?
    if [[ ${RUN_STATUS} -ne 0 ]]; then
      echo "Datacard compilation crashed with error code ${RUN_STATUS}. Exiting..."
      exit ${RUN_STATUS}
    fi

    echo ""
  done
done
