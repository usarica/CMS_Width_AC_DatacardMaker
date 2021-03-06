#!/bin/bash

date=$1
hypo=$2
year=$3

getBestHCGDirectory(){
  local dirmain=$1
  local chkdirs=(
  $(pwd) \
    /home/users/usarica/work/Width_AC_Run2/2L2Nu_MiniAOD/datacards/CMSSW_10_2_22/src \
    /afs/cern.ch/work/u/usarica/Offshell_2l2nu/datacards/hcgdcs \
  )
  for chkdir in "${chkdirs[@]}"; do
    if [[ -e ${chkdir}/${dirmain} ]]; then
      echo ${chkdir}/${dirmain}
      break
    fi
  done
}

inputcards="$(getBestHCGDirectory hig-19-009/Decompilation/Inputs/${hypo})"
inputcards=${inputcards}/Onshell_13TeV_${year}
outdirname=${date}_ZZ4L_${hypo}_Onshell_19009
if [[ ! -d ${inputcards} ]]; then
  echo "${inputcards} does not exist."
  exit 0
fi

tpldir=${inputcards/Inputs/Templates}
#dataopt="-r dummydir"
dataopt="-r ${inputcards/Inputs/Data}"

for channel in 2e2mu 4e 4mu; do
  echo "--------------"
  echo "Channel: ${channel}"
  echo "--------------"
  theCoords="KD1"
  theCatScheme=run2legacy4l

  if [[ "$hypo" == "SM" ]]; then
    python makeWidthDatacards.py -b --writeoutput -i $inputcards -t $tpldir ${extshapearg} ${dataopt} -a $outdirname --coord $theCoords --GHmodel 0 --CatScheme $theCatScheme --ac 0 --mLow 105 --mHigh 140 --channel $channel
  else
    python makeWidthDatacards.py -b --writeoutput -i $inputcards -t $tpldir ${dataopt} -a $outdirname --coord $theCoords --GHmodel 0 --CatScheme $theCatScheme --ac 2 --mLow 105 --mHigh 140 --channel $channel
  fi

  RUN_STATUS=$?
  if [[ ${RUN_STATUS} -ne 0 ]]; then
    echo "Datacard compilation crashed  with error code ${RUN_STATUS}. Exiting..."
    exit ${RUN_STATUS}
  fi

  echo ""
done
