#!/bin/bash

date=$1
tpltag=$2
dcinputstag=$3
datainputtag=$4

npara=1
if [[ "$5" != "" ]]; then
  npara=$5
fi

run(){
  local hasBoostedVH=$1
  local hypo=$2
  local year=$3

  cmd="./makeDatacards_ZWTo3L1Nu.run.sh $date $hypo $year $tpltag $dcinputstag $datainputtag $hasBoostedVH"
  $cmd
  RUN_STATUS=$?
  if [[ ${RUN_STATUS} -ne 0 ]]; then
    echo "$cmd crashed with status ${RUN_STATUS}."  
    exit ${RUN_STATUS}
  fi
}
job_limit(){
  local joblist=( )
  # Test for single positive integer input
  if [[ $# -eq 1 ]] && [[ $1 =~ ^[1-9][0-9]*$ ]]
  then
    # Check number of running jobs
    joblist=( $(jobs -rp) )
    while [[ ${#joblist[*]} -ge $1 ]]; do
      # Wait for any job to finish
      local command='wait '${joblist[0]}
      for job in ${joblist[@]:1}; do
        command+=' || wait '$job
      done
      eval ${command}
      joblist=( $(jobs -rp) )
    done
  fi
}

for hasBoostedVH in 0; do
  for hypo in SM a2 a3 L1; do

    for year in 2016 2017 2018; do
      if [[ $npara -eq 1 ]]; then
        run $hasBoostedVH $hypo $year
      else
        run $hasBoostedVH $hypo $year &
      fi
      job_limit $npara
    done
  done
done

wait

for hasBoostedVH in 0; do
  for hypo in SM a2 a3 L1; do
    theCatScheme="nj012_3l1nu"
    dirname=cards_${date}_ZW3L1Nu_InputTag_${dcinputstag}_TplTag_${tpltag}_${hypo}_${theCatScheme}/HCG

    if [[ ! -d $dirname ]]; then
      continue
    fi

    mkdir -p ${dirname}/13TeV
    for f in loadLib.C buildCards.run.sh; do
      cp utils/$f ${dirname}/13TeV/
    done

    cd ${dirname}/13TeV
    combineCards.py ../13TeV_201?/*.txt > hto3l1nu_allcats_13TeV.txt
    ./buildCards.run.sh hto3l1nu_allcats_13TeV.txt hto3l1nu_allcats_13TeV.root 13 $hypo
    if [[ "$hypo" == "SM" ]]; then
      ./buildCards.run.sh hto3l1nu_allcats_13TeV.txt hto3l1nu_allcats_13TeV_RVoverRF.root 13 "${hypo}.UseRVoverRF"
    fi
    for year in 2016 2017 2018; do
      combineCards.py ../13TeV_${year}/*.txt > hto3l1nu_allcats_13TeV_${year}.txt
      ./buildCards.run.sh hto3l1nu_allcats_13TeV_${year}.txt hto3l1nu_allcats_13TeV_${year}.root 13 $hypo
    done
    cd -

  done
done
