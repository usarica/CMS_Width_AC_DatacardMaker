#!/bin/bash

ddate=$1
npara=1
if [[ "$2" != "" ]]; then
  npara=$2
fi

run(){
  local date=$1
  local hypo=$2
  local year=$3

  ./makeDatacards_ZZTo4L_19009.run.sh $date $hypo $year
  RUN_STATUS=$?
  if [[ ${RUN_STATUS} -ne 0 ]]; then
    echo "Compilation of $hypo $year failed."
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

for hhypo in SM a3 a2 L1; do
  for yyear in 2016 2017 2018; do
    run $ddate $hhypo $yyear &
    job_limit $npara
  done
done

wait

for hhypo in SM a3 a2 L1; do
  for yyear in 2016 2017 2018; do
    dirname=cards_${ddate}_ZZ4L_${hhypo}_Onshell_19009
    dirname=${dirname}/HCG

    mkdir -p ${dirname}/13TeV
    for f in loadLib.C buildCards.run.sh; do
      cp utils/$f ${dirname}/13TeV/
    done

    cd ${dirname}/13TeV
    combineCards.py ../13TeV_201?/*.txt > hto4l_allcats_13TeV.txt
    cd -
  done
done
