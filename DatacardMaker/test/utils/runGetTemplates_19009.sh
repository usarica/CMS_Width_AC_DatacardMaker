#!/bin/bash

npara=$1

run(){
  local run_args="$1"
  root -b -l -q loadLib.C "${run_args}"
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


for f in $(ls HVV | grep .txt); do
  year=""
  if [[ "$f" == *"2016"* ]]; then
    year="2016"
  elif [[ "$f" == *"2017"* ]]; then
    year="2017"
  elif [[ "$f" == *"2018"* ]]; then
    year="2018"
  fi
  strout="Onshell_13TeV_${year}"
  for hypo in SM a2 a3 L1; do
    coutput_main="${hypo}/${strout}"

    run "getTemplates_19009.cc+(\"HVV/${f}\",\"${coutput_main}\",\"${hypo}\",-1,false,true)" &
    job_limit $npara
  done
done

wait
