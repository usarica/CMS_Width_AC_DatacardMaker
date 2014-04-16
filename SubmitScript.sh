#!/bin/bash

#for r in 0.25 0.5 0.75 1; do 
#    for itter in {0..39}
#      do
#      bsub -q 8nh -o lsflog_obs_${itter}.txt -e lsferr_obs_${itter}.err  ParallelizeObs.lsf.sh $((itter*5)) $((itter*5+4)) $r 2
#    done
#done

#for r in 5 10 15 20; do
#    for itter in {0..39}
#      do
#      bsub -q 8nh -o lsflog_obs_${itter}.txt -e lsferr_obs_${itter}.err  ParallelizeObs.lsf.sh $((itter*5)) $((itter*5+4)) $r 30
#    done
#done


for itter in {0..39}
  do
  bsub -q 8nh -o lsflog_obs_${itter}.txt -e lsferr_obs_${itter}.err  ParallelizeObs.lsf.sh $((itter*5)) $((itter*5+4)) 30
done


for itter in {0..39}
  do
  bsub -q 8nh -o lsflog_obs_${itter}.txt -e lsferr_obs_${itter}.err  Parallelize.lsf.sh $((itter*5)) $((itter*5+4)) 1 30
done