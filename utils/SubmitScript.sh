#!/bin/bash

for kind in {"Combined","4l","2l2n"}
do
  for itter in {0..39}
    do
    bsub -q 8nh -o lsflog_$kind_${itter}.txt -e lsferr_$kind_${itter}.err  Parallelize.lsf.sh $((itter*5)) $((itter*5+4)) $kind
  done
  for itter in {0..39}
    do
    bsub -q 8nh -o lsflogObs_$kind_${itter}.txt -e lsferrObs_$kind_${itter}.err  ParallelizeObs.lsf.sh $((itter*5)) $((itter*5+4)) $kind
  done
  for itter in {0..49}
    do
    bsub -q 8nh -o lsflogObs_lowR_${itter}.txt -e lsferrObs_lowR_${itter}.err  ParallelizeObs.lsf.sh $((itter*10)) $((itter*10+9)) $kind 2
  done
  for itter in {0..39}
    do
    bsub -q 8nh -o lsflog_lowR_${itter}.txt -e lsferr_lowR_${itter}.err  Parallelize.lsf.sh $((itter*5)) $((itter*5+4)) $kind 2
  done
done

