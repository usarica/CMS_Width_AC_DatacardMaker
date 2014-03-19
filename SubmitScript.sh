#!/bin/bash

mkdir toyParallel2D_093_95_postfit
cd toyParallel2D_093_95_postfit
for itter in {0..99}
do
  mkdir Par${itter}
  cd Par${itter}
  cp ../../Parallelize.lsf.sh .
  cp ../../cards_03_17_Moriond_093_2D/HCG/220/hzz4l_*S_8TeV.* .
  #bsub -q 8nh -o lsflog_${itter}.txt -e lsferr_${itter}.err -R "type=SLC5_64"  Parallelize.lsf.sh $itter
  bsub -q 1nh -o lsflog_${itter}.txt -e lsferr_${itter}.err Parallelize.lsf.sh $itter
  cd ..
done
