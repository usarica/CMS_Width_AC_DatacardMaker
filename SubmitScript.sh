#!/bin/bash

mkdir Parallel
cd Parallel
for itter in {0..99}
do
  mkdir Par${itter}
  cd Par${itter}
  cp ../../Parallelize.lsf.sh .
  cp ../../cards_03_05_Unblind_093_1DDgg/HCG/220/hzz4l_*S_8TeV.* .
  bsub -q 1nh -o lsflog_${itter}.txt -e lsferr_${itter}.err -R "type=SLC5_64"  Parallelize.lsf.sh $itter
  cd ..
done
