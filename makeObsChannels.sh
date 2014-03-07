#!/bin/bash


python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d 0 -b 

#cd test1D/HCG/240 or test2D/HCG/240
cd cards_"$1"/HCG/220
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_allS_8TeV.txt ;
cd .. ;
#for i in 2e2mu 4e 4mu all ; do
for i in all ; do
cp -r 220 220_$i;
cd 220_$i;

text2workspace.py -m 220 hzz4l_"$i"S_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_allS_8TeV.root ;
combine -M GenerateOnly hzz4l_allS_8TeV.root -m 220  --saveToys -V -v 1 -n Obs;
root -b -l -q ../../../utils/addToyDataset.C\(\"hzz4l_allS_8TeV.root\",\"higgsCombineObs.GenerateOnly.mH220.123456.root\",\"toy_asimov\",\"workspaceWithAsimovObs.root\"\);
combine -M MultiDimFit workspaceWithAsimovObs.root --algo=grid --points 100 -m 220 -n 2D_exp -D toys/toy_asimov -v 3 #--setPhysicsModelParameterRanges CMS_zz4l_GGsm=0.0001,1 ;
root -l -q ../../../utils/plotScan1D.C\(220,30\) ; #-q;
cd ..;
done
