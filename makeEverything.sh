#!/bin/bash

echo "python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d $2 -b "

python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d $2 -b 
python make_width2D_DCsandWSs.py -i SM_inputs_7TeV -a $1_7 -d $2 -b 

mv cards_"$1"_7/HCG/220/*.* cards_"$1"/HCG/220/.
rm -rf cards_"$1"_7

#cd test1D/HCG/240 or test2D/HCG/240
cd cards_"$1"/HCG/220
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_2e2muS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_4eS_7TeV.txt > hzz4l_allS_8TeV.txt

#text2workspace.py -m 240 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_allS_8TeV.root
text2workspace.py -m 220 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_allS_8TeV.root
#(add --stat if no systematics)
#(add -PO=GGsmVal=25 if you want to run with e.g. G/G_SM = 25 and not 1)

#combine -M GenerateOnly hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 --saveToys -V -v 7
#combine -M GenerateOnly hzz4l_allS_8TeV.root -m 220  --saveToys -V -v 1 -n Obs
#root -b -l -q ../../../utils/addToyDataset.C\(\"hzz4l_allS_8TeV.root\",\"higgsCombineTest.GenerateOnly.mH220.123456.root\",\"toy_asimov\",\"workspaceWithAsimov.root\"\)
#root -b -l -q ../../../utils/addToyDataset.C\(\"hzz4l_allS_8TeV.root\",\"higgsCombineObs.GenerateOnly.mH220.123456.root\",\"toy_asimov\",\"workspaceWithAsimovObs.root\"\)
#combine -M MultiDimFit workspaceWithAsimov.root --algo=grid --points 200 -m 220 -n 2D_exp -D toys/toy_asimov -v 3
#combine -M MultiDimFit workspaceWithAsimovObs.root --algo=grid --points 200 -m 220 -n 2D_obs -D toys/toy_asimov -v 3
#(add "-S 0 --fastScan" if no systematics)

combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 -V --algo=grid --points 50 -n 2D_exp -v 3
combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -V --algo=grid --points 50 -n 2D_obs -v 3

root -l -q ../../../utils/plotScan1D.C\(220,30\)  #-q


cd ..
cp -r 220 220_noSyst
cd 220_noSyst
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_allS_8TeV.txt

text2workspace.py -m 220 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_allS_8TeV.root
#combine -M GenerateOnly hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 --saveToys -V -v 7
#root -b -l -q ../../../utils/addToyDataset.C\(\"hzz4l_allS_8TeV.root\",\"higgsCombineTest.GenerateOnly.mH220.123456.root\",\"toy_asimov\",\"workspaceWithAsimov.root\"\)
#combine -M MultiDimFit workspaceWithAsimov.root --algo=grid --points 200 -m 220 -n 2D_exp -D toys/toy_asimov -v 3 -S 0 --fastScan

combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 -V --algo=grid --points 200 -n 2D_exp --setPhysicsModelParameters CMS_zz4l_GGsm=1,CMS_zz4l_mu=1 --freezeNuis CMS_eff_e,CMS_eff_m,lumi_8TeV,pdf_qqbar,CMS_QCDscale_VV,CMS_widthH_kbkg,CMS_zz4l_VBFscale_syst,CMS_zz4l_pdf_QCDscale_gg_syst,CMS_zz4l_mu,CMS_zz4l_ZXshape_syst,CMS_hzz2e2mu_Zjets,CMS_hzz4e_Zjets,CMS_hzz4mu_Zjets -v 3

root -l -q ../../../utils/plotScan1D.C\(220,30\)
