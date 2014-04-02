#!/bin/bash

echo "python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d $2 -b "

python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d $2 -b 
python make_width2D_DCsandWSs.py -i SM_inputs_7TeV -a $1_7 -d $2 -b 

mv cards_"$1"_7/HCG/220/*.* cards_"$1"/HCG/220/.
rm -rf cards_"$1"_7

cd cards_"$1"/HCG/
cp -r 220 220_mu
cd 220_mu
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_allS_8TeV.txt

text2workspace.py -m 220 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth --PO=RVRFfixed -o hzz4l_allS_8TeV.root

combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 -V --algo=grid --points 50 -n 2D_exp -v 3
combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -V --algo=grid --points 50 -n 2D_obs -v 3

root -l -q ../../../utils/plotScan1D.C\(220,30\)


cd ..
cp -r 220 220_mu_noSyst
cd 220_mu_noSyst
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_allS_8TeV.txt

text2workspace.py -m 220 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth --PO=RVRFfixed -o hzz4l_allS_8TeV.root

combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 -V --algo=grid --points 200 -n 2D_exp --setPhysicsModelParameters CMS_zz4l_GGsm=1,R=1 --freezeNuis CMS_eff_e,CMS_eff_m,lumi_8TeV,pdf_qqbar,QCDscale_VV,CMS_widthH_kbkg,CMS_zz4l_VBFscale_syst,QCDscale_ggH,pdf_gg,CMS_zz4l_ZXshape_syst,CMS_hzz2e2mu_Zjets,CMS_hzz4e_Zjets,CMS_hzz4mu_Zjets,R,RV,RF -v 3

root -l -q ../../../utils/plotScan1D.C\(220,30\)


cp -r 220 220_muVmuF
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_allS_8TeV.txt

text2workspace.py -m 220 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_allS_8TeV.root

combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 -V --algo=grid --points 200 -n 2D_exp -v 3

root -l -q ../../../utils/plotScan1D.C\(220,30\)


cd ..
cp -r 220 220_muVmuF_noSyst
cd 220_mu_noSyst
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_allS_8TeV.txt

text2workspace.py -m 220 hzz4l_allS_8TeV.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_allS_8TeV.root

combine -M MultiDimFit hzz4l_allS_8TeV.root -m 220 -t -1 --expectSignal=1 -V --algo=grid --points 200 -n 2D_exp --setPhysicsModelParameters CMS_zz4l_GGsm=1,RV=1,RF=1 --freezeNuis CMS_eff_e,CMS_eff_m,lumi_8TeV,pdf_qqbar,QCDscale_VV,CMS_widthH_kbkg,CMS_zz4l_VBFscale_syst,QCDscale_ggH,pdf_gg,CMS_zz4l_ZXshape_syst,CMS_hzz2e2mu_Zjets,CMS_hzz4e_Zjets,CMS_hzz4mu_Zjets,R,RV,RF -v 3

root -l -q ../../../utils/plotScan1D.C\(220,30\)
