#!/bin/bash

echo "python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d $2 -b "

python make_width2D_DCsandWSs.py -i SM_inputs_8TeV -a $1 -d $2 -b 
python make_width2D_DCsandWSs.py -i SM_inputs_7TeV -a $1_7 -d $2 -b 

mv cards_"$1"_7/HCG/220/*.* cards_"$1"/HCG/220/.
rm -rf cards_"$1"_7

cd cards_"$1"/HCG/220

if  [ "$2" == "2" ]
then 
ln -ns  ~/work/svnDatacards/HIG-14-002/hzz2l2nu/*.dat .
ln -ns  ~/work/svnDatacards/HIG-14-002/hzz2l2nu/hzz2l2v*.root .
fi
ln -ns ~/work/svnDatacards/HIG-14-002/hzz4l/*_0.* .
ln -ns ~/work/svnDatacards/HIG-14-002/hzz4l/*_1.* .

#combine cards ->hzz_all.txt, hzz4l_all.txt, hzz2l2n_all.txt

#combineCards.py hzz4l_4muS_7TeV_0.txt hzz4l_4muS_8TeV_0.txt hzz4l_4eS_7TeV_0.txt hzz4l_4eS_8TeV_0.txt hzz4l_2e2muS_7TeV_0.txt hzz4l_2e2muS_8TeV_0.txt hzz4l_4muS_7TeV_1.txt hzz4l_4muS_8TeV_1.txt hzz4l_4eS_7TeV_1.txt hzz4l_4eS_8TeV_1.txt hzz4l_2e2muS_7TeV_1.txt hzz4l_2e2muS_8TeV_1.txt > hzz4l_low.txt
combineCards.py hzz4l_2e2muS_7TeV_0.txt  hzz4l_2e2muS_8TeV_0.txt  hzz4l_4muS_7TeV_0.txt  hzz4l_4muS_8TeV_0.txt  hzz4l_4eS_7TeV_0.txt  hzz4l_4eS_8TeV_0.txt  hzz4l_2e2muS_7TeV_1.txt  hzz4l_2e2muS_8TeV_1.txt  hzz4l_4muS_7TeV_1.txt  hzz4l_4muS_8TeV_1.txt  hzz4l_4eS_7TeV_1.txt  hzz4l_4eS_8TeV_1.txt > hzz4l_low.txt

combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_2e2muS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_4eS_7TeV.txt hzz4l_2e2muS_7TeV_0.txt hzz4l_2e2muS_8TeV_0.txt hzz4l_4muS_7TeV_0.txt hzz4l_4muS_8TeV_0.txt hzz4l_4eS_7TeV_0.txt hzz4l_4eS_8TeV_0.txt hzz4l_2e2muS_7TeV_1.txt hzz4l_2e2muS_8TeV_1.txt hzz4l_4muS_7TeV_1.txt hzz4l_4muS_8TeV_1.txt hzz4l_4eS_7TeV_1.txt hzz4l_4eS_8TeV_1.txt > hzz4l_all.txt

#combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_2e2muS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_4eS_7TeV.txt hzz4l_2e2muS_7TeV_0.txt hzz4l_2e2muS_8TeV_0.txt hzz4l_4muS_7TeV_0.txt hzz4l_4muS_8TeV_0.txt hzz4l_4eS_7TeV_0.txt hzz4l_4eS_8TeV_0.txt hzz4l_2e2muS_7TeV_1.txt hzz4l_2e2muS_8TeV_1.txt hzz4l_4muS_7TeV_1.txt hzz4l_4muS_8TeV_1.txt hzz4l_4eS_7TeV_1.txt hzz4l_4eS_8TeV_1.txt card_combined.dat > hzz_all.txt


#workspaces
text2workspace.py -m 125.6 hzz4l_all.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz4l_all.root
if  [ "$2" == "2" ]
then 
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_2e2muS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_4eS_7TeV.txt hzz4l_2e2muS_7TeV_0.txt hzz4l_2e2muS_8TeV_0.txt hzz4l_4muS_7TeV_0.txt hzz4l_4muS_8TeV_0.txt hzz4l_4eS_7TeV_0.txt hzz4l_4eS_8TeV_0.txt hzz4l_2e2muS_7TeV_1.txt hzz4l_2e2muS_8TeV_1.txt hzz4l_4muS_7TeV_1.txt hzz4l_4muS_8TeV_1.txt hzz4l_4eS_7TeV_1.txt hzz4l_4eS_8TeV_1.txt card_combined.dat > hzz_all.txt
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_2e2muS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_4eS_7TeV.txt hzz4l_2e2muS_7TeV_0.txt hzz4l_2e2muS_8TeV_0.txt hzz4l_4muS_7TeV_0.txt hzz4l_4muS_8TeV_0.txt hzz4l_4eS_7TeV_0.txt hzz4l_4eS_8TeV_0.txt hzz4l_2e2muS_7TeV_1.txt hzz4l_2e2muS_8TeV_1.txt hzz4l_4muS_7TeV_1.txt hzz4l_4muS_8TeV_1.txt hzz4l_4eS_7TeV_1.txt hzz4l_4eS_8TeV_1.txt card_incl.dat > hzzincl_all.txt
combineCards.py hzz4l_low.txt card_combined.dat > hzz2l2n_all.txt
combineCards.py hzz4l_low.txt card_incl.dat > hzz2l2n_incl.txt
text2workspace.py -m 125.6 hzz2l2n_incl.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth --PO=is2l2nu -o hzz2l2n_all.root
text2workspace.py -m 125.6 hzz_all.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzz_all.root
text2workspace.py -m 125.6 hzzincl_all.txt -P HiggsAnalysis.CombinedLimit.HiggsWidth:higgswidth -o hzzincl_all.root
fi
