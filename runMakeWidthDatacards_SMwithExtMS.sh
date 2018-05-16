#!/bin/bash


todaysdate=$1
curdir=$(pwd)

# Make on-shell datacards
for data in 13TeV_2016 13TeV_2017; do
   outdirname=$todaysdate"_Onshell_SMwithExtMS_"$data
   outcardsname="cards_"$outdirname

   if [ ! -d $outcardsname ];then
      tpldir="templates2D/"$data"/"SM
      extmsdir="externalShapes/"$data
      python makeWidthDatacards.py -b -i Onshell_inputs_SM_"$data" --extMassShapes=$extmsdir -t $tpldir -a $outdirname --coord "mass:KD1" --GHmodel 0 --CatScheme 2 --ac 0 --mLow 105 --mHigh 140

      mkdir $outcardsname"/HCG/Scans/"

      cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
      cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
   fi
done

# Make combination datacards
outdirname=$todaysdate"_Combination_SMwithExtMS"
outcardsname="cards_"$outdirname

if [ ! -d $outcardsname ];then
   mkdir $outcardsname"/HCG/Scans/"

   for data in 13TeV_2016 13TeV_2017; do
      for region in Onshell; do
         incardsname="cards_"$todaysdate"_"$region"_SMwithExtMS_"$data"/HCG/"$data
         ln -sf $curdir"/"$incardsname $curdir"/"$outcardsname"/HCG/"$region"_"$data
      done
      for region in Offshell; do
         incardsname="cards_"$todaysdate"_"$region"_SM_"$data"/HCG/"$data
         ln -sf $curdir"/"$incardsname $curdir"/"$outcardsname"/HCG/"$region"_"$data
      done
   done

   pushd $curdir"/"$outcardsname"/HCG"

   combineCards.py \
      ch1=Offshell_13TeV_2016/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2016/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2016/hzz4e_Untagged.txt \
      ch4=Offshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
      ch7=Offshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2016/hzz4e_HadVHTagged.txt \
      ch10=Onshell_13TeV_2016/hzz2e2mu_Untagged.txt ch11=Onshell_13TeV_2016/hzz4mu_Untagged.txt ch12=Onshell_13TeV_2016/hzz4e_Untagged.txt \
      ch13=Onshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch14=Onshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch15=Onshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
      ch16=Onshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch17=Onshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch18=Onshell_13TeV_2016/hzz4e_HadVHTagged.txt \
      > hzz4l_Prop_All_13TeV_2016.txt
   combineCards.py \
      ch1=Offshell_13TeV_2017/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2017/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2017/hzz4e_Untagged.txt \
      ch4=Offshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
      ch7=Offshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2017/hzz4e_HadVHTagged.txt \
      ch10=Onshell_13TeV_2017/hzz2e2mu_Untagged.txt ch11=Onshell_13TeV_2017/hzz4mu_Untagged.txt ch12=Onshell_13TeV_2017/hzz4e_Untagged.txt \
      ch13=Onshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch14=Onshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch15=Onshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
      ch16=Onshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch17=Onshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch18=Onshell_13TeV_2017/hzz4e_HadVHTagged.txt \
      > hzz4l_Prop_All_13TeV_2017.txt
   combineCards.py hzz4l_Prop_All_13TeV_2016.txt hzz4l_Prop_All_13TeV_2017.txt > hzz4l_Prop_All_13TeV.txt
   rm -f hzz4l_Prop_All_13TeV.root
   text2workspace.py \
      -m 125 hzz4l_Prop_All_13TeV.txt -o hzz4l_Prop_All_13TeV.root \
      -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
      --PO="offshell" --PO="fai1fixed" --PO="sqrts=13"

   popd

   cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
   cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"


   chmod -R 755 "$outcardsname"
fi

