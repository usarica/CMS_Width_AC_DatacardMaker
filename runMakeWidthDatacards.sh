#!/bin/bash


todaysdate=$1
option=$2
curdir=$(pwd)

if [[ "$option" == *"offshell"* ]];then
echo "Doing off-shell cards"
# Make off-shell datacards
for p in SM L1 a2 a3; do
  for data in 13TeV_2016 13TeV_2017; do
    . makeDatacards.sh $todaysdate "Offshell" $p $data
  done
done
fi

if [[ "$option" == *"onshell"* ]];then
echo "Doing on-shell cards"
# Make on-shell datacards
for p in SM L1 L1ZGs a2 a3; do
  for data in 13TeV_2016 13TeV_2017; do
    . makeDatacards.sh $todaysdate "Onshell" $p $data
  done
done
fi

if [[ "$option" == *"combine"* ]];then
echo "Combining cards"
# Make combination datacards
for p in L1 L1ZGs a2 a3; do
   outdirname=$todaysdate"_Combination_"$p
   outcardsname="cards_"$outdirname

   if [ ! -d $outcardsname ];then
      mkdir $outcardsname"/HCG/Scans/"

      for data in 13TeV_2016 13TeV_2017; do
         for region in Onshell Offshell; do
            incardsname="cards_"$todaysdate"_"$region"_"$p"_"$data"/HCG/"$data
            if [ -d $curdir"/"$incardsname ];then
               ln -sf $curdir"/"$incardsname $curdir"/"$outcardsname"/HCG/"$region"_"$data
            fi
         done
      done

      cp utils/buildCards.sh $outcardsname"/HCG/"
      cp utils/buildCards.slurm.sh $outcardsname"/HCG/"

      pushd $curdir"/"$outcardsname"/HCG"

      if [[ "$p" == "L1ZGs" ]]; then
         combineCards.py \
            ch10=Onshell_13TeV_2016/hzz2e2mu_Untagged.txt ch11=Onshell_13TeV_2016/hzz4mu_Untagged.txt ch12=Onshell_13TeV_2016/hzz4e_Untagged.txt \
            ch13=Onshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch14=Onshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch15=Onshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
            ch16=Onshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch17=Onshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch18=Onshell_13TeV_2016/hzz4e_HadVHTagged.txt \
            > hzz4l_Prop_All_13TeV_2016.txt
         combineCards.py \
            ch10=Onshell_13TeV_2017/hzz2e2mu_Untagged.txt ch11=Onshell_13TeV_2017/hzz4mu_Untagged.txt ch12=Onshell_13TeV_2017/hzz4e_Untagged.txt \
            ch13=Onshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch14=Onshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch15=Onshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
            ch16=Onshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch17=Onshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch18=Onshell_13TeV_2017/hzz4e_HadVHTagged.txt \
            > hzz4l_Prop_All_13TeV_2017.txt

         combineCards.py hzz4l_Prop_All_13TeV_2016.txt hzz4l_Prop_All_13TeV_2017.txt > hzz4l_Prop_All_13TeV.txt
      else
         combineCards.py \
            ch10=Onshell_13TeV_2016/hzz2e2mu_Untagged.txt ch11=Onshell_13TeV_2016/hzz4mu_Untagged.txt ch12=Onshell_13TeV_2016/hzz4e_Untagged.txt \
            ch13=Onshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch14=Onshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch15=Onshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
            ch16=Onshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch17=Onshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch18=Onshell_13TeV_2016/hzz4e_HadVHTagged.txt \
            > hzz4l_Prop_All_Onshell_13TeV_2016.txt
         combineCards.py \
            ch10=Onshell_13TeV_2017/hzz2e2mu_Untagged.txt ch11=Onshell_13TeV_2017/hzz4mu_Untagged.txt ch12=Onshell_13TeV_2017/hzz4e_Untagged.txt \
            ch13=Onshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch14=Onshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch15=Onshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
            ch16=Onshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch17=Onshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch18=Onshell_13TeV_2017/hzz4e_HadVHTagged.txt \
            > hzz4l_Prop_All_Onshell_13TeV_2017.txt
         combineCards.py \
            ch1=Offshell_13TeV_2016/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2016/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2016/hzz4e_Untagged.txt \
            ch4=Offshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
            ch7=Offshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2016/hzz4e_HadVHTagged.txt \
            > hzz4l_Prop_All_Offshell_13TeV_2016.txt
         combineCards.py \
            ch1=Offshell_13TeV_2017/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2017/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2017/hzz4e_Untagged.txt \
            ch4=Offshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
            ch7=Offshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2017/hzz4e_HadVHTagged.txt \
            > hzz4l_Prop_All_Offshell_13TeV_2017.txt

         combineCards.py hzz4l_Prop_All_Onshell_13TeV_2016.txt hzz4l_Prop_All_Onshell_13TeV_2017.txt > hzz4l_Prop_All_Onshell_13TeV.txt
         combineCards.py hzz4l_Prop_All_Offshell_13TeV_2016.txt hzz4l_Prop_All_Offshell_13TeV_2017.txt > hzz4l_Prop_All_Offshell_13TeV.txt
         combineCards.py hzz4l_Prop_All_Onshell_13TeV.txt hzz4l_Prop_All_Offshell_13TeV.txt > hzz4l_Prop_All_13TeV.txt

         rm -f hzz4l_Prop_All_Onshell_13TeV.root
         . buildCards.sh hzz4l_Prop_All_Onshell_13TeV.txt hzz4l_Prop_All_Onshell_13TeV.root "13" $p
      fi

      rm -f hzz4l_Prop_All_13TeV.root
      . buildCards.sh hzz4l_Prop_All_13TeV.txt hzz4l_Prop_All_13TeV.root "13" $p

      if [[ "$p" != "SM" ]];then
         ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/HIG-18-002/HeshyDCs/f"$p Onshell_13TeV_2017_hronshell
         ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/HIG-18-002/HeshyDCs/f"$p Onshell_13TeV_2016_hronshell
         ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/RunII_Combination/CMSSW_9_4_3/src/HZZ4l/CreateWidthDatacards/cards_180817_Onshell_"$p"_7TeV_2011/HCG/7TeV_2011" Onshell_7TeV_2011
         ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/RunII_Combination/CMSSW_9_4_3/src/HZZ4l/CreateWidthDatacards/cards_180817_Onshell_"$p"_8TeV_2012/HCG/8TeV_2012" Onshell_8TeV_2012
         ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/RunII_Combination/CMSSW_9_4_3/src/HZZ4l/CreateWidthDatacards/cards_180817_Onshell_"$p"_13TeV_2015/HCG/13TeV_2015" Onshell_13TeV_2015

         combineCards.py \
            ch1=Onshell_7TeV_2011/hzz2e2mu_Inclusive.txt ch2=Onshell_7TeV_2011/hzz4mu_Inclusive.txt ch3=Onshell_7TeV_2011/hzz4e_Inclusive.txt \
            > hzz4l_Prop_All_7TeV_2011.txt
         combineCards.py \
            ch1=Onshell_8TeV_2012/hzz2e2mu_Inclusive.txt ch2=Onshell_8TeV_2012/hzz4mu_Inclusive.txt ch3=Onshell_8TeV_2012/hzz4e_Inclusive.txt \
            > hzz4l_Prop_All_8TeV_2012.txt
         combineCards.py \
            ch1=Onshell_13TeV_2015/hzz2e2mu_Inclusive.txt ch2=Onshell_13TeV_2015/hzz4mu_Inclusive.txt ch3=Onshell_13TeV_2015/hzz4e_Inclusive.txt \
            > hzz4l_Prop_All_13TeV_2015.txt
         combineCards.py hzz4l_Prop_All_7TeV_2011.txt hzz4l_Prop_All_8TeV_2012.txt hzz4l_Prop_All_13TeV_2015.txt > hzz4l_Prop_7813_15.txt

         if [[ "$p" != "L1ZGs" ]]; then
            combineCards.py \
               ch1=Offshell_13TeV_2016/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2016/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2016/hzz4e_Untagged.txt \
               ch4=Offshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
               ch7=Offshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2016/hzz4e_HadVHTagged.txt \
               ch10=Onshell_13TeV_2016_hronshell/hzz4l_4eS_Untagged_2016.lumi35.92.txt ch11=Onshell_13TeV_2016_hronshell/hzz4l_4eS_VBFtagged_2016.lumi35.92.txt ch12=Onshell_13TeV_2016_hronshell/hzz4l_4eS_VHHadrtagged_2016.lumi35.92.txt \
               ch13=Onshell_13TeV_2016_hronshell/hzz4l_2e2muS_Untagged_2016.lumi35.92.txt ch14=Onshell_13TeV_2016_hronshell/hzz4l_2e2muS_VBFtagged_2016.lumi35.92.txt ch15=Onshell_13TeV_2016_hronshell/hzz4l_2e2muS_VHHadrtagged_2016.lumi35.92.txt \
               ch16=Onshell_13TeV_2016_hronshell/hzz4l_4muS_Untagged_2016.lumi35.92.txt ch17=Onshell_13TeV_2016_hronshell/hzz4l_4muS_VBFtagged_2016.lumi35.92.txt ch18=Onshell_13TeV_2016_hronshell/hzz4l_4muS_VHHadrtagged_2016.lumi35.92.txt \
               > hzz4l_Prop_All_13TeV_2016_hronshell.txt
            combineCards.py \
               ch1=Offshell_13TeV_2017/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2017/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2017/hzz4e_Untagged.txt \
               ch4=Offshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
               ch7=Offshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2017/hzz4e_HadVHTagged.txt \
               ch10=Onshell_13TeV_2017_hronshell/hzz4l_4eS_Untagged_2017.lumi41.53.txt ch11=Onshell_13TeV_2017_hronshell/hzz4l_4eS_VBFtagged_2017.lumi41.53.txt ch12=Onshell_13TeV_2017_hronshell/hzz4l_4eS_VHHadrtagged_2017.lumi41.53.txt \
               ch13=Onshell_13TeV_2017_hronshell/hzz4l_2e2muS_Untagged_2017.lumi41.53.txt ch14=Onshell_13TeV_2017_hronshell/hzz4l_2e2muS_VBFtagged_2017.lumi41.53.txt ch15=Onshell_13TeV_2017_hronshell/hzz4l_2e2muS_VHHadrtagged_2017.lumi41.53.txt \
               ch16=Onshell_13TeV_2017_hronshell/hzz4l_4muS_Untagged_2017.lumi41.53.txt ch17=Onshell_13TeV_2017_hronshell/hzz4l_4muS_VBFtagged_2017.lumi41.53.txt ch18=Onshell_13TeV_2017_hronshell/hzz4l_4muS_VHHadrtagged_2017.lumi41.53.txt \
               > hzz4l_Prop_All_13TeV_2017_hronshell.txt
         else
            combineCards.py \
               ch10=Onshell_13TeV_2016_hronshell/hzz4l_4eS_Untagged_2016.lumi35.92.txt ch11=Onshell_13TeV_2016_hronshell/hzz4l_4eS_VBFtagged_2016.lumi35.92.txt ch12=Onshell_13TeV_2016_hronshell/hzz4l_4eS_VHHadrtagged_2016.lumi35.92.txt \
               ch13=Onshell_13TeV_2016_hronshell/hzz4l_2e2muS_Untagged_2016.lumi35.92.txt ch14=Onshell_13TeV_2016_hronshell/hzz4l_2e2muS_VBFtagged_2016.lumi35.92.txt ch15=Onshell_13TeV_2016_hronshell/hzz4l_2e2muS_VHHadrtagged_2016.lumi35.92.txt \
               ch16=Onshell_13TeV_2016_hronshell/hzz4l_4muS_Untagged_2016.lumi35.92.txt ch17=Onshell_13TeV_2016_hronshell/hzz4l_4muS_VBFtagged_2016.lumi35.92.txt ch18=Onshell_13TeV_2016_hronshell/hzz4l_4muS_VHHadrtagged_2016.lumi35.92.txt \
               > hzz4l_Prop_All_13TeV_2016_hronshell.txt
            combineCards.py \
               ch10=Onshell_13TeV_2017_hronshell/hzz4l_4eS_Untagged_2017.lumi41.53.txt ch11=Onshell_13TeV_2017_hronshell/hzz4l_4eS_VBFtagged_2017.lumi41.53.txt ch12=Onshell_13TeV_2017_hronshell/hzz4l_4eS_VHHadrtagged_2017.lumi41.53.txt \
               ch13=Onshell_13TeV_2017_hronshell/hzz4l_2e2muS_Untagged_2017.lumi41.53.txt ch14=Onshell_13TeV_2017_hronshell/hzz4l_2e2muS_VBFtagged_2017.lumi41.53.txt ch15=Onshell_13TeV_2017_hronshell/hzz4l_2e2muS_VHHadrtagged_2017.lumi41.53.txt \
               ch16=Onshell_13TeV_2017_hronshell/hzz4l_4muS_Untagged_2017.lumi41.53.txt ch17=Onshell_13TeV_2017_hronshell/hzz4l_4muS_VBFtagged_2017.lumi41.53.txt ch18=Onshell_13TeV_2017_hronshell/hzz4l_4muS_VHHadrtagged_2017.lumi41.53.txt \
               > hzz4l_Prop_All_13TeV_2017_hronshell.txt
         fi

         combineCards.py hzz4l_Prop_All_13TeV_2016_hronshell.txt hzz4l_Prop_All_13TeV_2017_hronshell.txt > hzz4l_Prop_All_13TeV_hronshell.txt
         . buildCards.sh hzz4l_Prop_All_13TeV_hronshell.txt hzz4l_Prop_All_13TeV_hronshell.root "13" $p

         combineCards.py hzz4l_Prop_7813_15.txt hzz4l_Prop_All_13TeV_hronshell.txt > hzz4l_Prop_All_7813TeV_hronshell.txt
         . buildCards.sh hzz4l_Prop_All_7813TeV_hronshell.txt hzz4l_Prop_All_7813TeV_hronshell.root "7,8,13" $p

         combineCards.py hzz4l_Prop_7813_15.txt hzz4l_Prop_All_13TeV.txt > hzz4l_Prop_All_7813TeV.txt
         . buildCards.sh hzz4l_Prop_All_7813TeV.txt hzz4l_Prop_All_7813TeV.root "7,8,13" $p
      fi

      popd

      cp utils/submitScan1D.sh $outcardsname"/HCG/Scans/"
      cp utils/scan1D.slurm.sh $outcardsname"/HCG/Scans/"
      cp utils/submitImpacts.sh $outcardsname"/HCG/Scans/"
      cp utils/doImpacts.slurm.sh $outcardsname"/HCG/Scans/"
      cp utils/submitPreApprovalChecks.sh $outcardsname"/HCG/Scans/"
      cp utils/doPreApprovalChecks.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
   fi
done
fi
