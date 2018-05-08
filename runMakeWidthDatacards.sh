#!/bin/bash


todaysdate=$1
curdir=$(pwd)

# Make off-shell datacards
for p in SM L1 a2 a3; do
  for data in 13TeV_2016 13TeV_2017; do
    outdirname=$todaysdate"_Offshell_"$p"_"$data
    outcardsname="cards_"$outdirname

    if [ ! -d $outcardsname ];then
      tpldir="templates2D/"$data"/"$p
      if [[ "$p" == "SM" ]];then
        python makeWidthDatacards.py -b -i Offshell_inputs_"$data" -t $tpldir -a $outdirname --coord "mass:KD1:KD2" --GHmodel 1 --CatScheme 2 --ac 0 --mLow 220 --mHigh 13000
      else
        python makeWidthDatacards.py -b -i Offshell_inputs_"$data" -t $tpldir -a $outdirname --coord "mass:KD1:KD2" --GHmodel 1 --CatScheme 2 --ac 2 --mLow 220 --mHigh 13000
      fi

      mkdir $outcardsname"/HCG/Scans/"

      cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
      cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
    fi
  done
done

# Make on-shell datacards
for p in SM L1 a2 a3; do
  for data in 13TeV_2016 13TeV_2017; do
    outdirname=$todaysdate"_Onshell_"$p"_"$data
    outcardsname="cards_"$outdirname

    if [ ! -d $outcardsname ];then
      tpldir="templates2D/"$data"/"$p
      if [[ "$p" == "SM" ]];then
        python makeWidthDatacards.py -b -i Onshell_inputs_AC_"$data" -t $tpldir -a $outdirname --coord "mass:KD1" --GHmodel 0 --CatScheme 2 --ac 0 --mLow 105 --mHigh 140
      else
        python makeWidthDatacards.py -b -i Onshell_inputs_AC_"$data" -t $tpldir -a $outdirname --coord "KD1:KD2:KD3" --GHmodel 0 --CatScheme 2 --ac 2 --mLow 105 --mHigh 140
      fi

      mkdir $outcardsname"/HCG/Scans/"

      cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
      cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
    fi
  done
done

# Make combination datacards
for p in SM L1 a2 a3; do
  outdirname=$todaysdate"_Combination_"$p
  outcardsname="cards_"$outdirname
  mkdir $outcardsname"/HCG/Scans/"

  for data in 13TeV_2016 13TeV_2017; do
    for region in Onshell Offshell; do
      incardsname="cards_"$todaysdate"_"$region"_"$p"_"$data"/HCG/"$data
      ln -sf $curdir"/"$incardsname $curdir"/"$outcardsname"/HCG/"$region"_"$data
      pushd $curdir"/"$outcardsname"/HCG"

      combineCards.py \
        ch1=Offshell_13TeV_2016/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2016/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2016/hzz4e_Untagged.txt \
        ch4=Offshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
        ch7=Offshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2016/hzz4e_HadVHTagged.txt \
        ch10=Offshell_13TeV_2016/hzz2e2mu_Untagged.txt ch11=Offshell_13TeV_2016/hzz4mu_Untagged.txt ch12=Offshell_13TeV_2016/hzz4e_Untagged.txt \
        ch13=Offshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch14=Offshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch15=Offshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
        ch16=Offshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch17=Offshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch18=Offshell_13TeV_2016/hzz4e_HadVHTagged.txt \
        > hzz4l_Prop_All_13TeV_2016.txt
      combineCards.py \
        ch1=Offshell_13TeV_2017/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2017/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2017/hzz4e_Untagged.txt \
        ch4=Offshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
        ch7=Offshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2017/hzz4e_HadVHTagged.txt \
        ch10=Offshell_13TeV_2017/hzz2e2mu_Untagged.txt ch11=Offshell_13TeV_2017/hzz4mu_Untagged.txt ch12=Offshell_13TeV_2017/hzz4e_Untagged.txt \
        ch13=Offshell_13TeV_2017/hzz2e2mu_JJVBFTagged.txt ch14=Offshell_13TeV_2017/hzz4mu_JJVBFTagged.txt ch15=Offshell_13TeV_2017/hzz4e_JJVBFTagged.txt \
        ch16=Offshell_13TeV_2017/hzz2e2mu_HadVHTagged.txt ch17=Offshell_13TeV_2017/hzz4mu_HadVHTagged.txt ch18=Offshell_13TeV_2017/hzz4e_HadVHTagged.txt \
        > hzz4l_Prop_All_13TeV_2017.txt
      combineCards.py hzz4l_Prop_All_13TeV_2016.txt hzz4l_Prop_All_13TeV_2017.txt > hzz4l_Prop_All_13TeV.txt
      if [[ "$p" != "SM" ]];then
        text2workspace.py \
          -m 125 hzz4l_Prop_All_13TeV.txt -o hzz4l_Prop_All_13TeV.root \
          -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
          --PO="offshell" --PO="allowPMF" --PO="sqrts=13"
      else
        text2workspace.py \
          -m 125 hzz4l_Prop_All_13TeV.txt -o hzz4l_Prop_All_13TeV.root \
          -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
          --PO="offshell" --PO="fai1fixed" --PO="sqrts=13"
      fi

      if [[ "$p" != "SM" ]];then
        ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/HIG-17-011/f"$p"/13_16TeV" Onshell_13TeV_2016_old
        ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/HIG-17-011/f"$p"/13_15TeV" Onshell_13TeV_2015
        ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/HIG-17-011/f"$p"/8TeV" Onshell_8TeV
        ln -sf "/work-zfs/lhc/usarica/hep/SpinWidthPaper_2015/HIG-17-011/f"$p"/7TeV" Onshell_7TeV

        combineCards.py \
          ch1=Offshell_13TeV_2016/hzz2e2mu_Untagged.txt ch2=Offshell_13TeV_2016/hzz4mu_Untagged.txt ch3=Offshell_13TeV_2016/hzz4e_Untagged.txt \
          ch4=Offshell_13TeV_2016/hzz2e2mu_JJVBFTagged.txt ch5=Offshell_13TeV_2016/hzz4mu_JJVBFTagged.txt ch6=Offshell_13TeV_2016/hzz4e_JJVBFTagged.txt \
          ch7=Offshell_13TeV_2016/hzz2e2mu_HadVHTagged.txt ch8=Offshell_13TeV_2016/hzz4mu_HadVHTagged.txt ch9=Offshell_13TeV_2016/hzz4e_HadVHTagged.txt \
          ch10=Onshell_13TeV_2016_old/hzz4l_4eS_Untagged_2016.lumi35.8671.txt ch11=Onshell_13TeV_2016_old/hzz4l_4eS_VBFtagged_2016.lumi35.8671.txt ch12=Onshell_13TeV_2016_old/hzz4l_4eS_VHHadrtagged_2016.lumi35.8671.txt \
          ch13=Onshell_13TeV_2016_old/hzz4l_2e2muS_Untagged_2016.lumi35.8671.txt ch14=Onshell_13TeV_2016_old/hzz4l_2e2muS_VBFtagged_2016.lumi35.8671.txt ch15=Onshell_13TeV_2016_old/hzz4l_2e2muS_VHHadrtagged_2016.lumi35.8671.txt \
          ch16=Onshell_13TeV_2016_old/hzz4l_4muS_Untagged_2016.lumi35.8671.txt ch17=Onshell_13TeV_2016_old/hzz4l_4muS_VBFtagged_2016.lumi35.8671.txt ch18=Onshell_13TeV_2016_old/hzz4l_4muS_VHHadrtagged_2016.lumi35.8671.txt \
          > hzz4l_Prop_All_13TeV_2016_old.txt
        head -n -2 hzz4l_Prop_All_13TeV_2016_old.txt > hzz4l_Prop_All_13TeV_2016_old.new.txt; mv hzz4l_Prop_All_13TeV_2016_old.new.txt hzz4l_Prop_All_13TeV_2016_old.txt

        combineCards.py hzz4l_Prop_All_13TeV_2016_old.txt hzz4l_Prop_All_13TeV_2017.txt > hzz4l_Prop_All_13TeV_old.txt
        text2workspace.py \
          -m 125 hzz4l_Prop_All_13TeV_old.txt -o hzz4l_Prop_All_13TeV_old.root \
          -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
          --PO="offshell" --PO="allowPMF" --PO="sqrts=13"
      fi

      popd

      cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
      cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"


      chmod -R 755 "$outcardsname"
    done
  done
done

