#!/bin/bash


todaysdate=$1

# Make off-shell datacards
for p in SM L1 a2 a3; do
  for data in 13TeV_2016 13TeV_2017; do
    outdirname=$todaysdate"_Offshell_"$p"_"$data
    outcardsname="cards_"$outdirname
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
  done
done

# Make on-shell datacards
for p in SM L1 a2 a3; do
  for data in 13TeV_2016 13TeV_2017; do
    outdirname=$todaysdate"_Onshell_"$p"_"$data
    outcardsname="cards_"$outdirname
    tpldir="templates2D/"$data"/"$p

    if [[ "$p" == "SM" ]];then
      python makeWidthDatacards.py -b -i Onshell_inputs_AC_"$data" -t $tpldir -a $outdirname --coord "mass:KD1:KD2" --GHmodel 0 --CatScheme 2 --ac 2 --mLow 105 --mHigh 140
    else
      python makeWidthDatacards.py -b -i Onshell_inputs_AC_"$data" -t $tpldir -a $outdirname --coord "KD1:KD2:KD3" --GHmodel 0 --CatScheme 2 --ac 2 --mLow 105 --mHigh 140
    fi

    mkdir $outcardsname"/HCG/Scans/"

    cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
    cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"

    chmod -R 755 "$outcardsname"
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
    ln -sf $incardsname $outcardsname"/HCG/"$data"_"$region

    cp utils/submitScanExp.sh $outcardsname"/HCG/Scans/"
    cp utils/scanExp.slurm.sh $outcardsname"/HCG/Scans/"

    chmod -R 755 "$outcardsname"
  done
  done
done

