#!/bin/bash

todaysdate=$1
massregion=$2
hypo=$3
data=$4

theCatScheme=2
theCoords="KD1:KD2:KD3"
if [[ "$massregion" == "Offshell" ]]; then
  theCoords="mass:KD1:KD2"
  if [[ "$data" == *"2011"* ]] || [[ "$data" == *"2012"* ]] || [[ "$data" == *"2015"* ]];then 
    theCatScheme=1
  fi

  if [[ "$data" == *"2011"* ]] || [[ "$data" == *"2012"* ]];then 
    theCoords="mass:KD1"
  fi
else
  theCoords="KD1:KD2:KD3"
  if [[ "$hypo" == "SM" ]];then
    theCoords="mass:KD1"
  fi
  if [[ "$data" == *"2011"* ]] || [[ "$data" == *"2012"* ]] || [[ "$data" == *"2015"* ]];then 
    if [[ "$hypo" == "SM" ]];then
      theCatScheme=1
    else
      theCatScheme=0
    fi
  fi

  if [[ "$data" == *"2011"* ]] || [[ "$data" == *"2012"* ]];then 
    if [[ "$hypo" == "SM" ]];then
      theCoords="mass:KD2:KD1"
    else
      theCoords="KD1:KD2:KD3"
    fi
  elif [[ "$data" == *"2015"* ]];then 
    theCoords="KD1:KD2:KD3"
  fi
fi

echo "Cat. scheme: "$theCatScheme
echo "Coordinates: "$theCoords

if [[ "$massregion" == "Offshell" ]]; then
   outdirname=$todaysdate"_Offshell_"$hypo"_"$data
   outcardsname="cards_"$outdirname
   datadirname="CMSdata/"$data"/"$hypo"/Offshell/"
   if [ ! -d $outcardsname ];then
      tpldir="templates2D/"$data"/"$hypo
      if [[ "$hypo" == "SM" ]];then
         python makeWidthDatacards.py -b --writeoutput --dopdfproj -i Offshell_inputs_"$data" -t $tpldir -r $datadirname -a $outdirname --coord $theCoords --GHmodel 1 --CatScheme $theCatScheme --ac 0 --mLow 220 --mHigh 13000
      else
         python makeWidthDatacards.py -b --writeoutput --dopdfproj -i Offshell_inputs_"$data" -t $tpldir -r $datadirname -a $outdirname --coord $theCoords --GHmodel 1 --CatScheme $theCatScheme --ac 2 --mLow 220 --mHigh 13000
      fi

      mkdir $outcardsname"/HCG/Scans/"

      cp utils/submitScan1D.sh $outcardsname"/HCG/Scans/"
      cp utils/scan1D.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
   fi
else
   outdirname=$todaysdate"_Onshell_"$hypo"_"$data
   outcardsname="cards_"$outdirname
   datadirname="CMSdata/"$data"/"$hypo"/Onshell/"

   if [ ! -d $outcardsname ];then
      tpldir="templates2D/"$data"/"$hypo
      if [[ "$hypo" == "SM" ]];then
         extmsdir="externalShapes/"$data"/"
         python makeWidthDatacards.py -b --writeoutput -i Onshell_inputs_SM_"$data" --extMassShapes=$extmsdir -t $tpldir -r $datadirname -a $outdirname --coord $theCoords --GHmodel 0 --CatScheme $theCatScheme --ac 0 --mLow 105 --mHigh 140
      else
         python makeWidthDatacards.py -b --writeoutput --dopdfproj -i Onshell_inputs_AC_"$data" -t $tpldir -r $datadirname -a $outdirname --coord $theCoords --GHmodel 0 --CatScheme $theCatScheme --ac 2 --mLow 105 --mHigh 140
      fi

      mkdir $outcardsname"/HCG/Scans/"

      cp utils/submitScan1D.sh $outcardsname"/HCG/Scans/"
      cp utils/scan1D.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
   fi
fi
