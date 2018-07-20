#!/bin/bash

#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=lrgmem
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
#SBATCH --mail-user=usarica1@jhu.edu

cd ${SLURM_SUBMIT_DIR}
echo "SLURM job running in: "$(pwd)

source /work-zfs/lhc/cms/cmsset_default.sh
module load boost/1.60.0
export LIBRARY_PATH=$LIBRARY_PATH:/cm/shared/apps/boost/1.60.0/lib
export CPATH=$CPATH:/cm/shared/apps/boost/1.60.0/include
eval `scram runtime -sh`

echo $CMSSW_VERSION

todaysdate=$1
massregion=$2
hypo=$3
data=$4

if [[ "$massregion" == "Offshell" ]]; then
   outdirname=$todaysdate"_Offshell_"$hypo"_"$data
   outcardsname="cards_"$outdirname
   datadirname="CMSdata/"$data"/"$hypo"/Offshell/"
   if [ ! -d $outcardsname ];then
      tpldir="templates2D/"$data"/"$hypo
      if [[ "$hypo" == "SM" ]];then
         python makeWidthDatacards.py -b --writeoutput --dopdfproj -i Offshell_inputs_"$data" -t $tpldir -r $datadirname -a $outdirname --coord "mass:KD1:KD2" --GHmodel 1 --CatScheme 2 --ac 0 --mLow 220 --mHigh 13000
      else
         python makeWidthDatacards.py -b --writeoutput --dopdfproj -i Offshell_inputs_"$data" -t $tpldir -r $datadirname -a $outdirname --coord "mass:KD1:KD2" --GHmodel 1 --CatScheme 2 --ac 2 --mLow 220 --mHigh 13000
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
         python makeWidthDatacards.py -b --writeoutput -i Onshell_inputs_SM_"$data" --extMassShapes=$extmsdir -t $tpldir -r $datadirname -a $outdirname --coord "mass:KD1" --GHmodel 0 --CatScheme 2 --ac 0 --mLow 105 --mHigh 140
      else
         python makeWidthDatacards.py -b --writeoutput --dopdfproj -i Onshell_inputs_AC_"$data" -t $tpldir -r $datadirname -a $outdirname --coord "KD1:KD2:KD3" --GHmodel 0 --CatScheme 2 --ac 2 --mLow 105 --mHigh 140
      fi

      mkdir $outcardsname"/HCG/Scans/"

      cp utils/submitScan1D.sh $outcardsname"/HCG/Scans/"
      cp utils/scan1D.slurm.sh $outcardsname"/HCG/Scans/"

      chmod -R 755 "$outcardsname"
   fi
fi
