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

wname=$1
outname=$2
poi=$3
rangel=$4
rangeh=$5
let npoints=$6
let firstpoint=$7
let lastpoint=$8
extarg=$9

cmdadd=""
if [ $firstpoint -ge 0 ];then
  cmdadd=$cmdadd" --firstPoint "$firstpoint
  outname=$outname"_"$firstpoint
fi

if [ $lastpoint -ge 0 ];then
  cmdadd=$cmdadd" --lastPoint "$lastpoint
  outname=$outname"_"$lastpoint
fi

if [[ "$poi" == "GGsm" ]];then
  echo "POI is GGsm"
  cmdadd=$cmdadd" --redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm="$rangel","$rangeh
elif [[ "$poi" == "fai1" ]];then
  echo "POI is fai1"
  cmdadd=$cmdadd" --redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=GGsm,kbkg_VBF --setParameterRanges CMS_zz4l_fai1="$rangel","$rangeh
fi


cmd="-M MultiDimFit "$wname" --algo=grid --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --alignEdges=1 -m 125 --saveNLL --saveSpecifiedNuis=all --saveSpecifiedFunc=R,RV,RF,R_13TeV,RV_13TeV,RF_13TeV -v 3 -S 1 -t -1 --points "$npoints" -n "$outname" "$cmdadd $extarg

echo "Command: combine "$cmd
combine $cmd
