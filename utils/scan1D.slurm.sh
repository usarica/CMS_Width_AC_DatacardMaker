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

source scan1D.run.sh $@
