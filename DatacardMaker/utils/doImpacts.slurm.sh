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

module load gcc/6.4.0

source /work-zfs/lhc/cms/cmsset_default.sh

eval `scram runtime -sh`

echo $CMSSW_VERSION

source doImpacts.run.sh $@
