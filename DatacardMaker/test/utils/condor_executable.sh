#!/bin/bash

getvarpaths(){
  for var in "$@";do
    tmppath=${var//:/ }
    for p in $(echo $tmppath);do
      if [[ -e $p ]];then
        echo $p
      fi
    done
  done  
}
searchfileinvar(){
  for d in $(getvarpaths $1);do
    for f in $(ls $d | grep $2);do
      echo "$d/$f"
    done
  done
}
getcmssw(){
  if [ -r "$OSGVO_CMSSW_Path"/cmsset_default.sh ]; then
    echo "sourcing environment: source $OSGVO_CMSSW_Path/cmsset_default.sh"
    source "$OSGVO_CMSSW_Path"/cmsset_default.sh
  elif [ -r "$OSG_APP"/cmssoft/cms/cmsset_default.sh ]; then
    echo "sourcing environment: source $OSG_APP/cmssoft/cms/cmsset_default.sh"
    source "$OSG_APP"/cmssoft/cms/cmsset_default.sh
  elif [ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
    echo "sourcing environment: source /cvmfs/cms.cern.ch/cmsset_default.sh"
    source /cvmfs/cms.cern.ch/cmsset_default.sh
  else
    echo "ERROR! Couldn't find $OSGVO_CMSSW_Path/cmsset_default.sh or /cvmfs/cms.cern.ch/cmsset_default.sh or $OSG_APP/cmssoft/cms/cmsset_default.sh"
    exit 1
  fi
}


CMSSWVERSION="$1"
SCRAMARCH="$2"
TARFILE=$3
NCPUS=$4
CONDORSITE=$5
CONDOROUTDIR=$6
NPOINTS=$7
FIRSTPOINT=$8
LASTPOINT=$9
OPTS=${10}

WSNAME=workspace.root
OUTNAME="Scan"

TOYSFILE=toysfile.root
if [[ ! -f ${TOYSFILE} ]]; then
  TOYSFILE=""
fi


export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "WSNAME: $WSNAME"
echo "TOYSFILE: $TOYSFILE"
echo "OUTNAME: $OUTNAME"
echo "CONDORSITE: $CONDORSITE"
echo "CONDOROUTDIR: $CONDOROUTDIR"
echo "NPOINTS: $NPOINTS"
echo "FIRSTPOINT: $FIRSTPOINT"
echo "LASTPOINT: $LASTPOINT"
echo "OPTS: $OPTS"


echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "whoami: $(whoami)"
echo "time: $(date +%s)"
echo "args: $@"
echo -e "\n--- end header output ---\n" #                       <----- section division

echo -e "\n--- begin memory specifications ---\n" #                     <----- section division
ulimit -s unlimited
ulimit -a
echo -e "\n--- end memory specifications ---\n" #                     <----- section division

INITIALDIR=$(pwd)
echo "Initial directory is ${INITIALDIR}"

mkdir -p rundir
cd rundir

getcmssw


# If the first file in the tarball filelist starts with CMSSW, it is a
# tarball made outside of the full CMSSW directory and must be handled
# differently
if [[ ! -s ${INITIALDIR}/${TARFILE} ]];then
  echo "Tar file ${INITIALDIR}/${TARFILE} does not exist"
elif [[ ! -z $(tar -tf ${INITIALDIR}/${TARFILE} | head -n 1 | grep "^CMSSW") ]]; then

  echo "This is a full cmssw tar file."

  mv ${INITIALDIR}/${TARFILE} $(pwd)/
  if [[ "${TARFILE}" == *".tgz" ]];then
    tar zxf ${TARFILE}
  else
    tar xf ${TARFILE}
  fi
  rm ${TARFILE}

  if [[ -e extras.more ]];then
    mv extras.more ${CMSSWVERSION}/extras.tar
  fi

  cd $CMSSWVERSION

  if [[ -e extras.more ]];then
    tar xf extras.tar
    rm extras.tar
  fi

  echo "Current directory ${PWD} =? ${CMSSWVERSION}"
  echo "Running ProjectRename"
  scramv1 b ProjectRename

else

  # Setup the CMSSW area
  echo "This is a selective CMSSW tar file."
  eval `scramv1 project CMSSW $CMSSWVERSION`
  cd $CMSSWVERSION

  mv ${INITIALDIR}/${TARFILE} $(pwd)/
  if [[ "${TARFILE}" == *".tgz" ]];then
    tar zxf ${TARFILE}
  else
    tar xf ${TARFILE}
  fi
  rm ${TARFILE}

  if [[ -e extras.more ]];then
    mv extras.more extras.tar
    tar xf extras.tar
    rm extras.tar
  fi

fi


# Setup the CMSSW environment
eval `scramv1 runtime -sh`
echo "CMSSW_BASE: ${CMSSW_BASE}"

# Make sure JHUGenMELA library location is recognized
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CMSSW_BASE}/src/JHUGenMELA/MELA/data/${SCRAM_ARCH}

cd src
eval $(./IvyFramework/IvyDataTools/setup.sh env)
scramv1 b -j 1
CMSSW_COMPILE_STATUS=$?
if [ $CMSSW_COMPILE_STATUS != 0 ];then
  echo "CMSSW compilation exited with error ${CMSSW_COMPILE_STATUS}. Printing the log:"
  cat compilation.log
  exit ${CMSSW_COMPILE_STATUS}
fi
rm -f compilation.log


mkdir -p subdir
cd subdir

mv ${INITIALDIR}/${WSNAME} ./


cmdadd=""
if [ $FIRSTPOINT -ge 0 ]; then
  cmdadd=$cmdadd" --firstPoint "$FIRSTPOINT
  OUTNAME=$OUTNAME"_"$FIRSTPOINT
fi
if [ $LASTPOINT -ge 0 ]; then
  cmdadd=$cmdadd" --lastPoint "$LASTPOINT
  OUTNAME=$OUTNAME"_"$LASTPOINT
fi
OUTNAME=$OUTNAME"_"$NPOINTS

useObs=0
if [[ "$OPTS" == *".Observed"* ]]; then
  useObs=1
  OPTS=${OPTS//.Observed}
fi

# Use the first one if you are starting from an initial fit with an embedded workspace
# and the second if you are also doing the initial fit in the same job.
if [[ $useObs -eq 0 ]]; then
  #cmdadd="$cmdadd --snapshotName MultiDimFit --skipInitialFit -D toys/toy_asimov"
  cmdadd="$cmdadd -t -1 --saveToys"
else
  #cmdadd="$cmdadd --snapshotName MultiDimFit --skipInitialFit --cminPreFit 1"
  cmdadd="$cmdadd --cminPreFit 1"
fi
if [[ $NPOINTS -eq 0 ]]; then
  cmdadd="$cmdadd --saveWorkspace --saveFitResult"
fi
if [[ "$OPTS" == *".Frequentist"* ]]; then
  OPTS=${OPTS//.Frequentist}
  cmdadd="${cmdadd} --toysFrequentist"
fi
if [[ "$OPTS" == *".LoadSnapshot"* ]]; then
  OPTS=${OPTS//.LoadSnapshot}
  cmdadd="${cmdadd} --snapshotName MultiDimFit"
fi

FITOPTS="-M MultiDimFit --algo grid --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd MINIMIZER_no_analytic -m 125 --alignEdges=1 --saveNLL --saveSpecifiedNuis=all --saveSpecifiedFunc=R,RV,RF,R_13TeV,RV_13TeV,RF_13TeV -v 3 --points ${NPOINTS}"
if [[ ${TOYSFILE} != "" ]]; then
  mv ${INITIALDIR}/${TOYSFILE} ./
  FITOPTS="${FITOPTS} --toysFile=${TOYSFILE}"
  if [[ ${useObs} -eq 1 ]]; then
    echo "Using a toys file with an observed fit is not possible."
    exit 1
  fi
fi

# If this option is not empty, make sure you precede with a comma as with any other parameter passed to --setParameters
LUMIOPTS=",LUMI_13TeV_2016=36.326450"

GGsm_range="0,6"
R_range="0,6"
RF_range="0,6"
RV_range="0,10"
fai_range=""
if [[ "$OPTS" == "a2."* ]] || [[ "$OPTS" == "a3."* ]]; then
  fai_range="-0.02,0.02"
elif [[ "$OPTS" == "L1."* ]]; then
  fai_range="-0.004,0.004"
fi

if [[ "$OPTS" == "SM.RF_RV" ]]; then
  RF_range="0,4"
  RV_range="0,4"
fi

rvoffshell_range="0,400"
if [[ "$OPTS" == "SM.R_RVFloated" ]]; then
  if [[ ${NPOINTS} -gt 1000 ]] && [[ ${NPOINTS} -lt 200000 ]]; then
    rvoffshell_range="0,50000"
  elif [[ ${NPOINTS} -ge 200000 ]] && [[ ${NPOINTS} -lt 2000000 ]]; then
    rvoffshell_range="0,500000"
  fi
fi

POIOPTS=""
if [[ "$OPTS" == "SM.GGsm" ]]; then
  POIOPTS="--redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm=${GGsm_range} --setParameters GGsm=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == "SM.RF" ]]; then
  POIOPTS="--redefineSignalPOIs=rf_offshell --freezeParameters=GGsm,R,RF,RV,r_offshell,CMS_zz4l_fai1,kbkg_VBF --setParameterRanges rf_offshell=${RF_range} --setParameters r_offshell=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == "SM.RV" ]]; then
  POIOPTS="--redefineSignalPOIs=rv_offshell --freezeParameters=GGsm,R,RF,RV,r_offshell,CMS_zz4l_fai1,kbkg_VBF --setParameterRanges rv_offshell=${RV_range} --setParameters r_offshell=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == "SM.RF_RV" ]]; then
  POIOPTS="--redefineSignalPOIs=rf_offshell,rv_offshell --freezeParameters=GGsm,R,RF,RV,r_offshell,CMS_zz4l_fai1,kbkg_VBF --setParameterRanges rf_offshell=${RF_range}:rv_offshell=${RV_range} --setParameters r_offshell=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == "SM.R_RVFloated" ]]; then
  POIOPTS="--redefineSignalPOIs=r_offshell --freezeParameters=GGsm,R,RF,RV,CMS_zz4l_fai1,kbkg_VBF --setParameterRanges r_offshell=${R_range}:rv_offshell=${rvoffshell_range} --setParameters r_offshell=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == "SM.R_RVFixed" ]]; then
  POIOPTS="--redefineSignalPOIs=r_offshell --freezeParameters=GGsm,R,RF,RV,rv_offshell,CMS_zz4l_fai1,kbkg_VBF --setParameterRanges r_offshell=${R_range} --setParameters r_offshell=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == *".GGsm_fai1Floated" ]]; then
  POIOPTS="--redefineSignalPOIs=GGsm --freezeParameters=kbkg_VBF --setParameterRanges GGsm=${GGsm_range}:CMS_zz4l_fai1=${fai_range} --setParameters GGsm=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == *".GGsm_fai1Fixed" ]]; then
  POIOPTS="--redefineSignalPOIs=GGsm --freezeParameters=CMS_zz4l_fai1,kbkg_VBF --setParameterRanges GGsm=${GGsm_range} --setParameters GGsm=1,CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == *".fai1_Onshell" ]]; then
  POIOPTS="--redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=kbkg_VBF --setParameterRanges CMS_zz4l_fai1=${fai_range} --setParameters CMS_zz4l_fai1=0${LUMIOPTS}"
elif [[ "$OPTS" == *".fai1_GGsmFloated" ]]; then
  POIOPTS="--redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=kbkg_VBF --setParameterRanges CMS_zz4l_fai1=${fai_range}:GGsm=${GGsm_range} --setParameters CMS_zz4l_fai1=0,GGsm=1${LUMIOPTS}"
elif [[ "$OPTS" == *".fai1_GGsmFixed" ]]; then
  POIOPTS="--redefineSignalPOIs=CMS_zz4l_fai1 --freezeParameters=GGsm,kbkg_VBF --setParameterRanges CMS_zz4l_fai1=${fai_range} --setParameters CMS_zz4l_fai1=0,GGsm=1${LUMIOPTS}"
fi

if [[ "${POIOPTS}" == "" ]]; then
  echo "No POI options for the run options ${OPTS}."
  exit 1
fi

cmd="${WSNAME} ${FITOPTS} ${POIOPTS} ${cmdadd} -n ${OUTNAME}"


RUNDIR=$(pwd)
echo "Running: combine ${cmd}"
combine ${cmd}

RUN_STATUS=$?
if [[ ${RUN_STATUS} -ne 0 ]]; then
  echo "Run has crashed with exit code ${RUN_STATUS}"
  exit 1
fi

if [[ -f higgsCombine${OUTNAME}.MultiDimFit.mH125.123456.root ]]; then
  mv higgsCombine${OUTNAME}.MultiDimFit.mH125.123456.root higgsCombine${OUTNAME}.MultiDimFit.mH125.root
fi
for f in $(find ./ -name *multidimfit*); do
  mv $f higgsCombine${OUTNAME}.MultiDimFit.mH125.FitResult.root
done

touch EXTERNAL_TRANSFER_LIST.LST
echo higgsCombine${OUTNAME}.MultiDimFit.mH125.root >> EXTERNAL_TRANSFER_LIST.LST
if [[ -e higgsCombine${OUTNAME}.MultiDimFit.mH125.FitResult.root ]]; then
  echo higgsCombine${OUTNAME}.MultiDimFit.mH125.FitResult.root >> EXTERNAL_TRANSFER_LIST.LST
fi


echo "Submission directory after running all steps before file transfers: ls -lrth"
ls -lrth

##################
# TRANSFER FILES #
##################
# In cases where transfer through the script fails
if [[ -f "EXTERNAL_TRANSFER_LIST.LST" ]];then
  echo -e "\n--- Begin EXTERNAL TRANSFER ---\n"
  while IFS='' read -r line || [[ -n "$line" ]]; do
    OUTFILENAME=${line}
    # If there is an instruction to compress, convert the file/directory name into a tar file.
    if [[ "${OUTFILENAME}" == "compress:"* ]]; then
      OUTFILENAME=${OUTFILENAME/'compress:'}
      if [[ "${OUTFILENAME}" == *"/" ]]; then
        OUTFILENAME=${OUTFILENAME%?}
      fi
      tar Jcf ${OUTFILENAME}.tar ${OUTFILENAME}
      OUTFILENAME=${OUTFILENAME}.tar
    fi
    # Begin copying the file
    echo "Copying output file ${OUTFILENAME}"
    copyFromCondorToSite.sh ${RUNDIR} ${OUTFILENAME} ${CONDORSITE} ${CONDOROUTDIR}
    TRANSFER_STATUS=$?
    if [ $TRANSFER_STATUS != 0 ]; then
      echo " - Transfer crashed with exit code ${TRANSFER_STATUS}"
    fi
  done < "EXTERNAL_TRANSFER_LIST.LST"
  echo -e "\n--- End EXTERNAL TRANSFER ---\n"
fi
##############


echo "Submission directory after running all steps: ls -lrth"
ls -lrth

echo "time at end: $(date +%s)"
