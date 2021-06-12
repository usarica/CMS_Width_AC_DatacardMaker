#!/bin/env python

import sys
import imp
import copy
import os
import filecmp
import shutil
import pickle
import math
import pprint
import subprocess
import socket
from datetime import date
from optparse import OptionParser
from CMSDataTools.AnalysisTree.TranslateStringBetweenPythonAndShell import *


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--batchqueue", type="string", help="Batch queue")
      self.parser.add_option("--batchscript", type="string", help="Name of the HTCondor script")
      self.parser.add_option("--tarfile", type="string", help="Name of the tar file to upload")
      self.parser.add_option("--wsfile", type="string", help="Name of the workspace file to upload")
      self.parser.add_option("--outdir", type="string", help="Name of the output directory")

      self.parser.add_option("--condorsite", type="string", help="Name of the HTCondor site")
      self.parser.add_option("--condoroutdir", type="string", help="Name of the HTCondor output directory")

      self.parser.add_option("--job_arg", type="string", action="append", help="Job arguments. Could be npoints, firstpoint, lastpoint etc. Must respect the order.", default=None)

      self.parser.add_option("--outlog", type="string", help="Name of the output log file")
      self.parser.add_option("--errlog", type="string", help="Name of the output error file")

      self.parser.add_option("--required_memory", type="string", default="2048M", help="Required RAM for the job")
      self.parser.add_option("--required_ncpus", type="int", default=1, help="Required number of CPUs for the job")
      self.parser.add_option("--job_flavor", type="string", default="tomorrow", help="Time limit for job (tomorrow = 1 day, workday = 8 hours, see https://batchdocs.web.cern.ch/local/submit.html#job-flavours for more)")

      self.parser.add_option("--extra_upload", type="string", action="append", help="Extra uploads in the format [link_name]:[file_name]", default=[])

      self.parser.add_option("--use_cloud", action="store_true", default=False, help="Use cloud computing")
      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")


      (self.opt,self.args) = self.parser.parse_args()
      optchecks=[
         "batchqueue",
         "batchscript",
         "job_arg",
         "tarfile",
         "wsfile",
         "outdir",
         "condorsite",
         "condoroutdir",
         "outlog",
         "errlog"
      ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      if self.opt.outdir.startswith("./"):
         self.opt.outdir = self.opt.outdir.replace(".",os.getcwd(),1)

      if not os.path.isfile(self.opt.batchscript):
         print("Batch script does not exist in current directory, will search for CMSSW_BASE/bin")
         if os.path.isfile(os.getenv("CMSSW_BASE")+"/bin/"+os.getenv("SCRAM_ARCH")+"/"+self.opt.batchscript):
            self.opt.batchscript = os.getenv("CMSSW_BASE")+"/bin/"+os.getenv("SCRAM_ARCH")+"/"+self.opt.batchscript
            print("\t- Found the batch script")
         else:
            sys.exit("Batch script {} does not exist. Exiting...".format(self.opt.batchscript))

      for theOpt in optchecks:
         print("Option {}={}".format(theOpt,getattr(self.opt, theOpt)))

      self.submitJobs()


   def produceCondorScript(self):
      currentdir = os.getcwd()
      currentCMSSWBASESRC = os.getenv("CMSSW_BASE")+"/src/" # Need the trailing '/'
      currendir_noCMSSWsrc = currentdir.replace(currentCMSSWBASESRC,'')

      hostname = socket.gethostname()

      scramver = os.getenv("SCRAM_ARCH")
      singularityver = "cms:rhel6"
      if "slc7" in scramver:
         singularityver = "cms:rhel7"

      os.system("ln -sf {} {}/workspace.root".format(self.opt.wsfile, self.opt.outdir))

      uploads=[]
      uploads.append(self.opt.tarfile)
      uploads.append("workspace.root")
      for extra_dict in self.opt.extra_upload:
         tmpvals = extra_dict.split(':')
         if len(tmpvals) != 2:
            raise RuntimeError("Incorrect format for extra upload string '{}'.".format(extra_dict))
         elif not os.path.isfile(tmpvals[1]):
            raise RuntimeError("Extra upload file {} does not exist.".format(tmpvals[1]))
         uploads.append(tmpvals[0])
         os.system("ln -sf {} {}/{}".format(os.path.abspath(tmpvals[1]), self.opt.outdir, tmpvals[0]))

      allowed_sites = None
      cloud_arg = ""
      if self.opt.use_cloud:
         allowed_sites = "T2_US_UCSD"
         cloud_arg = "+IS_CLOUD_JOB=True"
      else:
         allowed_sites = "T2_US_UCSD,T2_US_Caltech,T2_US_MIT,T3_US_UCR,T3_US_Baylor,T3_US_Colorado,T3_US_NotreDame,T3_US_Cornell,T3_US_Rice,T3_US_Rutgers,T3_US_UCD,T3_US_TAMU,T3_US_TTU,T3_US_FIU,T3_US_FIT,T3_US_UMD,T3_US_OSU,T3_US_OSG,T3_US_UMiss,T3_US_PuertoRico"

      strrequirements = r'(HAS_SINGULARITY=?=True) && !(( regexp("(mh-epyc7662-1)\..*",TARGET.Machine) || regexp("(mh-epyc7662-5)\..*",TARGET.Machine) || regexp("(mh-epyc7662-6)\..*",TARGET.Machine) || regexp("(mh-epyc7662-9)\..*",TARGET.Machine) || regexp("(mh-epyc7662-10)\..*",TARGET.Machine) || regexp("(sdsc-84)\..*",TARGET.Machine) || regexp("(sdsc-3)\..*",TARGET.Machine) || regexp("(cabinet-0-0-29)\..*",TARGET.Machine) || regexp("(cabinet-0-0-23)\..*",TARGET.Machine) || regexp("(cabinet-0-0-21)\..*",TARGET.Machine) || regexp("(cabinet-11-11-3)\..*",TARGET.Machine) )=?=True)'
      if "uscms.org" in hostname:
         strrequirements = r'(HAS_SINGULARITY=?=True) || (NODE_MOUNTS_CVMFS =?= true)'


      strjobargs=' '.join(self.opt.job_arg)

      scriptargs = {
         "home" : os.path.expanduser("~"),
         "uid" : os.getuid(),
         "batchScript" : self.opt.batchscript,
         "outDir" : self.opt.outdir,
         "outLog" : self.opt.outlog,
         "errLog" : self.opt.errlog,
         "QUEUE" : self.opt.batchqueue,
         "SINGULARITYVERSION" : singularityver,
         "REQMEM" : self.opt.required_memory,
         "NCPUS" : self.opt.required_ncpus,
         "JOBFLAVOR" : self.opt.job_flavor,
         "UPLOADS" : ",".join(uploads),
         "CMSSWVERSION" : os.getenv("CMSSW_VERSION"),
         "SCRAMARCH" : scramver,
         "SUBMITDIR" : currendir_noCMSSWsrc,
         "TARFILE" : self.opt.tarfile,
         "JOBARGS" : strjobargs,
         "DESIRED_SITES" : allowed_sites,
         "CLOUD_ARG" : cloud_arg,
         "CONDORSITE" : self.opt.condorsite,
         "CONDOROUTDIR" : self.opt.condoroutdir,
         "REQUIREMENTS" : strrequirements
      }

      scriptcontents = """
universe={QUEUE}
+DESIRED_Sites="{DESIRED_SITES}"
executable              = {batchScript}
arguments               = {CMSSWVERSION} {SCRAMARCH} {TARFILE} {NCPUS} {CONDORSITE} {CONDOROUTDIR} {JOBARGS}
Initialdir              = {outDir}
output                  = {outLog}.$(ClusterId).$(ProcId).txt
error                   = {errLog}.$(ClusterId).$(ProcId).err
log                     = $(ClusterId).$(ProcId).log
request_memory          = {REQMEM}
request_cpus            = {NCPUS}
+JobFlavour             = "{JOBFLAVOR}"
x509userproxy           = {home}/x509up_u{uid}
#https://www-auth.cs.wisc.edu/lists/htcondor-users/2010-September/msg00009.shtml
periodic_remove         = JobStatus == 5
transfer_executable=True
transfer_input_files    = {UPLOADS}
transfer_output_files = ""
+Owner = undefined
+project_Name = "cmssurfandturf"
notification=Never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
Requirements = {REQUIREMENTS}
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/{SINGULARITYVERSION}"
{CLOUD_ARG}


queue

"""
      scriptcontents = scriptcontents.format(**scriptargs)

      self.condorScriptName = "condor.sub"
      condorScriptFile = open(self.opt.outdir+"/"+self.condorScriptName,'w')
      condorScriptFile.write(scriptcontents)
      condorScriptFile.close()


   def submitJobs(self):
      self.produceCondorScript()

      jobcmd = "cd {}; condor_submit {}; cd -".format(self.opt.outdir, self.condorScriptName)
      if self.opt.dryRun:
         print("Job command: '{}'".format(jobcmd))
      else:
         ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = BatchManager()
