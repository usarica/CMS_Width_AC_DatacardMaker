#!/bin/env python

import os
import pprint
import argparse
import subprocess
import math
import multiprocessing as mp


def makeDirectory(subDirName):
   ret = os.system("mkdir -p {} > /dev/null".format(subDirName))
   if ret != 0:
      raise RuntimeError("makeDirectory: Error in creating {}.".format(subDirName))


def getVOMSProxy():
   gridproxy = None
   if os.getenv("X509_USER_PROXY") is None or not os.getenv("X509_USER_PROXY"):
      gridproxycheckfiles = [
         "{home}/x509up_u{uid}".format(home=os.path.expanduser("~"), uid=os.getuid()),
         "/tmp/x509up_u{uid}".format(uid=os.getuid())
         ]
      for gridproxycheckfile in gridproxycheckfiles:
         if os.path.exists(gridproxycheckfile):
            gridproxy = gridproxycheckfile
            break
   else:
      gridproxy = os.getenv("X509_USER_PROXY")
   if gridproxy is None or not os.path.exists(gridproxy):
      sys.exit("Cannot find a valid grid proxy")
   return gridproxy


def run_single(args,grid_user,jobmaindir,batchscript,condorsite,run_arg,run_mode):
   wsfile = args.wsfile
   run_args = []
   extra_uploads = args.extra_upload
   condoroutdir="/hadoop/cms/store/user/{}/Offshell_2L2Nu/Worker/output/WidthACFits/{}".format(grid_user,args.date)
   jobdir = None

   if run_mode.lower() == "nll":
      npoints = args.npoints
      ipoint = int(run_arg)
      jobdir = jobmaindir+"/job_{}_{}".format(ipoint,npoints)
      run_args = [ npoints, ipoint, ipoint ]
   elif run_mode.lower() == "impact_fits":
      jobdir = jobmaindir+"/job_{}".format(run_arg)
      run_args = [ run_arg ]
   elif run_mode.lower() == "toygen":
      jobdir = jobmaindir+"/job_{}".format(run_arg)
      run_args = [ run_arg ]
   elif run_mode.lower() == "gof":
      jobdir = jobmaindir+"/job_{}".format(run_arg[0])
      run_args = [ run_arg[0] ]
      extra_uploads.append("toysfile.root:{}".format(run_arg[1]))
   elif run_mode.lower() == "singlefit":
      jobdir = jobmaindir+"/job_{}".format(run_arg[0])
      run_args = [ run_arg[0] ]
      extra_uploads.append("toysfile.root:{}".format(run_arg[1]))


   if len(run_args) == 0:
      raise RuntimeError("Job arguments list is empty.")

   # Add extra arguments for the executable if there are any.
   for extra_arg in args.extra_args:
      run_args.append(extra_arg)

   print("Creating {}".format(jobdir))
   makeDirectory(jobdir+"/Logs")
   ret = os.system("ln -sf {} {}/".format(batchscript, jobdir))
   if ret != 0:
      print("Could not link the batch script in {}".format(jobdir))

   ret = os.system("ln -sf {}/cmswidthac.tar {}/".format(jobmaindir, jobdir))
   if ret != 0:
      print("Could not link the main tar file in {}".format(jobdir))

   jobargs = {
      "batchqueue" : "vanilla",
      "batchscript" : batchscript,
      "tarfile" : "cmswidthac.tar",
      "wsfile" : wsfile,
      "outdir" : jobdir,
      "condorsite" : condorsite,
      "condoroutdir" : condoroutdir,
      "outlog" : "Logs/log_job",
      "errlog" : "Logs/err_job",
      "required_memory" : "4096M",
      "job_flavor" : "tomorrow"
   }
   runCmd = "configureCMSWidthACCondorJob.py"
   if not args.direct_submit:
      runCmd = runCmd + " --dry"
   for key,val in jobargs.iteritems():
      runCmd = runCmd + " --{}={}".format(key,val)
   for fdict in extra_uploads:
      runCmd = runCmd + " --extra_upload={}".format(fdict)
   for run_arg in run_args:
      runCmd = runCmd + " --job_arg={}".format(run_arg)
   if args.use_cloud:
      runCmd = runCmd + " --use_cloud --required_ncpus=1 --required_memory=2048M" # Ignore memory specification
   else:
      runCmd = runCmd + " --required_memory={}".format(args.required_memory)

   if args.dry:
      print("Job command: '{}'".format(runCmd))
   else:
      os.system(runCmd + " > /dev/null")
      #os.system(runCmd)
      if ret != 0:
         raise RuntimeError("Could not run command '{}'".format(runCmd))


def getLeftRightAlignedPoint(n,N,c):
   return int(-float(c)*float(N)*math.log(1.-float(n)/float(N)*(1.-math.exp(-1./float(c)))))

def run(args):
   CMSSWBASE = os.getenv("CMSSW_BASE")
   SCRAMARCH = os.getenv("SCRAM_ARCH")
   batchscript = os.path.abspath(args.executable)
   args.wsfile = os.path.abspath(args.wsfile)
   condorsite="t2.ucsd.edu"
   nthreads = min(args.nthreads, mp.cpu_count())

   gridproxy = getVOMSProxy()
   grid_user = subprocess.check_output("voms-proxy-info -identity -file={} | cut -d '/' -f6 | cut -d '=' -f2".format(gridproxy), shell=True)
   if not grid_user:
      grid_user = os.environ.get("USER")
   grid_user = grid_user.strip()

   curdir=os.getcwd()
   jobmaindir=curdir+"/tasks/{}".format(args.date)
   makeDirectory(jobmaindir)
   if not os.path.exists("{}/cmswidthac.tar".format(jobmaindir)) or args.remake:
      if args.precompiled_tar is not None and args.precompiled_tar!="" and os.path.exists(args.precompiled_tar):
         os.system("cp {} {}/cmswidthac.tar".format(args.precompiled_tar, jobmaindir))
      else:
         os.system("cd {}; createCMSWidthACDatacardTarball.sh; cd -;".format(jobmaindir))

   use_likelihood = (args.npoints is not None)
   use_impacts = (args.impact_parsfile is not None)
   use_gentoys = (args.generate_ntoys is not None)
   use_GoF = args.run_GoF
   use_SingleFit = args.run_SingleFit

   run_args = []
   run_mode = None
   if use_likelihood:
      run_mode = "nll"
      if args.point_distribution == 'uniform':
         for ipoint in range(0, args.npoints):
            if (args.firstpoint is not None and ipoint<args.firstpoint) or (args.lastpoint is not None and ipoint>args.lastpoint):
               continue
            run_args.append(ipoint)
      elif args.point_distribution == 'left-aligned':
         npoints_original = args.npoints
         isOdd = False
         npoints_mult = int(100)
         npoints_c = 1./3.
         if args.npoints % 2 == 1:
            isOdd = True
            args.npoints = args.npoints - 1
         args.npoints = args.npoints * npoints_mult
         for ipoint in range(0, args.npoints):
            if ipoint % npoints_mult != 0:
               continue
            run_args.append(getLeftRightAlignedPoint(ipoint,args.npoints,npoints_c))
         if isOdd:
            run_args.append(args.npoints)
            args.npoints = args.npoints+1
      elif args.point_distribution == 'right-aligned':
         npoints_original = args.npoints
         isOdd = False
         npoints_mult = int(100)
         npoints_c = -1./3.
         if args.npoints % 2 == 1:
            isOdd = True
            args.npoints = args.npoints - 1
         args.npoints = args.npoints * npoints_mult
         for ipoint in range(0, args.npoints):
            if ipoint % npoints_mult != 0:
               continue
            run_args.append(getLeftRightAlignedPoint(ipoint,args.npoints,npoints_c))
         if isOdd:
            run_args.append(args.npoints)
            args.npoints = args.npoints+1
      elif args.point_distribution == 'centered':
         npoints_original = args.npoints
         isOdd = False
         npoints_mult = int(100)
         npoints_c = -1./3.
         if args.npoints % 2 == 1:
            isOdd = True
            args.npoints = args.npoints - 1
         args.npoints = args.npoints * npoints_mult
         for ipoint in range(0, args.npoints/2):
            if ipoint % npoints_mult != 0:
               continue
            run_args.append(getLeftRightAlignedPoint(ipoint,args.npoints/2,npoints_c))
            run_args.append(getLeftRightAlignedPoint(ipoint,args.npoints/2,-npoints_c)+args.npoints/2)
         if isOdd:
            run_args.append(args.npoints)
            args.npoints = args.npoints+1

   elif use_impacts:
      run_mode = "impact_fits"
      with open(args.impact_parsfile,'rb') as fin:
         for par in fin:
            par = par.strip()
            if par != "":
               run_args.append(par)

   elif use_gentoys:
      run_mode = "toygen"
      seed_begin = 100
      for itoy in range(0,args.generate_ntoys):
         run_args.append(seed_begin + itoy)

   elif use_GoF:
      run_mode = "gof"
      itoy=0
      for fname in os.listdir(args.toysdir):
         if ".root" in fname:
            run_args.append([itoy, args.toysdir + '/' + fname])
            itoy = itoy+1

   elif use_SingleFit:
      # Workflow is the same as a GoF test, but keep the case separate in case it needs to specialize later on
      run_mode = "singlefit"
      itoy=0
      for fname in os.listdir(args.toysdir):
         if ".root" in fname:
            run_args.append([itoy, args.toysdir + '/' + fname])
            itoy = itoy+1

   print("Running in mode {}...".format(run_mode))

   pool = mp.Pool(nthreads)
   #starmap does not exist until Python3...
   #pool.starmap_async(run_single, [(njdir, refcsubs) for njdir, refcsubs in sub_data])
   [ pool.apply_async(run_single, args=(args, grid_user, jobmaindir, batchscript, condorsite, run_arg, run_mode)) for run_arg in run_args ]
   pool.close()
   pool.join()


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--date", type=str, help="Tag for the fits", required=True)
   parser.add_argument("--nthreads", help="Number of threads", type=int, required=False, default=4)
   parser.add_argument("--remake", action="store_true", help="Remake all jobs", required=False)
   parser.add_argument("--skip_existing", action="store_true", help="Skip existing folders in addition to tar files to save time", required=False)
   parser.add_argument("--dry", action="store_true", help="Test run without creation of jobs", required=False)
   parser.add_argument("--direct_submit", action="store_true", help="Submut jobs directly", required=False)
   parser.add_argument("--executable", type=str, help="Executable to run on condor", required=True)
   parser.add_argument("--wsfile", type=str, help="Workspace file", required=True)
   parser.add_argument("--extra_upload", action="append", help="Extra uploads. Must be in the format [link_name]:[file_name]", required=False, default=[])
   parser.add_argument("--extra_args", action="append", help="Extra arguments to the executable. Can be any string.", required=False, default=[])
   parser.add_argument("--npoints", type=int, help="Number of points to submit", required=False, default=None)
   parser.add_argument("--firstpoint", type=int, help="Index of the first point (default=None)", required=False, default=None)
   parser.add_argument("--lastpoint", type=int, help="Index of the last point (default=None)", required=False, default=None)
   parser.add_argument("--point_distribution", type=str, help="Distribution of points, can be 'left-aligned', 'centered', 'right-aligned', or 'uniform' (default)", required=False, default='uniform')
   parser.add_argument("--impact_parsfile", type=str, help="File that lists the parameters to run impacts", required=False, default=None)
   parser.add_argument("--generate_ntoys", type=int, help="Number of toys to generate", required=False, default=None)
   parser.add_argument("--run_GoF", action="store_true", help="Run the GoF tests, requires --toysdir as well.", required=False)
   parser.add_argument("--run_SingleFit", action="store_true", help="Run a single toy fit, requires --toysdir as well.", required=False)
   parser.add_argument("--toysdir", type=str, help="Directory of single toys", required=False, default=None)
   parser.add_argument("--required_memory", type=str, help="Required RAM for the job", required=False, default="2048M")
   parser.add_argument("--precompiled_tar", type=str, help="Precompiled tar file", required=False, default=None)
   parser.add_argument("--use_cloud", action="store_true", help="Use cloud computing submission", required=False)

   args = parser.parse_args()
   if args.npoints is None and args.impact_parsfile is None and args.generate_ntoys is None and not args.run_GoF and not args.run_SingleFit:
      raise RuntimeError("You must specify a likelihood, an impacts, or a toy generation run, or run the GoF tests/single fits.")
   elif args.npoints is not None and args.npoints < 0:
      raise RuntimeError("You must specify npoints>=0 for the likelihood run.")
   elif args.generate_ntoys is not None and args.generate_ntoys <= 0:
      raise RuntimeError("Number of toys to generate should be greater than 0.")
   elif (args.run_GoF or args.run_SingleFit) and args.toysdir is None:
      raise RuntimeError("GoF and single fit runs need a toys directory.")

   if (args.firstpoint is not None or args.lastpoint is not None) and args.point_distribution.lower() != 'uniform':
      raise RuntimeError("First and last points are only implemented for the uniform point distribution at the moment.")

   args.point_distribution = args.point_distribution.lower()

   if args.dry:
      print("Dry mode is specified.")

   run(args)
