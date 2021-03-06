#!/usr/bin/python

import pwd
import commands
import optparse
import shlex
from CMS_Width_AC_DatacardMaker.DatacardMaker.WidthDatacardMaker import *


def loadIncludes():
   # ---------------- LOAD THE INCLUDES ----------------- ##
   ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
   ROOT.gSystem.AddIncludePath("-Iutils/")
   ROOT.gROOT.ProcessLine(".L utils/tdrstyle.cc")
   ROOT.gSystem.Load("libRooFit")
   ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
   # ---------------- SET PLOTTING STYLE ---------------- ##
   ROOT.setTDRStyle(True)
   ROOT.gStyle.SetPalette(1)
   ROOT.gStyle.SetPadLeftMargin(0.16)

class makeWidthDatacards:
   def __init__(self):
      # parse the arguments and options
      self.parseOptions()

      loadIncludes()

      dirName = 'cards_' + self.opt.appendName
      subdir = ['HCG', 'figs']
      for d in subdir:
         self.makeDirectory(dirName + '/' + d)
      self.creationLoop(dirName)

      sys.exit()


   def parseOptions(self):
      parser = optparse.OptionParser()

      parser.add_option('-b', action='store_true',
                         dest='noX', default=True, help='no X11 windows')

      parser.add_option('-i', '--input', dest='inputDir',
                         type='string', default="", help='Inputs directory')
      parser.add_option('-t', '--templates', dest='templateDir',
                         type='string', default="", help='Directory of templates')
      parser.add_option('--extMassShapes', type='string', default=None, help='Directory of external mass shapes')
      parser.add_option('-a', '--append', dest='appendName',
                         type='string', default="", help='Append name for cards directory')

      parser.add_option('--coord', dest='coordinates',
                         type='string', default="KD1:KD2:KD3", help='Template dimensions (mass:KD1:KD2/3/int etc. Ignore the CMS_zz(*)_ prefix.)')

      parser.add_option('--GHmodel', type='int', dest='GHmodel', default=1,
                         help='GH model. 0: No GH in muF or muV, 1: Add GH/GHSM as a multiplicative factor, 2: Add GH as a multiplicative factor, -1: Add GH/GHSM as a divisive factor, -2: Add GH as a divisive factor')
      parser.add_option('--GHrefval', type='float', dest='GHrefval', default=4.07, help='GH MC reference')

      parser.add_option('-c', '--CatScheme', type='string', dest='iCatScheme', default='vbfvhcat',
                         help='Categorization scheme. See CategoryHelper.py for the index reference and naming conventions.')

      parser.add_option("--channel", dest="customChannels", type="string", action="append", help="Channels to run (default=all turned on)")
      parser.add_option("--category", dest="customCategories", type="string", action="append", help="Categories to run (default=all turned on)")

      parser.add_option('-r', '--datadir', type='string', dest='customDataDir', default="",
                         help='customDataDir: Reference data folder per measurement')

      parser.add_option('--ac', '--AnomCoupl', type='int', dest='anomCouplIndex', default="0",
                         help='anomCouplIndex: 0: SM-only, 1: Full parameterization, 2: Parameterization with phi=0 or pi')

      parser.add_option('--mPOLE', '--mH', type='float', dest='mPOLE', default="125",
                         help='mPOLE: Pole mass for the Higgs')

      parser.add_option('--mLow', type='float', dest='mLow', default="220",
                         help='mLow: Low mass boundary (def=220)')
      parser.add_option('--mHigh', type='float', dest='mHigh', default="1600",
                         help='mHigh: High mass boundary (def=1600)')

      parser.add_option('--checkpdfs', type="string", action="append", help='Check specified pdfs explicitly')
      parser.add_option('--writeoutput', action='store_true', default=False, help='Write output to file')
      parser.add_option('--dopdfproj', action='store_true', default=False, help='Do pdf projections for cross checks')

       # store options and arguments as global variables
      (self.opt, self.args) = parser.parse_args()

      if (self.opt.coordinates == ''):
         print('makeWidthDatacards: Please pass template dimensions! Exiting...')
         sys.exit()
      else:
         self.opt.coordList = self.opt.coordinates.split(":")
         if self.opt.extMassShapes is not None and self.opt.coordList[0]!="mass":
            print("makeWidthDatacards: First coordinate has to be mass when external mass shapes are given. Listed coordinates:")
            print(self.opt.coordList)
            sys.exit()
         self.opt.dimensions = len(self.opt.coordList)

      if (self.opt.appendName == ''):
         print('makeWidthDatacards: Please pass an append name for the cards directory! Exiting...')
         sys.exit()

      if (self.opt.inputDir == ''):
         print('makeWidthDatacards: Please pass an input directory! Exiting...')
         sys.exit()

      if (self.opt.templateDir == ''):
         print('makeWidthDatacards: Please pass an input directory for the templates! Exiting...')
         sys.exit()

      print("makeWidthDatacards: Categorization scheme: {0}".format(self.opt.iCatScheme))
      print("makeWidthDatacards: AC scheme: {0}".format(self.opt.anomCouplIndex))
      print("makeWidthDatacards: GH model: {0}".format(self.opt.GHmodel))


   # define make directory function
   def makeDirectory(self,subDirName):
      if (not os.path.exists(subDirName)):
         cmd = 'mkdir -p ' + subDirName
         status, output = commands.getstatusoutput(cmd)
         if status != 0:
            print('makeWidthDatacards::makeDirectory: Error in creating submission dir',subDirName,'.')
            sys.exit()
      else:
         print('makeWidthDatacards::makeDirectory: Directory',subDirName,'already exists.')


   # define function for processing of os command
   def processCmd(self,cmd):
      status, output = commands.getstatusoutput(cmd)
      if status != 0:
         print("makeWidthDatacards::processCmd: Error in processing command:\n   '" + cmd + "' \nExiting...")
         sys.exit()


   def creationLoop(self,theOutputDir):
      finalstates = [ "4mu" , "4e" , "2e2mu", "2mu2e", "2e2nu", "2mu2nu", "enumunu", "enuenu", "munumunu" ]
      CatHelper = CategoryHelper(self.opt.iCatScheme)
      for iCat in range(0,CatHelper.nCategories):
         catname=CatHelper.catNameList[iCat]
         if self.opt.customCategories is not None:
            if not catname in self.opt.customCategories:
               continue
         for ifs in finalstates:
            if self.opt.customChannels is not None:
               if not ifs in self.opt.customChannels:
                  continue
            inputCardFile = self.opt.inputDir + "/inputs_" + ifs + "_" + catname + ".txt"
            if os.path.isfile(inputCardFile):
               theInputCard = InputCardReader(inputCardFile)

               pathToDatacards = "{0}/HCG/{1}/".format(theOutputDir, theInputCard.theSqrtsPeriod)
               print("makeWidthDatacards::creationLoop: Path to datacards:",pathToDatacards)
               self.makeDirectory(pathToDatacards)

               pathToPlots = "{0}/figs/{1}/hto{2}_{3}/".format(theOutputDir, theInputCard.theSqrtsPeriod, theInputCard.decayChanName, catname)
               print("makeWidthDatacards::creationLoop: Path to plots:",pathToPlots)
               self.makeDirectory(pathToPlots)

               stdout_original=None
               stderr_original=None
               if self.opt.writeoutput:
                  try:
                     stdout_original=sys.stdout
                     sys.stdout = open(pathToPlots + "output.log", 'w')
                  except IOError: print("makeWidthDatacards::creationLoop ERROR: Could not open output.log")
                  try:
                     stderr_original=sys.stderr
                     sys.stderr = open(pathToPlots + "errors.log", 'w')
                  except IOError: print("makeWidthDatacards::creationLoop ERROR: Could not open errors.log")

               SystHelper = SystematicsHelper(theInputCard)
               theEqnsMaker = EquationsMaker(self.opt,theInputCard)
               theMaker = WidthDatacardMaker(self.opt,theInputCard,theEqnsMaker,CatHelper,SystHelper,iCat,theOutputDir)

               if self.opt.writeoutput:
                  if stdout_original is not None: sys.stdout=stdout_original
                  if stderr_original is not None: sys.stderr=stderr_original
            else:
               print("makeWidthDatacards::creationLoop: Input card {} does not exist.".format(inputCardFile))


# run the create_RM_cfg() as main()
if __name__ == "__main__":
   maker = makeWidthDatacards()


