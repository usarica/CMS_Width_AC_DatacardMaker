#!/usr/bin/python
import pwd
import commands
import optparse
import shlex
from WidthDatacardMaker import *


def loadIncludes():
   # ---------------- LOAD THE INCLUDES ----------------- ##
   ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
   ROOT.gSystem.AddIncludePath("-Iinclude/")
   ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
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
                         type='string', default="",    help='Inputs directory')
      parser.add_option('-t', '--templates', dest='templateDir',
                         type='string', default="", help='Directory of templates')
      parser.add_option('--extMassShapes', type='string', default=None, help='Directory of external mass shapes')
      parser.add_option('-a', '--append', dest='appendName',
                         type='string', default="",    help='Append name for cards directory')

      parser.add_option('--coord', dest='coordinates',
                         type='string', default="KD1:KD2:KD3",    help='Template dimensions (mass:KD1:KD2/3/int etc. Ignore CMS_zz4l_ prefix.)')

      parser.add_option('--GHmodel', type='int', dest='GHmodel', default=1,
                         help='GH model. 0: No GH in muF or muV, 1: Add GH/GHSM as a multiplicative factor, 2: Add GH as a multiplicative factor, -1: Add GH/GHSM as a divisive factor, -2: Add GH as a divisive factor')
      parser.add_option('--GHrefval', type='float', dest='GHrefval', default=4.07, help='GH MC reference')

      parser.add_option('-c', '--CatScheme', type='int', dest='iCatScheme', default=1,
                         help='Categorization scheme. 0: Inclusive, 1: VBF-tagged and untagged, 2: VBF-, VH-tagged and untagged')

      parser.add_option('-r', '--datadir', type='string', dest='dataDirAppend', default="",
                         help='dataDirAppend: Reference CMSdata folder per measurement')

      parser.add_option('--ac', '--AnomCoupl', type='int', dest='anomCouplIndex', default="0",
                         help='anomCouplIndex: 0: SM-only, 1: Full parameterization, 2: Parameterization with phi=0 or pi')

      parser.add_option('--mPOLE', '--mH', type='float', dest='mPOLE', default="125",
                         help='mPOLE: Pole mass for the Higgs')

      parser.add_option('--mLow', type='float', dest='mLow', default="220",
                         help='mLow: Low m4l boundary (def=220)')
      parser.add_option('--mHigh', type='float', dest='mHigh', default="1600",
                         help='mHigh: High m4l boundary (def=1600)')

       # store options and arguments as global variables
      (self.opt, self.args) = parser.parse_args()

      if (self.opt.coordinates == ''):
         print 'Please pass template dimensions! Exiting...'
         sys.exit()
      else:
         self.opt.coordList = self.opt.coordinates.split(":")
         if self.opt.extMassShapes is not None and self.opt.coordList[0]!="mass":
            print "First coordinate has to be mass when external mass shapes are given."
            sys.exit()
         self.opt.dimensions = len(self.opt.coordList)

      if (self.opt.appendName == ''):
         print 'Please pass an append name for the cards directory! Exiting...'
         sys.exit()

      if (self.opt.inputDir == ''):
         print 'Please pass an input directory! Exiting...'
         sys.exit()

      if (self.opt.templateDir == ''):
         print 'Please pass an input directory for the templates! Exiting...'
         sys.exit()


   # define make directory function
   def makeDirectory(self,subDirName):
      if (not os.path.exists(subDirName)):
         cmd = 'mkdir -p ' + subDirName
         status, output = commands.getstatusoutput(cmd)
         if status != 0:
            print 'Error in creating submission dir ' + subDirName + '.'
            sys.exit()
      else:
         print 'Directory ' + subDirName + ' already exists.'


   # define function for processing of os command
   def processCmd(self,cmd):
      status, output = commands.getstatusoutput(cmd)
      if status != 0:
         print 'Error in processing command:\n   [' + cmd + '] \nExiting...'
         sys.exit()


   def creationLoop(self,theOutputDir):
      CatHelper = CategoryHelper(self.opt.iCatScheme)

      for iCat in range(0,CatHelper.nCategories):
         finalstates = [ "4mu" , "4e" , "2e2mu" ]
         for ifs in finalstates:
            inputCardDir = self.opt.inputDir + "/inputs_" + ifs + "_" + CatHelper.catNameList[iCat] + ".txt"
            theInputCard = InputCardReader(inputCardDir)

            pathToDatacards = "{0}/HCG/{1}/".format(theOutputDir, theInputCard.theSqrtsPeriod)
            print "Path to datacards:",pathToDatacards
            self.makeDirectory(pathToDatacards)

            pathToPlots = "{0}/figs/{1}/hzz{2}_{3}/".format(theOutputDir, theInputCard.theSqrtsPeriod, theInputCard.decayChanName, CatHelper.catNameList[iCat])
            print "Path to plots:",pathToPlots
            self.makeDirectory(pathToPlots)

            SystHelper = SystematicsHelper(theInputCard)
            theEqnsMaker = EquationsMaker(self.opt,theInputCard)
            theMaker = WidthDatacardMaker(self.opt,theInputCard,theEqnsMaker,CatHelper,SystHelper,iCat,theOutputDir)


# run the create_RM_cfg() as main()
if __name__ == "__main__":
   maker = makeWidthDatacards()


