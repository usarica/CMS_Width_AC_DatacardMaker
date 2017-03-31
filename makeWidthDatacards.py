#!/usr/bin/python
import pwd
import commands
import optparse
import shlex
from WidthDatacardMaker import *


def parseOptions():

   usage = ('usage: %prog [options] datasetList\n'
            + '%prog -h for help')
   parser = optparse.OptionParser(usage)

   parser.add_option('-b', action='store_true',
                      dest='noX', default=True, help='no X11 windows')

   parser.add_option('-i', '--input', dest='inputDir',
                      type='string', default="",    help='Inputs directory')
   parser.add_option('-t', '--templates', dest='templateDir',
                      type='string', default="", help='Directory of templates')
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
   global opt, args
   (opt, args) = parser.parse_args()

   if (opt.coordinates == ''):
      print 'Please pass template dimensions! Exiting...'
      sys.exit()
   else:
      opt.coordList = opt.coordinates.split(":")
      opt.dimensions = len(opt.coordList)

   if (opt.appendName == ''):
      print 'Please pass an append name for the cards directory! Exiting...'
      sys.exit()

   if (opt.inputDir == ''):
      print 'Please pass an input directory! Exiting...'
      sys.exit()

   if (opt.templateDir == ''):
      print 'Please pass an input directory for the templates! Exiting...'
      sys.exit()


# define make directory function
def makeDirectory(subDirName):
   if (not os.path.exists(subDirName)):
      cmd = 'mkdir -p ' + subDirName
      status, output = commands.getstatusoutput(cmd)
      if status != 0:
         print 'Error in creating submission dir ' + subDirName + '.'
         sys.exit()
   else:
      print 'Directory ' + subDirName + ' already exists.'


# define function for processing of os command
def processCmd(cmd):
   status, output = commands.getstatusoutput(cmd)
   if status != 0:
      print 'Error in processing command:\n   [' + cmd + '] \nExiting...'
      sys.exit()


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


def creationLoop(theOutputDir):
   global opt, args

   CatHelper = CategoryHelper(opt.iCatScheme)

   for iCat in range(0,CatHelper.nCategories):
      finalstates = [ "4mu" , "4e" , "2e2mu" ]
      for ifs in finalstates:
         inputCardDir = opt.inputDir + "/inputs_" + ifs + "_" + CatHelper.catNameList[iCat] + ".txt"
         theInputCard = InputCardReader(inputCardDir)

         pathToDatacards = "{0}/HCG/{1:.0f}TeV/".format(theOutputDir, theInputCard.sqrts)
         print "Path to datacards:",pathToDatacards
         makeDirectory(pathToDatacards)

         SystHelper = SystematicsHelper(theInputCard)
         theEqnsMaker = EquationsMaker(opt,theInputCard)
         theMaker = WidthDatacardMaker(opt,theInputCard,theEqnsMaker,CatHelper,SystHelper,iCat,theOutputDir)


# the main procedure
def makeWidthDatacards():
   loadIncludes()

   # parse the arguments and options
   global opt, args
   parseOptions()

   dirName = 'cards_' + opt.appendName
   subdir = ['HCG', 'figs']
   for d in subdir:
      makeDirectory(dirName + '/' + d)
   creationLoop(dirName)

   sys.exit()


# run the create_RM_cfg() as main()
if __name__ == "__main__": makeWidthDatacards()


