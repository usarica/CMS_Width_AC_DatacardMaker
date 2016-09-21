#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.08.30
# by Matt Snowball
#-----------------------------------------------
import sys
import os
import pwd
import commands
import optparse
import shlex
import re
import math
from ROOT import *
import ROOT
from array import array
from InputCardReader import *
from CategoryHelper import *
from SystematicsHelper import *
from WidthDatacardMaker import *


def parseOptions():

   usage = ('usage: %prog [options] datasetList\n'
            + '%prog -h for help')
   parser = optparse.OptionParser(usage)

   parser.add_option('-b', action='store_true',
                      dest='noX', default=True, help='no X11 windows')

   parser.add_option('-i', '--input', dest='inputDir',
                      type='string', default="",    help='inputs directory')
   parser.add_option('-t', '--templates', dest='templateDir',
                      type='string', default="", help='directory for templates')
   parser.add_option('-a', '--append', dest='appendName',
                      type='string', default="",    help='append name for cards dir')

   parser.add_option('-d', '--dimension', type='int', dest='dimensions',
                      default="3", help='Template dimensions>0')
   parser.add_option('-p', '--projDim', type='int', dest='ProjDim',
                      default="-1", help='-1->2D/3D(m4l,KD,KD2), 0->1D(m4l), 1->1D(KD), 2->1D(KD2)')
   parser.add_option('--NoBkgSigInterf', type='int', dest='iBkgSigOnly', default=0,
                      help='Bkg-sig interference. 0: No interference, 1: Add interference terms')
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

   # Make a boolean flag out of the int version
   opt.isBkgSigOnly = (opt.iBkgSigOnly!=0)

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
         print 'Error in creating submission dir ' + subDirName + '. Exiting...'
         sys.exit()
   else:
      print 'Directory ' + subDirName + ' already exists. Exiting...'


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

   for iCat in range(0,self.CatHelper.nCategories):
      finalstates = [ "4mu","4e","2e2mu"]
      for ifs in finalstates:
         inputCardDir = opt.inputDir + "/inputs_" + ifs + "_" + CatHelper.catNameList[iCat] + ".txt"
         theInputCard = InputCardReader(inputCardDir)
         theInputCard.readInputs()

         SystHelper = SystematicsHelper(theInputCard)

         theMaker = WidthDatacardClass(opt,theInputCard,CatHelper,SystHelper,iCat,theOutputDir)



   if (opt.useDjet == 0):
      myReader4e = inputReader(opt.inputDir + "/inputs_4e.txt")
      myReader4e.readInputs()
      theInputs4e = myReader4e.getInputs()

      myReader4mu = inputReader(opt.inputDir + "/inputs_4mu.txt")
      myReader4mu.readInputs()
      theInputs4mu = myReader4mu.getInputs()

      myReader2e2mu = inputReader(opt.inputDir + "/inputs_2e2mu.txt")
      myReader2e2mu.readInputs()
      theInputs2e2mu = myReader2e2mu.getInputs()

   if (opt.useDjet == 1):
      myReader4e_0 = inputReader(opt.inputDir + "_tagged/inputs_4e_0.txt")
      myReader4e_0.readInputs()
      theInputs4e_0 = myReader4e_0.getInputs()

      myReader4mu_0 = inputReader(opt.inputDir + "_tagged/inputs_4mu_0.txt")
      myReader4mu_0.readInputs()
      theInputs4mu_0 = myReader4mu_0.getInputs()

      myReader2e2mu_0 = inputReader(opt.inputDir + "_tagged/inputs_2e2mu_0.txt")
      myReader2e2mu_0.readInputs()
      theInputs2e2mu_0 = myReader2e2mu_0.getInputs()

      myReader4e_1 = inputReader(opt.inputDir + "_tagged/inputs_4e_1.txt")
      myReader4e_1.readInputs()
      theInputs4e_1 = myReader4e_1.getInputs()

      myReader4mu_1 = inputReader(opt.inputDir + "_tagged/inputs_4mu_1.txt")
      myReader4mu_1.readInputs()
      theInputs4mu_1 = myReader4mu_1.getInputs()

      myReader2e2mu_1 = inputReader(
      opt.inputDir + "_tagged/inputs_2e2mu_1.txt")
      myReader2e2mu_1.readInputs()
      theInputs2e2mu_1 = myReader2e2mu_1.getInputs()


   if (opt.useDjet == 0):
      makeDirectory(directory + '/HCG')
      myClass.makeCardsWorkspaces(opt, directory, theInputs4e, 0)
      myClass.makeCardsWorkspaces(opt, directory, theInputs4mu, 0)
      myClass.makeCardsWorkspaces(opt, directory, theInputs2e2mu, 0)
   if (opt.useDjet == 1):
      makeDirectory(directory + '_tagged/HCG')
      myClass.makeCardsWorkspaces(opt, directory + '_tagged', theInputs4e_0, 1)
      myClass.makeCardsWorkspaces(opt, directory + '_tagged', theInputs4mu_0, 1)
      myClass.makeCardsWorkspaces(opt, directory + '_tagged', theInputs2e2mu_0, 1)
      myClass.makeCardsWorkspaces(opt, directory + '_tagged', theInputs4e_1, 2)
      myClass.makeCardsWorkspaces(opt, directory + '_tagged', theInputs4mu_1, 2)
      myClass.makeCardsWorkspaces(opt, directory + '_tagged', theInputs2e2mu_1, 2)



# the main procedure
def makeWidthDatacards():
   loadIncludes()

   # parse the arguments and options
   global opt, args
   parseOptions()

   dirName = 'cards_' + opt.appendName
   if (opt.iCatScheme != 0): dirName = dirName + '_tagged'
   subdir = ['HCG', 'figs']
   for d in subdir:
      makeDirectory(dirName + '/' + d)
   creationLoop(dirName)

   sys.exit()


# run the create_RM_cfg() as main()
if __name__ == "__main__": makeWidthDatacards()


