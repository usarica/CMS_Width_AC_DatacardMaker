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
from width_datacardClass2DNew import *
from inputReader import *

# define function for parsing options


def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--input', dest='inputDir',
                      type='string', default="",    help='inputs directory')
    parser.add_option('-t', '--templates', dest='templateDir',
                      type='string', default="", help='directory for templates')
    parser.add_option('-a', '--append', dest='appendName',
                      type='string', default="",    help='append name for cards dir')
    parser.add_option('-b', action='store_true',
                      dest='noX', default=True, help='no X11 windows')
    parser.add_option('-d', '--dimension', type='int', dest='dimensions',
                      default="2", help='0->1D(KD), 1->1D(m4l), 2->2D(m4l,KD)')
    parser.add_option('-j', '--Djet', type='int', dest='useDjet', default=0,
                      help='useDjet cut of 0.5 for VBF categorization (default:0)')
    parser.add_option('-l', '--uselegacy', type='int', dest='useLegacy', default="0",
                      help='useLegacy: Option to read templates from the legacy (2014 paper) template files 1=Legacy')
    parser.add_option('-r', '--datadir', type='string', dest='dataDirAppend', default="",
                      help='dataDirAppend: Reference CMSdata folder per measurement')
    parser.add_option('-c', '--AnomCoupl', type='int', dest='anomalousCouplingIndex', default="0",
                      help='anomalousCouplingIndex: 0: SM-only, 1: fLQ')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.appendName == ''):
        print 'Please pass an append name for the cards directory! Exiting...'
        sys.exit()

    if (opt.inputDir == ''):
        print 'Please pass an input directory! Exiting...'
        sys.exit()

    if (opt.templateDir == ''):
        print 'Please pass an input directory for the templates! Exiting...'
        sys.exit()

    if (opt.useDjet != 0 and opt.useDjet != 1):
        print 'The input ' + opt.useDjet + ' is not supported for useDjet. Use either 0 (no cut) or 1 (use cut). Exiting...'
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
#        sys.exit()


# define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status != 0:
        print 'Error in processing command:\n   [' + cmd + '] \nExiting...'
        sys.exit()


def creationLoop(directory):
    global opt, args

    startMass = [220]

    myClass = width_datacardClass()
    myClass.loadIncludes()
    myClass.setDimensions(opt.dimensions)

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

        myReader2e2mu_0 = inputReader(
            opt.inputDir + "_tagged/inputs_2e2mu_0.txt")
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

    a = 0
    while (a < len(startMass)):

        mh = startMass[a]
        mhs = str(mh).replace('.0', '')

        print mh

        if (opt.useDjet == 0):
            makeDirectory(directory + '/HCG/' + mhs)
            makeDirectory(directory + '/HCG_XSxBR/' + mhs)
            myClass.makeCardsWorkspaces(
                mh, opt, directory, theInputs4e)
            myClass.makeCardsWorkspaces(
                mh, opt, directory, theInputs4mu)
            myClass.makeCardsWorkspaces(
                mh, opt, directory, theInputs2e2mu)

        if (opt.useDjet == 1):
            makeDirectory(directory + '_tagged/HCG/' + mhs)
            makeDirectory(directory + '_tagged/HCG_XSxBR/' + mhs)

            myClass.makeCardsWorkspaces(
                mh, opt, directory + '_tagged', theInputs4e_0, 1)
            myClass.makeCardsWorkspaces(
                mh, opt, directory + '_tagged', theInputs4mu_0, 1)
            myClass.makeCardsWorkspaces(
                mh, opt, directory + '_tagged', theInputs2e2mu_0, 1)

            myClass.makeCardsWorkspaces(
                mh, opt, directory + '_tagged', theInputs4e_1, 2)
            myClass.makeCardsWorkspaces(
                mh, opt, directory + '_tagged', theInputs4mu_1, 2)
            myClass.makeCardsWorkspaces(
                mh, opt, directory + '_tagged', theInputs2e2mu_1, 2)

        a += 1


# the main procedure
def make_width_DCsandWSs():

    # parse the arguments and options
    global opt, args
    parseOptions()

    if (opt.appendName != ''):
        dirName = 'cards_' + opt.appendName

    subdir = ['HCG', 'HCG_XSxBR', 'figs']

    for d in subdir:
        if (opt.useDjet == 0):
            makeDirectory(dirName + '/' + d)
        if (opt.useDjet == 1):
            makeDirectory(dirName + '_tagged/' + d)

    creationLoop(dirName)

    sys.exit()


# run the create_RM_cfg() as main()
if __name__ == "__main__":
    make_width_DCsandWSs()
