#!/usr/bin/python
import os
import re
import math
import collections
#from ROOT import *
from array import array

## ---------------------------------------------------------------
## card reader class
## ---------------------------------------------------------------

class InputCardReader:

   def __init__(self, inputPath):
      if not os.path.exists(inputPath):
         raise RuntimeError("File {0} does not exist!!!".format(inputPath))

      # input file
      self.theInput = inputPath
      # decay channel
      self.decayChan = None
      # lumi
      self.lumi = None
      # sqrts
      self.sqrts = None

      # list of
      # [ channel name, rate, lumi, (int) iBkg ]
      self.channels = []

      # list of
      # [ parameter name, value ]
      self.parameters = []

      # list of either
      # [ systematics name, systematics type=lnN, [ [channel, 1+sigma] ] ] --> e.g. systematics mySyst lnN ggH:1.02 qqH:1.04
      # or
      # [ systematics name, systematics type=lnN, [ [channel, 1+sigma, 1-sigma] ] ] --> e.g. systematics mySyst lnN ggH:1.02:0.96 qqH:1.04:0.90
      # or
      # [ systematics name, systematics type=param, [central value, 1 sigma error, (optional) parameter minimum, (optional) parameter maximum] ] --> e.g. systematics mySyst param 0:1:-3:3
      # or
      # [ systematics name, systematics type=template, [ [ channel, systematics up appendix name, systematics dn appendix name ] ] --> e.g. systematics mySyst template ggZZ_offshell:QCDUp:QCDDn VBF_offhshell:QCDUp:QCDDn
      self.systematics = []


   def parseBoolString(self,theString):
      return theString[0].upper()=='T'

   def isGoodEntry(self, var):
      if (var is None):
         return False
      elif (var == []):
         return False
      else:
         return True

   def readInputs(self):
      for line in open(self.theInput,'r'):
         f = line.split()

         if len(f) < 1: continue
         if f[0].startswith("#"): continue

         if f[0].lower().startswith("lumi"):
            self.lumi = float(f[1])

         if f[0].lower().startswith("sqrts"):
            self.sqrts = float(f[1])

         if f[0].lower().startswith("decay"):
            if f[1] == "4mu": self.decayChan = 1
            elif f[1] == "4e": self.decayChan = 2
            elif f[1] == "2e2mu": self.decayChan = 3
            elif f[1] == "2mu2e": self.decayChan = 4
            else: raise RuntimeError("Unknown decay channel {0}, choices are 4mu, 4e, 2e2mu or 2mu2e".format(f[1]))

         if f[0].lower().startswith("channel"):
            channame = f[1]
            chanrate = 1.0
            chanlumi = -1.0
            iBkg = 0
            if len(f)>2:
               chanrate = float(f[2])
               if len(f)>3:
                  chanlumi = float(f[3])
            chanlist = [ channame, chanrate, chanlumi, iBkg ]
            self.channels.append(chanlist)

         if f[0].lower().startswith("parameter"):
            parname = f[1]
            parvalue = float(f[2])
            parlist = [ parname, parvalue ]
            self.parameters.append(parlist)

         if f[0].lower().startswith("systematics"):
            parname = f[1]
            partype = f[2]
            parconfig = []
            if partype == "lnN":
               if len(f)<4:
                  raise RuntimeError("{0} uncertainty for systematic {1} is not given any process!".format(partype, parname))
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":")
                  for itc in range(1,len(tmpconfig)): # convert strings to floats
                     tmpconfig[itc] = float(tmpconfig[itc])
                  parconfig.append(tmpconfig)
            elif partype == "param":
               if len(f)!=4:
                  raise RuntimeError("{0} uncertainty for systematic {1} has to consist of 4 whitespace-separated string!".format(partype, parname))
               tmpconfig = f[4].split(":")
               for itc in range(0,len(tmpconfig)): # convert strings to floats
                  tmpconfig[itc] = float(tmpconfig[itc])
               parconfig = tmpconfig
            elif partype.lower() == "template":
               if len(f)<4:
                  raise RuntimeError("{0} uncertainty name strings for systematic {1} is not given any process!".format(partype, parname))
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":") # Leave strings as strings
                  if ((tmpconfig[0] == "") or (tmpconfig[1] == "") or (tmpconfig[2] == "")):
                     raise RuntimeError("{0} uncertainty does not specify any process or template appendix names!".format(parname, tmpconfig[0]))
                  parconfig.append(tmpconfig)

            parlist = [ parname, partype, parconfig ]
            self.parameters.append(parlist)

   def getInputs(self):

      theDict = {}

      if not self.isGoodEntry(self.sqrts):
         raise RuntimeError("{0} is not set. Check inputs!".format("sqrts"))
      if not self.isGoodEntry(self.lumi):
         raise RuntimeError("{0} is not set. Check inputs!".format("lumi"))
      if not self.isGoodEntry(self.decayChan):
         raise RuntimeError("{0} is not set. Check inputs!".format("decay"))

      if not self.isGoodEntry(self.channels):
         raise RuntimeError("{0} is empty. Check inputs!".format("channels"))
      if not self.isGoodEntry(self.parameters):
         raise RuntimeError("{0} is empty. Check inputs!".format("parameters"))
      if not self.isGoodEntry(self.systematics):
         raise RuntimeError("{0} is empty. Check inputs!".format("systematics"))

      # Make the dictionaries of each channel, parameter, systematic
      for par in self.channels:
         theDict[par[0]] = par[1:]
      for par in self.parameters:
         theDict[par[0]] = par[1:]
      for par in self.systematics:
         theDict[par[0]] = par[1:]

      return theDict
