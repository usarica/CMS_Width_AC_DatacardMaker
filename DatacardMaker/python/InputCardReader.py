#!/usr/bin/python
import os
import re
import math
import collections
from array import array
from CMS_Width_AC_DatacardMaker.DatacardMaker.WidthHelperFunctions import GetDataPeriodString

## ---------------------------------------------------------------
## card reader class
## ---------------------------------------------------------------

class InputCardReader:

   def __init__(self, inputPath):
      if not os.path.exists(inputPath):
         raise RuntimeError("File {0} does not exist.".format(inputPath))

      # input file
      self.theInput = inputPath
      # decay channel
      self.decayChanName = ""
      self.decayChan = None
      # lumi
      self.lumi = None
      # sqrts
      self.sqrts = None
      # Data period
      self.dataperiod = None
      # Category
      self.catName = ""

      # list of
      # [ channel name, rate, lumi, (int) iBkg ]
      self.channels = []

      # list of
      # [ parameter name, value ]
      self.parameters = []

      # list of either
      # [ systematics name, systematics type=lnN, [ [channel, +1 sigma] ] ] --> e.g. systematic mySyst lnN ggH:1.02 qqH:1.04
      # or
      # [ systematics name, systematics type=lnN, [ [channel, -1 sigma, +1 sigma] ] ] --> e.g. systematic mySyst lnN ggH:0.96:1.02 qqH:0.90:1.04
      # or
      # [ systematics name, systematics type=param, [central value, 1 sigma error, (optional) parameter minimum, (optional) parameter maximum] ] --> e.g. systematic mySyst param 0:1:-3:3
      # or
      # [ systematics name, systematics type=template, [ [ channel, systematics up appendix name, systematics dn appendix name ] ] --> e.g. systematic mySyst template ggZZ_offshell:QCDUp:QCDDn VBF_offhshell:QCDUp:QCDDn
      self.systematics = []

      self.readInputs()


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

         if f[0].lower().startswith("period"):
            try:
               self.dataperiod = float(f[1])
            except:
               self.dataperiod = str(f[1])

         if f[0].lower().startswith("decay"):
            if f[1] == "4l": self.decayChan = -1
            elif f[1] == "4mu": self.decayChan = 1
            elif f[1] == "4e": self.decayChan = 2
            elif f[1] == "2e2mu": self.decayChan = 3
            elif f[1] == "2mu2e": self.decayChan = 4
            elif f[1] == "2e2nu": self.decayChan = 5 # ZZ
            elif f[1] == "2mu2nu": self.decayChan = 6 # ZZ
            elif f[1] == "enumunu": self.decayChan = 7 # WW
            elif f[1] == "enuenu": self.decayChan = 8 # WW
            elif f[1] == "munumunu": self.decayChan = 9 # WW
            elif f[1] == "2l1e": self.decayChan = 10 # ZW
            elif f[1] == "2l1mu": self.decayChan = 11 # ZW
            else: raise RuntimeError("Unknown decay channel {0}, choices are 4l, 4mu, 4e, 2e2mu, 2mu2e, 2e2nu, 2mu2nu, enumunu, enuenu, munumunu, 2l1e, 2l1mu".format(f[1]))
            self.decayChanName = f[1]

         if f[0].lower().startswith("category"):
            self.catName = f[1]

         if f[0].lower().startswith("channel"):
            channame = f[1]
            chanrate = 1.0
            chanlumi = -1.0
            iBkg = 0 # 0 bkg-only, 1 for sig and 2 for BSI
            chanopts=[]
            if len(f)>2:
               chanrate = float(f[2])
               if len(f)>3:
                  chanlumi = float(f[3])
                  if len(f)>4:
                     iBkg = int(f[4])
                     if len(f)>5:
                        strchanopts=f[5]
                        if strchanopts.lower().startswith("options:"): # Options could be "Options:Conditional=KD1;FileNameAlias=ggZZ"
                           strchanopts=strchanopts[8:]
                           chanopts=strchanopts.split(';')

            chanlist = [ channame, chanrate, chanlumi, iBkg, chanopts ]
            self.channels.append(chanlist)

         if f[0].lower().startswith("parameter"):
            parname = f[1]
            parvalue = float(f[2])
            parlist = [ parname, parvalue ]
            self.parameters.append(parlist)

         if f[0].lower().startswith("systematic"):
            parname = f[1]
            partype = f[2]
            parconfig = []
            paropts = []
            if partype == "lnN":
               if len(f)<4:
                  raise RuntimeError("{0} uncertainty for systematic {1} is not given any process!".format(partype, parname))
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":")
                  for itc in range(1,len(tmpconfig)): # convert strings to floats, tmpconfig[0] is the process name
                     tmpconfig[itc] = float(tmpconfig[itc])
                  parconfig.append(tmpconfig)
            elif partype == "param":
               if len(f)!=4:
                  raise RuntimeError("{0} uncertainty for systematic {1} has to consist of 4 whitespace-separated strings!".format(partype, parname))
               tmpconfig = f[3].split(":")
               for itc in range(0,len(tmpconfig)): # convert strings to floats
                  tmpconfig[itc] = float(tmpconfig[itc])
               parconfig = tmpconfig
            elif partype == "quadN":
               if len(f)<4:
                  raise RuntimeError("{0} uncertainty for systematic {1} is not given any process!".format(partype, parname))
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":")
                  for itc in range(1,len(tmpconfig)): # convert strings to floats, tmpconfig[0] is the process name or range specification
                     tmpconfig[itc] = float(tmpconfig[itc])
                  parconfig.append(tmpconfig)
            elif partype.lower() == "template":
               if len(f)<4:
                  raise RuntimeError("{0} uncertainty name strings for systematic {1} is not given any process!".format(partype, parname))
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":") # Leave strings as strings
                  if (tmpconfig[0] == ""):
                     raise RuntimeError("{0} uncertainty does not specify any process or option flag!".format(parname))
                  if (tmpconfig[0].lower() == "options"): # Can specify something like options:normonly=ggZZ_offshell,VVZZ_offshell etc.
                     if len(tmpconfig)!=2:
                        raise RuntimeError("{0} uncertainty options cannot be split by ':'!".format(parname))
                     else:
                        paroptslist = tmpconfig[1].split(';')
                        for paroptraw in paroptslist:
                           paroptpair = paroptraw.split('=')
                           if len(paroptpair)==1:
                              paroptpair.append("all")
                           paropts.append(paroptpair)
                  else:
                     if ((tmpconfig[1] == "") or (tmpconfig[2] == "")):
                        raise RuntimeError("{0} uncertainty does not specify any central value or 1-sigma range!".format(parname))
                     parconfig.append(tmpconfig)

            parlist = [ parname, partype, parconfig, paropts ]
            self.systematics.append(parlist)
      if not self.isGoodEntry(self.sqrts):
         raise RuntimeError("{0} is not set. Check inputs!".format("sqrts"))
      else:
         if self.sqrts>=13:
            if not self.isGoodEntry(self.dataperiod):
               raise RuntimeError("{0} is not set. Check inputs!".format("period"))
         self.theSqrtsPeriod=GetDataPeriodString(self.sqrts,self.dataperiod)
         self.theSqrts=GetDataPeriodString(self.sqrts,None)
         if self.sqrts==13:
            self.theSqrts_2015_2016=GetDataPeriodString(self.sqrts,"2015_2016") # Needed for the special unc. for 2015 and 2016 corelated
            self.theSqrts_2017_2018=GetDataPeriodString(self.sqrts,"2017_2018") # Needed for the special unc. for 2017 and 2018 corelated
         else:
            self.theSqrts_2015_2016=None
            self.theSqrts_2017_2018=None

      if not self.isGoodEntry(self.lumi):
         raise RuntimeError("{0} is not set. Check inputs!".format("lumi"))
      if not self.isGoodEntry(self.decayChan):
         raise RuntimeError("{0} is not set. Check inputs!".format("decay"))

      if not self.isGoodEntry(self.channels):
         raise RuntimeError("{0} is empty. Check inputs!".format("channels"))
      #if not self.isGoodEntry(self.parameters):
      #   raise RuntimeError("{0} is empty. Check inputs!".format("parameters"))
      #if not self.isGoodEntry(self.systematics):
      #   raise RuntimeError("{0} is empty. Check inputs!".format("systematics"))

   def getInputs(self):

      theDict = {}

      # Make the dictionaries of each channel, parameter, systematic
      for par in self.channels:
         theDict[par[0]] = par[1:]
      for par in self.parameters:
         theDict[par[0]] = par[1:]
      for par in self.systematics:
         theDict[par[0]] = par[1:]

      return theDict

   def getNSigProcs(self):
      ctr = 0
      for proc in self.channels:
         if proc[3]>0:
            ctr += 1
      return ctr

   def getNBkgProcs(self):
      ctr = 0
      for proc in self.channels:
         if proc[3]==0:
            ctr += 1
      return ctr


