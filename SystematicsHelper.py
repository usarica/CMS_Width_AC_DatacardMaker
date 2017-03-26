#! /usr/bin/env python
import os
import re
import math
from ROOT import *
import ROOT
from array import array

# Some helper functions that don't need to belong to the systematics class
def FloatToString(inputValue):
   return ('%.10f' % inputValue).rstrip('0').rstrip('.')


## ------------------------------------
##  systematics class
## ------------------------------------

class SystematicsHelper:
   def __init__(self,theInputs):
      self.sqrts = theInputs.sqrts
      self.channels = theInputs.channels
      self.systematics = theInputs.systematics
      self.systVars = []

      self.buildParamSystVars()

   def buildParamSystVars(self):
      for syst in self.systematics:
         systname = syst[0]
         systtpye = syst[1]
         systconfig = syst[2]
         theSystVar = None
         if (systtype == "param"):
            if len(systconfig)==2:
               theSystVar = ROOT.RooRealVar( systname, systname, systconfig[0], -7, 7)
            elif len(systconfig)==4:
               theSystVar = ROOT.RooRealVar( systname, systname, systconfig[0], systconfig[2], systconfig[3])
            else:
               raise RuntimeError("Systematics variable {} configuration does not have size 2 or 4!".format(systname))
         elif (systtype == "template"):
            theSystVar = ROOT.RooRealVar( systname, systname, 0, -3, 3)
         if theSystVar is not None:
            self.systVars.append(theSystVar)

   def getVariableDict(self):
      theDict = ()
      for par in self.systVars:
         theDict[par.GetName()] = par
      return theDict

   def writeSystematics(self, theFile):
      for syst in self.systematics:
         systname = syst[0]
         systtype = syst[1]
         systtype_ALT = systtype
         if systtype == "template":
            systtype_ALT = "param"
         systconfig = syst[2]
         systline = ""
         if (systtype == "param"):
            if len(systconfig)==4:
               systline = "{0} {1} [{2},{3}]".format(FloatToString(systconfig[0]), FloatToString(systconfig[1]), FloatToString(systconfig[2]), FloatToString(systconfig[3]))
            elif len(systconfig)==2:
               systline = "{0} {1} [-7,7]".format(FloatToString(systconfig[0]), FloatToString(systconfig[1]))
            else:
               raise RuntimeError("Systematics variable {} configuration does not have size 2 or 4!".format(systname))
         elif (systtype == "template"):
            systline = "0 1 [-3,3]"
         elif (systtype == "lnN"):
            for ch in self.channels:
               channameToFind = ch[0]
               chanConfig = None
               addLine = "- "
               for ic in systconfig:
                  if ( (ic[0] == channameToFind) or (ic[0] == "") ):
                     chanConfig = ic
                     break
               if chanConfig is not None:
                  if len(chanConfig)==2:
                     addLine = "{0} ".format(FloatToString(chanConfig[1]))
                  elif len(chanConfig)==3:
                     addLine = "{0}/{1} ".format(FloatToString(chanConfig[1],chanConfig[2]))
                  else:
                     raise RuntimeError("Systematics variable {} configuration for channel {} does not have size 2 or 3!".format(systname,channameToFind))
               systline = "{0}{1}".format(systline,addLine)
         else:
            raise RuntimeError("SystematicsHelper::writeSystematics does not support systematics variable type {} in variable {}.".format(systtype,systname))
         systline = "{0} {1} {2}\n".format(systname,systtype_ALT,systline)
         theFile.write(systline)


