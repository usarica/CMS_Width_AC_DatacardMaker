#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array

class ExtendedTemplate:

   def __init__(self, origTemplate, nDimensions, ProjDim, varX, varY, varZ, condDim=None):
      self.dimensions = nDimensions
      self.ProjDim = -1
      self.condDim = -1
      #self.ProjDim = ProjDim
      #if condDim is None:
      #   self.condDim = -1
      #else if condDim<self.dimensions:
      #   self.condDim = condDim

      self.obslist = ROOT.RooArgList()
      theVars = [ varX, varY, varZ ]
      for v in range(0,len(theVars)):
         if v<self.dimensions and theVars[v] is not None:
            self.obslist.add(theVars[v])


      self.origTemplate = origTemplate # Just to store
      TemplateName = self.origTemplate.GetName() # Template name contains extra suffix if already cloned
      HistFuncName = "{}_HistFunc".format(TemplateName)
      self.getHistFunc(HistFuncName)

      # Construct the rate
      integral = float(self.theTemplate.IntegralWidth()) # Ensure that the precision is floating-point
      RateName = "{}_Rate".format(TemplateName)
      self.theRate = ROOT.RooConstVar(RateName,RateName,integral)


   def getHistFunc(self,newname):
      if self.dimensions==3:
         self.theTemplate = FastHisto3D_f(self.origTemplate)
         self.theHistFunc = FastHisto3DFunc_f(newname,"",self.obslist,self.theTemplate)
      elif self.dimensions==2:
         self.theTemplate = FastHisto2D_f(self.origTemplate)
         self.theHistFunc = FastHisto2DFunc_f(newname,"",self.obslist,self.theTemplate)
      else:
         self.theTemplate = FastHisto_f(self.origTemplate)
         self.theHistFunc = FastHistoFunc_f(newname,"",self.obslist,self.theTemplate)
