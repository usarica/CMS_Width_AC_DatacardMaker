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

   def __init__(self, origTemplate, nDimensions, varX, varY, varZ, condDim=None):
      self.dimensions = nDimensions
      self.condDim = -1
      if condDim is not None:
         self.condDim=condDim

      self.origTemplate = origTemplate # Just to store

      self.obslist = ROOT.RooArgList()
      theVars = [ varX, varY, varZ ]
      for v in range(0,len(theVars)):
         if theVars[v] is not None:
            self.obslist.add(theVars[v])
            # FIXME: Revise for variable binning
            axis=None
            if v==0:
               axis=self.origTemplate.GetXaxis()
            elif v==1:
               axis=self.origTemplate.GetYaxis()
            else:
               axis=self.origTemplate.GetZaxis()
            aBins = axis.GetNbins()
            aLow = axis.GetXmin()
            aHigh = axis.GetXmax()
            theVars[v].setBins(aBins)
            theVars[v].setRange(aLow,aHigh)

      TemplateName = self.origTemplate.GetName() # Template name contains extra suffix if already cloned
      HistFuncName = "{}_HistFunc".format(TemplateName)
      self.getHistFunc(HistFuncName)

      # Construct the rate
      integral = float(self.theTemplate.IntegralWidth()) # Ensure that the precision is floating-point
      RateName = "{}_Rate".format(TemplateName)
      self.theRate = ROOT.RooConstVar(RateName,RateName,integral)


   def getHistFunc(self,newname):
      isNormX = (self.condDim>0 and self.condDim%2==0)
      isNormY = (self.condDim>0 and self.condDim%3==0)
      isNormZ = (self.condDim>0 and self.condDim%5==0)
      if self.dimensions==3:
         self.theTemplate = FastHisto3D_f(self.origTemplate)
         self.theHistFunc = FastHisto3DFunc_f(newname,"",self.obslist,self.theTemplate, isNormX, isNormY, isNormZ)
      elif self.dimensions==2:
         self.theTemplate = FastHisto2D_f(self.origTemplate)
         self.theHistFunc = FastHisto2DFunc_f(newname,"",self.obslist,self.theTemplate, isNormX, isNormY)
      else:
         self.theTemplate = FastHisto_f(self.origTemplate)
         self.theHistFunc = FastHistoFunc_f(newname,"",self.obslist,self.theTemplate, isNormX)
