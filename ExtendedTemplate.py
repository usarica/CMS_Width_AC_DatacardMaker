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

   def __init__(self, origTemplate, nDimensions, ProjDim, varX, varY, varZ):
      self.dimensions = nDimensions
      self.ProjDim = ProjDim

      self.arglist = ROOT.RooArgList()
      self.argset = ROOT.RooArgSet()

      self.origTemplate = origTemplate # Just to store
      TemplateName = self.origTemplate.GetName()

      # Construct the rate and integral first, dealing with templates is a little bit more involved
      self.integral = float(self.origTemplate.Integral("width")) # Ensure that the precision is floating-point
      RateName = "{0}_Rate".format(TemplateName)
      self.theRate = ROOT.RooRealVar(RateName,RateName,self.integral)
      self.theRate.setConstant(True)

      # Construct the actual template,
      # FIXME: X and Y projections do not take bin widths into account, should sum manually!
      if (self.ProjDim==0 and self.dimensions>1):
         self.theTemplate = self.origTemplate.ProjectionX("_ProjX")
         self.arglist.add(varX)
         self.argset.add(varX)
      elif (self.ProjDim==1 and self.dimensions>1):
         self.theTemplate = self.origTemplate.ProjectionY("_ProjY")
         self.arglist.add(varY)
         self.argset.add(varY)
      elif (self.ProjDim==2 and self.dimensions>2):
         self.theTemplate = self.origTemplate.ProjectionZ("_ProjZ")
         self.arglist.add(varZ)
         self.argset.add(varZ)
      else:
         self.theTemplate = self.origTemplate
         if self.dimensions>0:
            self.arglist.add(self.varX)
            self.argset.add(self.varX)
         if self.dimensions>1:
            self.arglist.add(self.varY)
            self.argset.add(self.varY)
         if self.dimensions>2:
            self.arglist.add(self.varZ)
            self.argset.add(self.varZ)

      DataHistName = "{0}_DataHist".format(TemplateName)
      self.theDataHist = ROOT.RooDataHist(DataHistName, DataHistName, self.arglist, self.theTemplate)
      HistFuncName = "{0}_HistFunc".format(TemplateName)
      self.theHistFunc = ROOT.RooHistFunc(HistFuncName, HistFuncName, self.argset, self.theDataHist)

