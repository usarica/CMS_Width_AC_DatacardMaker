#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array

class ExtendedTemplate:

   def __init__(self, origTemplate, nDimensions, varX, varY, varZ, condDim=None):
      self.dimensions = nDimensions
      self.condDim = -1
      if condDim is not None:
         self.condDim=condDim

      self.origTemplate = origTemplate # Just to store
      TemplateName = self.origTemplate.GetName() # Template name contains extra suffix if already cloned

      self.obslist = ROOT.RooArgList()
      theVars = [ varX, varY, varZ ]
      for v in range(0,len(theVars)):
         if theVars[v] is not None:
            self.obslist.add(theVars[v])
            axis=None
            if v==0:
               axis=self.origTemplate.GetXaxis()
            elif v==1:
               axis=self.origTemplate.GetYaxis()
            else:
               axis=self.origTemplate.GetZaxis()
            axisBinning = self.getBinningFromHistogram(axis)
            print "ExtendedTemplate: {} axis {} has {} bins.".format(TemplateName,v,axisBinning.numBins())
            theVars[v].setBinning(axisBinning)
            theVars[v].setVal((axisBinning.lowBound()+axisBinning.highBound())/2.)

      HistFuncName = "{}_HistFunc".format(TemplateName)
      self.getHistFunc(HistFuncName)

      # Construct the rate
      integral = float(self.theTemplate.IntegralWidth()) # Ensure that the precision is floating-point
      RateName = "{}_Rate".format(TemplateName)
      self.theRate = ROOT.RooConstVar(RateName,RateName,integral)
      print "Template",TemplateName,"has integral =",self.theRate.getVal()


   def getBinningFromHistogram(self,axis):
      aBins = axis.GetNbins()
      aLow = axis.GetXmin()
      aHigh = axis.GetXmax()
      aBinning = ROOT.RooBinning(aLow, aHigh)
      for ix in range(1,aBins):
         aBinning.addBoundary(axis.GetBinLowEdge(ix+1))
      return aBinning


   def getHistFunc(self,newname):
      isNormX = (self.condDim>0 and self.condDim%2==0)
      isNormY = (self.condDim>0 and self.condDim%3==0)
      isNormZ = (self.condDim>0 and self.condDim%5==0)
      if self.dimensions==3:
         self.theTemplate = ROOT.FastHisto3D_f(self.origTemplate, isNormX, isNormY, isNormZ)
         self.theHistFunc = ROOT.FastHisto3DFunc_f(newname,"",self.obslist,self.theTemplate)
      elif self.dimensions==2:
         self.theTemplate = ROOT.FastHisto2D_f(self.origTemplate, isNormX, isNormY)
         self.theHistFunc = ROOT.FastHisto2DFunc_f(newname,"",self.obslist,self.theTemplate)
      else:
         self.theTemplate = ROOT.FastHisto_f(self.origTemplate, isNormX)
         self.theHistFunc = ROOT.FastHistoFunc_f(newname,"",self.obslist,self.theTemplate)
