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
      self.ProjDim = ProjDim
      if condDim is None:
         self.condDim = -1
      else if condDim<self.dimensions:
         self.condDim = condDim

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
         if self.condDim==0:
            for ix in range(1,self.theTemplate.GetNbinsX()+1):
               dbin = self.theTemplate.GetXAxis().GetBinWidth(ix)
               if self.dimensions==3:
                  integral = self.theTemplate.Integral(ix,ix,1,self.theTemplate.GetNbinsY(),1,self.theTemplate.GetNbinsZ(),"width")/dbin # /dbin is to divide d_condDim
                  if integral!=0.0:
                     for iy in range(1,self.theTemplate.GetNbinsY()+1):
                        for iz in range(1,self.theTemplate.GetNbinsZ()+1):
                           self.theTemplate.SetBinContent(ix,iy,iz,self.theTemplate.GetBinContent(ix,iy,iz)/integral)
                           self.theTemplate.SetBinError(ix,iy,iz,self.theTemplate.GetBinError(ix,iy,iz)/integral)
               elif self.dimensions==2:
                  integral = self.theTemplate.Integral(ix,ix,1,self.theTemplate.GetNbinsY(),"width")
                  if integral!=0.0:
                     for iy in range(1,self.theTemplate.GetNbinsY()+1):
                        self.theTemplate.SetBinContent(ix,iy,self.theTemplate.GetBinContent(ix,iy)/integral)
                        self.theTemplate.SetBinError(ix,iy,self.theTemplate.GetBinError(ix,iy)/integral)
         elif self.condDim==1:
            for iy in range(1,self.theTemplate.GetNbinsY()+1):
               dbin = self.theTemplate.GetYAxis().GetBinWidth(iy)
               if self.dimensions==3:
                  integral = self.theTemplate.Integral(1,self.theTemplate.GetNbinsX(),iy,iy,1,self.theTemplate.GetNbinsZ(),"width")/dbin # /dbin is to divide d_condDim
                  if integral!=0.0:
                     for ix in range(1,self.theTemplate.GetNbinsX()+1):
                        for iz in range(1,self.theTemplate.GetNbinsZ()+1):
                           self.theTemplate.SetBinContent(ix,iy,iz,self.theTemplate.GetBinContent(ix,iy,iz)/integral)
                           self.theTemplate.SetBinError(ix,iy,iz,self.theTemplate.GetBinError(ix,iy,iz)/integral)
               elif self.dimensions==2:
                  integral = self.theTemplate.Integral(1,self.theTemplate.GetNbinsY(),iy,iy,"width")
                  if integral!=0.0:
                     for ix in range(1,self.theTemplate.GetNbinsX()+1):
                        self.theTemplate.SetBinContent(ix,iy,self.theTemplate.GetBinContent(ix,iy)/integral)
                        self.theTemplate.SetBinError(ix,iy,self.theTemplate.GetBinError(ix,iy)/integral)
         elif self.condDim==2:
            for iz in range(1,self.theTemplate.GetNbinsZ()+1):
               dbin = self.theTemplate.GetZAxis().GetBinWidth(iz)
               integral = self.theTemplate.Integral(1,self.theTemplate.GetNbinsX(),1,self.theTemplate.GetNbinsY(),iz,iz,"width")/dbin # /dbin is to divide d_condDim
               if integral!=0.0:
                  for ix in range(1,self.theTemplate.GetNbinsX()+1):
                     for iy in range(1,self.theTemplate.GetNbinsY()+1):
                        self.theTemplate.SetBinContent(ix,iy,iz,self.theTemplate.GetBinContent(ix,iy,iz)/integral)
                        self.theTemplate.SetBinError(ix,iy,iz,self.theTemplate.GetBinError(ix,iy,iz)/integral)

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

