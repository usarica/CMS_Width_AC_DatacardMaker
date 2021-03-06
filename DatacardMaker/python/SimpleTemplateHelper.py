#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array
from CMS_Width_AC_DatacardMaker.DatacardMaker.CategoryHelper import CategoryHelper
from CMS_Width_AC_DatacardMaker.DatacardMaker.ExtendedTemplate import ExtendedTemplate

class SimpleTemplateHelper:
   def __init__(self, options, theMaker, theEqnsMaker, theCategorizer, theProcess, templateFileName, iCat, systName):
      self.condDim = 0
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.theSqrtsPeriod = theMaker.theSqrtsPeriod
      self.channel = theMaker.channel
      self.theChannelName = theMaker.theChannelName
      self.workspace = theMaker.workspace
      self.theInputs = theMaker.theInputCard.getInputs()

      # RooRealVars from the datacard maker class
      self.mass = theMaker.mass
      self.KD1 = theMaker.KD1
      self.KD2 = theMaker.KD2
      self.KD3 = theMaker.KD3

      self.templateDir = options.templateDir
      self.dimensions = options.dimensions # Number of template dimensions>0

      self.iCatScheme = theCategorizer.iCatScheme
      self.catNameList = theCategorizer.catNameList
      self.nCategories = theCategorizer.nCategories
      self.iCat = iCat
      if self.iCat>=self.nCategories:
         sys.exit("self.iCat={} >= self.nCategories={}!".format(self.iCat,self.nCategories))

      self.procname = theProcess[0]
      self.proctype = theProcess[3]
      self.procopts = theProcess[4]
      self.isBkgOnly = self.proctype==0
      self.isSigOnly = self.proctype==1
      self.forceOnshell = False
      self.procTplAlias=self.procname
      self.condDim = 0

      if not self.isBkgOnly and not self.isSigOnly:
         raise RuntimeError("SimpleTemplateHelper only supports bkg-only or sig-only templates, so please revise process {}".format(self.procname))

      for procopt in self.procopts:
         procoptl=procopt.lower()
         if "conditional" in procoptl:
            # Note: kd1/2/3 here refer only to self.KD1/2/3, so KD1 could actually be mass in the EquationMaker.
            self.condDim = 2**int("kd1" in procoptl) * 3**int("kd2" in procoptl) * 5**int("kd3" in procoptl)
         if "templatenamealias" in procoptl:
            self.procTplAlias = procopt.split('=')[1]
         if "forceonshell" in procoptl:
            self.forceOnshell = True

      if self.condDim==1:
         self.condDim=0
      if self.condDim>0:
         self.condVars = ROOT.RooArgSet()
         if self.condDim%2==0:
            self.condVars.add(self.KD1)
         if self.condDim%3==0:
            self.condVars.add(self.KD2)
         if self.condDim%5==0:
            self.condVars.add(self.KD3)

      # RooRealVars from the equations maker class
      if self.forceOnshell:
         self.muF = theEqnsMaker.rrvars["muF_onshell"] # Could itself be a RooFormulaVar (e.g. Off-shell: muF = R*RF*x. On-shell: muF = R*RF)
         self.muV = theEqnsMaker.rrvars["muV_onshell"] # Could itself be a RooFormulaVar (e.g. Off-shell: muV = R*RV*x. On-shell: muV = R*RV)
      else:
         self.muF = theEqnsMaker.rrvars["muF"] # Could itself be a RooFormulaVar (e.g. Off-shell: muF = R*RF*x. On-shell: muF = R*RF)
         self.muV = theEqnsMaker.rrvars["muV"] # Could itself be a RooFormulaVar (e.g. Off-shell: muV = R*RV*x. On-shell: muV = R*RV)
      self.kbkg_gg = theEqnsMaker.rrvars["kbkg_gg"]
      self.kbkg_VBF = theEqnsMaker.rrvars["kbkg_VBF"]

      self.templateFileName = templateFileName
      self.systName = systName
      self.templateSuffix = "{0}_{1}_{2}_{3}".format(self.systName,self.catNameList[self.iCat],self.theChannelName,self.theSqrtsPeriod)

      self.templateFile = None

      self.theTpl = None
      self.thePdf = None
      self.theRate = None
      self.thePdf_extras = []
      self.theRate_extras = []


# Import the pdf
   def importToWorkspace(self):
      for xrate in self.theRate_extras:
         getattr(self.workspace, 'import')(xrate, ROOT.RooFit.RecycleConflictNodes())
      if self.theRate is not None:
         getattr(self.workspace, 'import')(self.theRate, ROOT.RooFit.RecycleConflictNodes())
      for xpdf in self.thePdf_extras:
         getattr(self.workspace, 'import')(xpdf, ROOT.RooFit.RecycleConflictNodes())
      if self.thePdf is not None:
         getattr(self.workspace, 'import')(self.thePdf, ROOT.RooFit.RecycleConflictNodes())


# Open the template files
   def openFile(self):
      print "Opening file ",self.templateFileName
      self.templateFile = ROOT.TFile.Open(self.templateFileName, "read")
      if self.templateFile is None:
         raise RuntimeError("SimpleTemplateHelper file {} is None!".format(self.templateFileName))
      elif self.templateFile.IsZombie():
         raise RuntimeError("SimpleTemplateHelper could not open file {}!".format(self.templateFileName))
# Close the template files
   def close(self):
      if self.templateFile is not None:
         if self.templateFile.IsOpen():
            self.templateFile.Close()


   def getThePdf(self):
      return self.thePdf
   def getTheRate(self):
      return self.theRate


# Get shapes for each category
   def getTemplates(self,templatePrefix="T"):
      self.openFile()

      self.templatePrefix = templatePrefix
      self.templatePrefix = "{}_{}".format(self.templatePrefix,self.procname)
      self.templateInputPrefix = templatePrefix
      self.templateInputPrefix = "{}_{}".format(self.templateInputPrefix,self.procTplAlias)

#---------- TEMPLATES AND PDFS -------------
   # Construct the p.d.f.s
      # Construct the templates
      theTpl_uncond = ExtendedTemplate(
               self.templateFile.Get(self.templateInputPrefix).Clone("{}_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions,
               self.KD1, self.KD2, self.KD3
            )
      if self.condDim>0:
         self.theRate_extras.append(theTpl_uncond)
         condDimSuffix="condDim"
         if self.condDim%2==0: condDimSuffix=condDimSuffix+"0"
         if self.condDim%3==0: condDimSuffix=condDimSuffix+"1"
         if self.condDim%5==0: condDimSuffix=condDimSuffix+"2"
         self.theTpl = ExtendedTemplate(
                  self.templateFile.Get(self.templateInputPrefix+"_"+condDimSuffix).Clone("{}_{}_{}".format(self.templatePrefix,self.templateSuffix,condDimSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3,
                  self.condDim
               )
         self.theTpl.theRate=theTpl_uncond.theRate # Rate should still correspond to unconditional tpl
      else:
         self.theTpl = theTpl_uncond
      PdfName = "{}Pdf_{}".format(self.procname, self.templateSuffix)
      RateName = "{}TotalRate_{}".format(self.procname, self.templateSuffix)
      RateRawName = "{}TotalRateRaw_{}".format(self.procname, self.templateSuffix)
      self.theTpl.theRate.SetName(RateRawName)
      # Construct the rate, RFV if necessary
      if "gg" in self.procname.lower():
         if self.isBkgOnly:
            self.theRate = ROOT.RooFormulaVar(RateName, "@0*abs(@1)", ROOT.RooArgList(self.theTpl.theRate, self.kbkg_gg))
         else:
            self.theRate = ROOT.RooFormulaVar(RateName, "@0*abs(@1)", ROOT.RooArgList(self.theTpl.theRate, self.muF))
      elif "tth" in self.procname.lower() or "bbh" in self.procname.lower() \
         or "ttvv" in self.procname.lower() or "bbvv" in self.procname.lower() \
         or "ttww" in self.procname.lower() or "bbww" in self.procname.lower() \
         or "ttzz" in self.procname.lower() or "bbzz" in self.procname.lower():
         if self.isSigOnly:
            self.theRate = ROOT.RooFormulaVar(RateName, "@0*abs(@1)", ROOT.RooArgList(self.theTpl.theRate, self.muF))
         else:
            self.theRate = self.theTpl.theRate
      elif "vbf" in self.procname.lower() \
         or "zh" in self.procname.lower() or "wh" in self.procname.lower() \
         or "vbs" in self.procname.lower() \
         or "zvv" in self.procname.lower() or "wvv" in self.procname.lower() or "vvvv" in self.procname.lower() \
         or "zzz" in self.procname.lower() or "wzz" in self.procname.lower() or "vvzz" in self.procname.lower() \
         or "zww" in self.procname.lower() or "www" in self.procname.lower() or "vvww" in self.procname.lower():
         if self.isBkgOnly:
            self.theRate = ROOT.RooFormulaVar(RateName, "@0*abs(@1)", ROOT.RooArgList(self.theTpl.theRate, self.kbkg_VBF))
         else:
            self.theRate = ROOT.RooFormulaVar(RateName, "@0*abs(@1)", ROOT.RooArgList(self.theTpl.theRate, self.muV))
      else:
         self.theRate = self.theTpl.theRate
      # Rename the rate
      self.theRate.SetName(RateName)

      self.thePdf = ROOT.RooRealFlooredSumPdf(
         PdfName, PdfName,
         ROOT.RooArgList(self.theTpl.theHistFunc),ROOT.RooArgList()
      )
      print self.thePdf.GetName(),"value =",self.thePdf.getVal()
      print self.theRate.GetName(),"rate =",self.theRate.getVal()

