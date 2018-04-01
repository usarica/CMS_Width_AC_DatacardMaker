#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array
from CategoryHelper import CategoryHelper
from ExtendedTemplate import ExtendedTemplate

class BkgTemplateHelper:
   def __init__(self, options, theMaker, theCategorizer, theProcess, templateFileName, iCat, systName):
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
      self.procopts = theProcess[4]
      self.condDim = 0
      for procopt in self.procopts:
         if "conditional" in procopt:
            self.condDim = 2**int("kd1" in procopt) * 3**int("kd2" in procopt) * 5**int("kd3" in procopt)
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

      self.templateFileName = templateFileName
      self.systName = systName
      self.templateSuffix = "{0}_{1}_{2}_{3}".format(self.systName,self.catNameList[self.iCat],self.theChannelName,self.theSqrtsPeriod)

      self.templateFile = None

      self.bkgTpl = None
      self.bkgPdf = None
      self.bkgPdf_extras = []


# Import the pdf
   def importToWorkspace(self):
      if self.bkgTpl is not None:
         getattr(self.workspace, 'import')(self.bkgTpl, ROOT.RooFit.RecycleConflictNodes())
      for xpdf in self.bkgPdf_extras:
         getattr(self.workspace, 'import')(xpdf, ROOT.RooFit.RecycleConflictNodes())
      if self.bkgPdf is not None:
         getattr(self.workspace, 'import')(self.bkgPdf, ROOT.RooFit.RecycleConflictNodes())


# Open the template files
   def openFile(self):
      print "Opening file ",self.templateFileName
      self.templateFile = ROOT.TFile.Open(self.templateFileName, "read")
      if self.templateFile is None:
         raise RuntimeError("BkgTemplateHelper file {} is None!".format(self.templateFileName))
      elif self.templateFile.IsZombie():
         raise RuntimeError("BkgTemplateHelper could not open file {}!".format(self.templateFileName))
# Close the template files
   def close(self):
      if self.templateFile is not None:
         if self.templateFile.IsOpen():
            self.templateFile.Close()


   def getThePdf(self):
      return self.bkgPdf
   def getTheRate(self):
      return self.bkgTpl.theRate


# Get shapes for each category
   def getTemplates(self,templatePrefix="T"):
      self.openFile()

      self.templatePrefix = templatePrefix
      self.templatePrefix = "{}_{}".format(self.templatePrefix,self.procname)

#---------- TEMPLATES AND PDFS -------------
   # Construct the p.d.f.s
      if not("gg" in self.procname.lower() or "vbf" in self.procname.lower() or "vbs" in self.procname.lower() or "qq" in self.procname.lower() or "zx" in self.procname.lower() or "zjets" in self.procname.lower()):
         print "BkgTemplateHelper::getTemplates({}): Process might not be a background!".format(self.procname)

      # Construct the templates
      self.bkgTpl = ExtendedTemplate(
               self.templateFile.Get(self.templatePrefix).Clone("{}_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions,
               self.KD1, self.KD2, self.KD3,
               self.condDim
            )
      PdfName = "{}_{}".format(self.procname, self.templateSuffix)

      if self.condDim>0:
         HistPdfName = "{}_others_{}".format(self.procname, self.templateSuffix)
         self.bkgHistPdf = ROOT.RooRealFlooredSumPdf(
            HistPdfName, HistPdfName,
            ROOT.RooArgList(self.bkgTpl.theHistFunc),ROOT.RooArgList()
         )
         self.bkgPdf_extras.append(bkgHistPdf)
         MassPdfName = "{}_mass_{}".format(self.procname, self.templateSuffix)
         MassPdfInName = "{}_mass".format(self.templatePrefix)
         bkgMassPdf = self.templateFile.Get(MassPdfInName).Clone(MassPdfName)
         self.bkgPdf_extras.append(bkgMassPdf)
         self.bkgPdf = ROOT.RooProdPdf(
            PdfName, PdfName,
            ROOT.RooArgSet( bkgMassPdf ),
            ROOT.RooFit.Conditional(
               ROOT.RooArgSet( bkgHistPdf ),
               self.condVars
            )
         )
      else:
         self.bkgPdf = ROOT.RooRealFlooredSumPdf(
            PdfName, PdfName,
            ROOT.RooArgList(self.bkgTpl.theHistFunc),ROOT.RooArgList()
         )
      print self.bkgPdf.GetName(),"value =",self.bkgPdf.getVal()
      print self.bkgTpl.theRate.GetName(),"rate =",self.bkgTpl.theRate.getVal()



