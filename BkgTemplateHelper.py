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
   def __init__(self, options, theMaker, theCategorizer, procname, templateFileName, iCat, systName):
      self.condDim = 0
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
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

      self.procname = procname

      self.templateFileName = templateFileName
      self.systName = systName
      self.templateSuffix = "{0}_{1}_{2}_{3:.0f}TeV".format(self.systName,self.catNameList[self.iCat],self.theChannelName,self.sqrts)

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
      # qq bkg
      if("gg" in self.procname.lower() or "vbf" in self.procname.lower() or "vbs" in self.procname.lower()):
         # Construct the templates
         self.bkgTpl = ExtendedTemplate(
                  self.templateFile.Get(self.templatePrefix).Clone("{}_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         PdfName = "{}_{}".format(self.procname, self.templateSuffix)
         self.bkgPdf = ROOT.RooRealFlooredSumPdf(
            PdfName, PdfName,
            ROOT.RooArgList(self.bkgTpl.theHistFunc),ROOT.RooArgList()
         )

      # qq bkg
      elif("qq" in self.procname.lower()):
         # Construct the templates
         self.bkgTpl = ExtendedTemplate(
                  self.templateFile.Get(self.templatePrefix).Clone("{}_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         PdfName = "{}_{}".format(self.procname, self.templateSuffix)
         self.bkgPdf = ROOT.RooRealFlooredSumPdf(
            PdfName, PdfName,
            ROOT.RooArgList(self.bkgTpl.theHistFunc),ROOT.RooArgList()
         )

      # Z+X bkg
      elif("zx" in self.procname.lower() or "zjets" in self.procname.lower()):
         self.condDim = (self.KD1==self.mass)*2 + (self.KD2==self.mass)*3 + (self.KD3==self.mass)*5
         print "Zjets template condDim =",self.condDim
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

            val_mean_3P1F = float(self.theInputs['zjetsShape_mean_3P1F'])
            val_sigma_3P1F = float(self.theInputs['zjetsShape_sigma_3P1F'])
            val_norm_3P1F = float(self.theInputs['zjetsShape_norm_3P1F'])
            var_mean_3P1F = None
            var_sigma_3P1F = None
            var_norm_3P1F = None
            pdf_3P1F = None
            if val_norm_3P1F>0.:
               tmpname = "{}_mass_3p1f_mean_{}".format(self.procname, self.templateSuffix)
               var_mean_3P1F = ROOT.RooRealVar(tmpname, "mean landau Zjet 3p1f", val_mean_3P1F)
               tmpname = "{}_mass_3p1f_sigma_{}".format(self.procname, self.templateSuffix)
               var_sigma_3P1F = ROOT.RooRealVar(tmpname, "sigma landau Zjet 3p1f", val_sigma_3P1F)
               tmpname = "{}_mass_3p1f_norm_{}".format(self.procname, self.templateSuffix)
               var_norm_3P1F = ROOT.RooRealVar(tmpname, "norm landau Zjet 3p1f", val_norm_3P1F)
               pdf_3P1F = ROOT.RooLandau(
                  "{}_mass_3p1f_{}".format(self.procname, self.templateSuffix),
                  "{}_mass_3p1f_{}".format(self.procname, self.templateSuffix),
                  self.mass,
                  var_mean_3P1F, var_sigma_3P1F
               )

            val_mean_2P2F = float(self.theInputs['zjetsShape_mean_2P2F'])
            val_sigma_2P2F = float(self.theInputs['zjetsShape_sigma_2P2F'])
            val_norm_2P2F = float(self.theInputs['zjetsShape_norm_2P2F'])
            val_pol0_2P2F = float(self.theInputs['zjetsShape_pol0_2P2F'])
            val_pol1_2P2F = float(self.theInputs['zjetsShape_pol1_2P2F'])
            var_mean_2P2F = None
            var_sigma_2P2F = None
            var_norm_2P2F = None
            var_pol0_2P2F = None
            var_pol1_2P2F = None
            pdf_2P2F = None
            if val_norm_2P2F>0.:
               tmpname = "{}_mass_2p2f_mean_{}".format(self.procname, self.templateSuffix)
               var_mean_2P2F = ROOT.RooRealVar(tmpname, "mean landau Zjet 2p2f", val_mean_2P2F)
               tmpname = "{}_mass_2p2f_sigma_{}".format(self.procname, self.templateSuffix)
               var_sigma_2P2F = ROOT.RooRealVar(tmpname, "sigma landau Zjet 2p2f", val_sigma_2P2F)
               tmpname = "{}_mass_2p2f_norm_{}".format(self.procname, self.templateSuffix)
               var_norm_2P2F = ROOT.RooRealVar(tmpname, "norm landau Zjet 2p2f", val_norm_2P2F)
               if val_pol1_2P2F==0.:
                  pdf_2P2F = ROOT.RooLandau(
                     "{}_mass_2p2f_{}".format(self.procname, self.templateSuffix),
                     "{}_mass_2p2f_{}".format(self.procname, self.templateSuffix),
                     self.mass,
                     var_mean_2P2F, var_sigma_2P2F
                  )
               else:
                  tmpname = "{}_mass_2p2f_pol0_{}".format(self.procname, self.templateSuffix)
                  var_pol0_2P2F = ROOT.RooRealVar(tmpname, "pol0 landau Zjet 2p2f", val_pol0_2P2F)
                  tmpname = "{}_mass_2p2f_pol1_{}".format(self.procname, self.templateSuffix)
                  var_pol1_2P2F = ROOT.RooRealVar(tmpname, "pol1 landau Zjet 2p2f", val_pol1_2P2F)
                  pdf_2P2F = ROOT.RooGenericPdf(
                     "{}_mass_2p2f_{}".format(self.procname, self.templateSuffix),
                     "{}_mass_2p2f_{}".format(self.procname, self.templateSuffix),
                     "(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",
                     ROOT.RooArgList(
                        self.mass,
                        var_mean_2P2F, var_sigma_2P2F,
                        var_pol0_2P2F, var_pol1_2P2F
                     )
                  )

            # For 2e2mu if needed
            val_mean_2P2F_2 = float(self.theInputs['zjetsShape_mean_2P2F_2'])
            val_sigma_2P2F_2 = float(self.theInputs['zjetsShape_sigma_2P2F_2'])
            val_norm_2P2F_2 = float(self.theInputs['zjetsShape_norm_2P2F_2'])
            var_mean_2P2F_2 = None
            var_sigma_2P2F_2 = None
            var_norm_2P2F_2 = None
            pdf_2P2F_2 = None
            if val_norm_2P2F_2>0.:
               tmpname = "{}_mass_2p2f_2_mean_{}".format(self.procname, self.templateSuffix)
               var_mean_2P2F_2 = ROOT.RooRealVar(tmpname, "mean landau Zjet 2p2f_2", val_mean_2P2F_2)
               tmpname = "{}_mass_2p2f_2_sigma_{}".format(self.procname, self.templateSuffix)
               var_sigma_2P2F_2 = ROOT.RooRealVar(tmpname, "sigma landau Zjet 2p2f_2", val_sigma_2P2F_2)
               tmpname = "{}_mass_2p2f_2_norm_{}".format(self.procname, self.templateSuffix)
               var_norm_2P2F_2 = ROOT.RooRealVar(tmpname, "norm landau Zjet 2p2f_2", val_norm_2P2F_2)
               pdf_2P2F_2 = ROOT.RooLandau(
                  "{}_mass_2p2f_2_{}".format(self.procname, self.templateSuffix),
                  "{}_mass_2p2f_2_{}".format(self.procname, self.templateSuffix),
                  self.mass,
                  var_mean_2P2F_2, var_mean_2P2F_2
               )

            MassPdfName = "{}_mass_{}".format(self.procname, self.templateSuffix)
            pdflist = ROOT.RooArgList()
            coeflist = ROOT.RooArgList()
            if ((pdf_3P1F is not None) and (var_norm_3P1F is not None)):
               pdflist.add(pdf_3P1F)
               coeflist.add(var_norm_3P1F)
               self.bkgPdf_extras.append(pdf_3P1F)
            if ((pdf_2P2F is not None) and (var_norm_2P2F is not None)):
               pdflist.add(pdf_2P2F)
               coeflist.add(var_norm_2P2F)
               self.bkgPdf_extras.append(pdf_2P2F)
            if ((pdf_2P2F_2 is not None) and (var_norm_2P2F_2 is not None)):
               pdflist.add(pdf_2P2F_2)
               coeflist.add(var_norm_2P2F_2)
               self.bkgPdf_extras.append(pdf_2P2F_2)
            bkgMassPdf = ROOT.RooAddPdf(MassPdfName,MassPdfName,pdflist,coeflist)
            self.bkgPdf_extras.append(bkgMassPdf)

            self.bkgPdf = ROOT.RooProdPdf(
               PdfName, PdfName,
               ROOT.RooArgSet( bkgMassPdf ),
               ROOT.RooFit.Conditional(
                   ROOT.RooArgSet( bkgHistPdf ),
                   ROOT.RooArgSet( self.mass )
               )
            )
         else:
            self.bkgPdf = ROOT.RooRealFlooredSumPdf(
               PdfName, PdfName,
               ROOT.RooArgList(self.bkgTpl.theHistFunc),ROOT.RooArgList()
            )
      print self.bkgPdf.GetName(),"value =",self.bkgPdf.getVal()
      print self.bkgTpl.theRate.GetName(),"rate =",self.bkgTpl.theRate.getVal()



