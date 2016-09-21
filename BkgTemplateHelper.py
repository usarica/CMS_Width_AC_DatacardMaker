#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
import CategoryHelper
import ExtendedTemplate

class BkgTemplateHelper:
   def __init__(self, options, theMaker, theCategorizer, strBkgType, templateFileName, iCat, systName):
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.channel = theMaker.channel
      self.workspace = theMaker.workspace
      self.theInputs = theMaker.theInputs.getInputs()

      # RooRealVars from the datacard maker class
      self.varm4l = theEqnsMaker.varm4l
      self.varKD = theEqnsMaker.varKD
      self.varKD2 = theEqnsMaker.varKD2

      self.templateDir = options.templateDir
      self.dimensions = options.dimensions # Number of template dimensions>0
      self.ProjDim = options.ProjDim # The projected variable, -1 means do not project

      self.iCatScheme = theCategorizer.iCatScheme
      self.catNameList = theCategorizer.catNameList
      self.nCategories = theCategorizer.nCategories
      self.iCat = iCat
      if self.iCat>=self.nCategories:
         sys.exit("self.iCat={} >= self.nCategories={}!".format(self.iCat,self.nCategories))

      self.strBkgType = strBkgType
      self.condDim = None
      if((self.strBkgType.lower().startswith() == "zx") or (self.strBkgType.lower().startswith() == "zjets")):
         self.condDim = 0

      self.templateFileName = templateFileName
      self.systName = systName
      self.templateSuffix = "{0}_{1}_{2:.0f}_{3:.0f}TeV".format(self.systName,self.catNameList[self.iCat],self.channel,self.sqrts)

      self.templateFile = None

      self.bkg_T_2 = None
      self.bkgPdf = None
      self.bkgPdf_extras = []


# Import the pdf
   def importToWorkspace(self):
      if bkg_T_2 is not None:
         getattr(self.workspace, 'import')(self.bkg_T_2, ROOT.RooFit.RecycleConflictNodes())
      for xpdf in self.bkgPdf_extras:
         getattr(self.workspace, 'import')(xpdf, ROOT.RooFit.RecycleConflictNodes())
      if self.bkgPdf is not None:
         getattr(self.workspace, 'import')(self.bkgPdf, ROOT.RooFit.RecycleConflictNodes())


# Close the template files
   def close(self):
      self.templateFile.Close()


   def getThePdf(self):
      return self.bkgPdf
   def getTheRate(self):
      return self.bkg_T_2.theRate


# Get shapes for each category
   def getTemplates(self,templatePrefix="T_2D"):
      self.templateFile = ROOT.TFile.Open(self.templateFileName, "read")

      self.templatePrefix = templatePrefix
      self.templatePrefix = "{}_{}".format(self.templatePrefix,self.strBkgType)

#---------- TEMPLATES AND PDFS -------------
   # Construct the templates
      self.bkg_T_2 =
         ExtendedTemplate(
               self.templateFile.Get(self.templatePrefix).Clone("{}_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions, self.ProjDim,
               self.varm4l, self.varKD, self.varKD2,
               self.condDim
            )

   # Construct the p.d.f.s
      # qq bkg
      if(self.strBkgType.lower().startswith() == "qq"):
         PdfName = "qqZZ_OffshellPdf_{}".format(self.templateSuffix)
         self.bkgPdf = ROOT.RooHistPdf(PdfName,PdfName,self.bkg_T_2.argset,self.bkg_T_2.theDataHist,0)

      # Z+X bkg
      elif((self.strBkgType.lower().startswith() == "zx") or (self.strBkgType.lower().startswith() == "zjets")):
         PdfName = "zjets_OffshellPdf_{}".format(self.templateSuffix)
         if self.ProjDim==0: # If projection on dim-0 is requested, just use the (unconditional) template already projected
            self.bkgPdf = ROOT.RooHistPdf(PdfName,PdfName,self.bkg_T_2.argset,self.bkg_T_2.theDataHist,0)
         else: # If projection on dim-0 is not requested, use the product of the mass pdf with the template conditional over dim-0
            HistPdfName = "zjets_OffshellPdf_others_{}".format(self.templateSuffix)
            bkgHistPdf = ROOT.RooHistPdf(HistPdfName,HistPdfName,self.bkg_T_2.argset,self.bkg_T_2.theDataHist,0)
            self.bkgPdf_extras.append(bkgHistPdf)

            val_mean_3P1F = float(self.theInputs['zjetsShape_mean_3P1F'])
            val_sigma_3P1F = float(self.theInputs['zjetsShape_sigma_3P1F'])
            val_norm_3P1F = float(self.theInputs['zjetsShape_norm_3P1F'])
            var_mean_3P1F = None
            var_sigma_3P1F = None
            var_norm_3P1F = None
            pdf_3P1F = None
            if val_norm_3P1F>0.:
               tmpname = "zjets_OffshellPdf_mass_3p1f_mean_{}".format(self.templateSuffix)
               var_mean_3P1F = ROOT.RooRealVar(tmpname, "mean landau Zjet 3p1f", val_mean_3P1F)
               tmpname = "zjets_OffshellPdf_mass_3p1f_sigma_{}".format(self.templateSuffix)
               var_sigma_3P1F = ROOT.RooRealVar(tmpname, "sigma landau Zjet 3p1f", val_sigma_3P1F)
               tmpname = "zjets_OffshellPdf_mass_3p1f_norm_{}".format(self.templateSuffix)
               var_norm_3P1F = ROOT.RooRealVar(tmpname, "norm landau Zjet 3p1f", val_norm_3P1F)
               pdf_3P1F = ROOT.RooLandau(
                  "zjets_OffshellPdf_mass_3p1f_{}".format(self.templateSuffix),
                  "zjets_OffshellPdf_mass_3p1f_{}".format(self.templateSuffix),
                  self.varm4l,
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
               tmpname = "zjets_OffshellPdf_mass_2p2f_mean_{}".format(self.templateSuffix)
               var_mean_2P2F = ROOT.RooRealVar(tmpname, "mean landau Zjet 2p2f", val_mean_2P2F)
               tmpname = "zjets_OffshellPdf_mass_2p2f_sigma_{}".format(self.templateSuffix)
               var_sigma_2P2F = ROOT.RooRealVar(tmpname, "sigma landau Zjet 2p2f", val_sigma_2P2F)
               tmpname = "zjets_OffshellPdf_mass_2p2f_norm_{}".format(self.templateSuffix)
               var_norm_2P2F = ROOT.RooRealVar(tmpname, "norm landau Zjet 2p2f", val_norm_2P2F)
               if val_pol1_2P2F==0.:
                  pdf_2P2F = ROOT.RooLandau(
                     "zjets_OffshellPdf_mass_2p2f_{}".format(self.templateSuffix),
                     "zjets_OffshellPdf_mass_2p2f_{}".format(self.templateSuffix),
                     self.varm4l,
                     var_mean_2P2F, var_sigma_2P2F
                  )
               else:
                  tmpname = "zjets_OffshellPdf_mass_2p2f_pol0_{}".format(self.templateSuffix)
                  var_pol0_2P2F = ROOT.RooRealVar(tmpname, "pol0 landau Zjet 2p2f", val_pol0_2P2F)
                  tmpname = "zjets_OffshellPdf_mass_2p2f_pol1_{}".format(self.templateSuffix)
                  var_pol1_2P2F = ROOT.RooRealVar(tmpname, "pol1 landau Zjet 2p2f", val_pol1_2P2F)
                  pdf_2P2F = ROOT.RooGenericPdf(
                     "zjets_OffshellPdf_mass_2p2f_{}".format(self.templateSuffix),
                     "zjets_OffshellPdf_mass_2p2f_{}".format(self.templateSuffix),
                     "(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",
                     ROOT.RooArgList(
                        self.varm4l,
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
               tmpname = "zjets_OffshellPdf_mass_2p2f_2_mean_{}".format(self.templateSuffix)
               var_mean_2P2F_2 = ROOT.RooRealVar(tmpname, "mean landau Zjet 2p2f_2", val_mean_2P2F_2)
               tmpname = "zjets_OffshellPdf_mass_2p2f_2_sigma_{}".format(self.templateSuffix)
               var_sigma_2P2F_2 = ROOT.RooRealVar(tmpname, "sigma landau Zjet 2p2f_2", val_sigma_2P2F_2)
               tmpname = "zjets_OffshellPdf_mass_2p2f_2_norm_{}".format(self.templateSuffix)
               var_norm_2P2F_2 = ROOT.RooRealVar(tmpname, "norm landau Zjet 2p2f_2", val_norm_2P2F_2)
               pdf_2P2F_2 = ROOT.RooLandau(
                  "zjets_OffshellPdf_mass_2p2f_2_{}".format(self.templateSuffix),
                  "zjets_OffshellPdf_mass_2p2f_2_{}".format(self.templateSuffix),
                  self.varm4l,
                  var_mean_2P2F_2, var_mean_2P2F_2
               )

            MassPdfName = "zjets_OffshellPdf_mass_{}".format(self.templateSuffix)
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
                   ROOT.RooArgSet( self.varm4l )
               )
            )



