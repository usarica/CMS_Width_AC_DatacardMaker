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

   def __init__(self, options, theMaker, theCategorizer, strBkgType, templateFileAppendName, systName):
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.channel = theMaker.channel
      self.workspace = theMaker.workspace
      self.theInputs = theMaker.theInputs

      # RooRealVars from the datacard maker class
      self.varm4l = theMaker.CMS_zz4l_widthMass
      self.varKD = theMaker.CMS_zz4l_widthKD
      self.varKD2 = theMaker.CMS_zz4l_widthKD2

      self.templateDir = options.templateDir
      self.dimensions = optionds.dimensions # Number of template dimensions>0
      self.ProjDim = options.ProjDim # The projected variable, -1 means do not project

      self.iCatScheme = theCategorizer.iCatScheme
      self.catNameList = theCategorizer.catNameList
      self.nCategories = theCategorizer.nCategories

      self.strBkgType = strBkgType
      self.condDim = None
      if(self.strBkgType == "ZX"):
         self.condDim = 0
      self.templateFileAppendName = templateFileAppendName
      self.systName = systName
      self.templateFileList = []

# Close the template files
   def close(self):
      for theFile in self.templateFileList:
         theFile.Close()

# Get shapes for each category
   def getTemplates(self):
      templateNameMain = "{0}/{1}".format(self.templateDir, self.templateFileAppendName)
      templateNameList = []
      for icat in range(0,self.nCategories):
         if(self.iCatScheme == 1): # icat==0: VBF, ==1: Non-VBF
            if(self.catNameList[icat] == ""):
               templateNameList.append("{0}{1}".format(templateNameMain,".root"))
            else:
               templateNameList.append("{0}{1}{2}{3}".format(templateNameMain,"_",self.catNameList[icat],".root"))
      for fname in templateNameList:
         print "Extracting template from {0}".format(fname)
         self.templateFileList.append(ROOT.TFile(fname, "read"))


#---------- TEMPLATES AND PDFS -------------

      self.bkg_T_2_list = []
      for icat in range(0,self.nCategories):
         self.bkg_T_2_list.append(
            ExtendedTemplate(
                  self.templateFileList[icat].Get("T_2D_{0}".format(self.strBkgType)).Clone("T_2D_{4}_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts,self.strBkgType)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2,
                  self.condDim
               )
            )

      self.bkgPdf_list = []
      self.bkgPdf_extras_list = []

      # qq bkg
      if(self.strBkgType == "qqZZ"):
         for icat in range(0,self.nCategories):
            PdfName = "qqZZ_OffshellPdf_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
            bkgPdf = ROOT.RooHistPdf(PdfName,PdfName,self.bkg_T_2_list[icat].argset,self.bkg_T_2_list[icat].theDataHist,0)
            self.bkgPdf_list.append(bkgPdf)

      # Z+X bkg
      elif(self.strBkgType == "ZX"):
         for icat in range(0,self.nCategories):
            PdfName = "zjets_OffshellPdf_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
            if self.ProjDim==0: # If projection on dim-0 is requested, just use the (unconditional) template already projected
               bkgPdf = ROOT.RooHistPdf(PdfName,PdfName,self.bkg_T_2_list[icat].argset,self.bkg_T_2_list[icat].theDataHist,0)
               self.bkgPdf_list.append(bkgPdf)
            else: # If projection on dim-0 is not requested, use the product of the mass pdf with the template conditional over dim-0
               pdfextraslist = []

               HistPdfName = "zjets_OffshellPdf_others_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
               bkgHistPdf = ROOT.RooHistPdf(HistPdfName,HistPdfName,self.bkg_T_2_list[icat].argset,self.bkg_T_2_list[icat].theDataHist,0)
               pdfextraslist.append(bkgHistPdf)

               val_mean_3P1F = float(self.theInputs['zjetsShape_mean_3P1F'])
               val_sigma_3P1F = float(self.theInputs['zjetsShape_sigma_3P1F'])
               val_norm_3P1F = float(self.theInputs['zjetsShape_norm_3P1F'])
               var_mean_3P1F = None
               var_sigma_3P1F = None
               var_norm_3P1F = None
               pdf_3P1F = None
               if val_norm_3P1F>0.:
                  tmpname = "zjets_OffshellPdf_mass_3p1f_mean_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_mean_3P1F = ROOT.RooRealVar(tmpname, "mean landau Zjet 3p1f", val_mean_3P1F)
                  tmpname = "zjets_OffshellPdf_mass_3p1f_sigma_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_sigma_3P1F = ROOT.RooRealVar(tmpname, "sigma landau Zjet 3p1f", val_sigma_3P1F)
                  tmpname = "zjets_OffshellPdf_mass_3p1f_norm_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_norm_3P1F = ROOT.RooRealVar(tmpname, "norm landau Zjet 3p1f", val_norm_3P1F)
                  pdf_3P1F = ROOT.RooLandau(
                     "zjets_OffshellPdf_mass_3p1f_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
                     "zjets_OffshellPdf_mass_3p1f_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
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
                  tmpname = "zjets_OffshellPdf_mass_2p2f_mean_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_mean_2P2F = ROOT.RooRealVar(tmpname, "mean landau Zjet 2p2f", val_mean_2P2F)
                  tmpname = "zjets_OffshellPdf_mass_2p2f_sigma_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_sigma_2P2F = ROOT.RooRealVar(tmpname, "sigma landau Zjet 2p2f", val_sigma_2P2F)
                  tmpname = "zjets_OffshellPdf_mass_2p2f_norm_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_norm_2P2F = ROOT.RooRealVar(tmpname, "norm landau Zjet 2p2f", val_norm_2P2F)
                  if val_pol1_2P2F==0.:
                     pdf_2P2F = ROOT.RooLandau(
                        "zjets_OffshellPdf_mass_2p2f_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
                        "zjets_OffshellPdf_mass_2p2f_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
                        self.varm4l,
                        var_mean_2P2F, var_sigma_2P2F
                     )
                  else:
                     tmpname = "zjets_OffshellPdf_mass_2p2f_pol0_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                     var_pol0_2P2F = ROOT.RooRealVar(tmpname, "pol0 landau Zjet 2p2f", val_pol0_2P2F)
                     tmpname = "zjets_OffshellPdf_mass_2p2f_pol1_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                     var_pol1_2P2F = ROOT.RooRealVar(tmpname, "pol1 landau Zjet 2p2f", val_pol1_2P2F)
                     pdf_2P2F = ROOT.RooGenericPdf(
                        "zjets_OffshellPdf_mass_2p2f_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
                        "zjets_OffshellPdf_mass_2p2f_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
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
                  tmpname = "zjets_OffshellPdf_mass_2p2f_2_mean_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_mean_2P2F_2 = ROOT.RooRealVar(tmpname, "mean landau Zjet 2p2f_2", val_mean_2P2F_2)
                  tmpname = "zjets_OffshellPdf_mass_2p2f_2_sigma_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_sigma_2P2F_2 = ROOT.RooRealVar(tmpname, "sigma landau Zjet 2p2f_2", val_sigma_2P2F_2)
                  tmpname = "zjets_OffshellPdf_mass_2p2f_2_norm_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
                  var_norm_2P2F_2 = ROOT.RooRealVar(tmpname, "norm landau Zjet 2p2f_2", val_norm_2P2F_2)
                  pdf_2P2F_2 = ROOT.RooLandau(
                     "zjets_OffshellPdf_mass_2p2f_2_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
                     "zjets_OffshellPdf_mass_2p2f_2_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts),
                     self.varm4l,
                     var_mean_2P2F_2, var_mean_2P2F_2
                  )

               MassPdfName = "zjets_OffshellPdf_mass_{0}_{1}_{2}_{3}TeV".format(self.systName,self.catNameList[icat],self.channel,self.sqrts)
               pdflist = ROOT.RooArgList()
               coeflist = ROOT.RooArgList()
               if ((pdf_3P1F is not None) and (var_norm_3P1F is not None)):
                  pdflist.add(pdf_3P1F)
                  coeflist.add(var_norm_3P1F)
                  pdfextraslist.append(pdf_3P1F)
               if ((pdf_2P2F is not None) and (var_norm_2P2F is not None)):
                  pdflist.add(pdf_2P2F)
                  coeflist.add(var_norm_2P2F)
                  pdfextraslist.append(pdf_2P2F)
               if ((pdf_2P2F_2 is not None) and (var_norm_2P2F_2 is not None)):
                  pdflist.add(pdf_2P2F_2)
                  coeflist.add(var_norm_2P2F_2)
                  pdfextraslist.append(pdf_2P2F_2)
               bkgMassPdf = ROOT.RooAddPdf(MassPdfName,MassPdfName,pdflist,coeflist)
               pdfextraslist.append(bkgMassPdf)

               self.bkgPdf_extras_list.append(pdfextraslist)

               bkgPdf = ROOT.RooProdPdf(
                  PdfName, PdfName,
                  ROOT.RooArgSet( bkgMassPdf ),
                  ROOT.RooFit.Conditional(
                      ROOT.RooArgSet( bkgHistPdf ),
                      ROOT.RooArgSet( self.varm4l )
                  )
               )
               self.bkgPdf_list.append(bkgPdf)



