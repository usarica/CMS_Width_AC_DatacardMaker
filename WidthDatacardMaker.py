#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from InputCardReader import *
from CategoryHelper import *
from SystematicsHelper import *
from BSITemplateHelper import *
from BkgTemplateHelper import *


class WidthDatacardMaker:
   def __init__(self,options,theInputCard,theEqnsMaker,theCategorizer,theSystematizer,iCat,theOutputDir):
      self.options = options
      self.templateDir = self.options.templateDir
      self.dataAppendDir = self.options.dataDirAppend
      self.iCatScheme = self.options.iCatScheme
      self.mH = self.options.mPOLE
      self.low_M = self.options.mLow
      self.high_M = self.options.mHigh

      self.theInputCard = theInputCard
      self.theInputs = self.theInputCard.getInputs()
      self.sqrts = self.theInputCard.sqrts
      self.channel = self.theInputCard.channel

      self.theCategorizer = theCategorizer
      self.iCat = iCat
      self.theSystematizer = theSystematizer
      self.theOutputDir = theOutputDir

      self.workspace = ROOT.RooWorkspace("w", "w")

      # RooRealVars from the equations maker class
      self.muF = theEqnsMaker.muF
      self.muV = theEqnsMaker.muV
      self.kbkg_gg = theEqnsMaker.kbkg_gg
      self.kbkg_VBF = theEqnsMaker.kbkg_VBF
      self.fai1 = theEqnsMaker.fai1
      self.phiai1 = theEqnsMaker.phiai1
      self.phia1_gg = theEqnsMaker.phia1_gg
      self.phia1_VBF = theEqnsMaker.phia1_VBF
      self.varm4l = theEqnsMaker.varm4l
      self.varKD = theEqnsMaker.varKD
      self.varKD2 = theEqnsMaker.varKD2

      # LEFT HERE

      one = ROOT.RooRealVar("one", "one", 1.0)
      one.setConstant(True)

      self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts), "LUMI_{0:.0f}".format(self.sqrts), self.lumi)
      self.LUMI.setConstant(True)


      #----------------- Djet systematics ----------------------#

      CMS_zz4l_ggzz_djet_syst = w.factory("Djetscale_ggzz[-3,3]")
      CMS_zz4l_vbf_djet_syst = w.factory("Djetscale_vbf_offshell[-3,3]")
      CMS_zz4l_qqzz_djet_syst = w.factory("Djetscale_qqzz[-3,3]")
      CMS_zz4l_zjets_djet_syst = w.factory("Djetscale_zjets[-1,1]") # Djet>0.5 uncertainty is 99% for Z+X
      if useDjet == 0:
      CMS_zz4l_ggzz_djet_syst.setConstant(True)
      CMS_zz4l_vbf_djet_syst.setConstant(True)
      CMS_zz4l_qqzz_djet_syst.setConstant(True)
      CMS_zz4l_zjets_djet_syst.setConstant(True)

      morphVarListggZZ_djet = ROOT.RooArgList()
      morphVarListVBF_djet = ROOT.RooArgList()
      morphVarListqqZZ_djet = ROOT.RooArgList()
      morphVarListzjets_djet = ROOT.RooArgList()

      MorphNormList_ggZZ_Djet = ROOT.RooArgList()
      MorphNormList_VBF_Djet = ROOT.RooArgList()
      MorphNormList_qqZZ_Djet = ROOT.RooArgList()
      MorphNormList_zjets_Djet = ROOT.RooArgList()

      morphVarListggZZ_djet.add(CMS_zz4l_ggzz_djet_syst)
      ggzz_djetsyst_norm_nominal_name = "ggzz_djetsyst_norm_nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggzz_djetsyst_norm_nominal = RooRealVar(ggzz_djetsyst_norm_nominal_name, ggzz_djetsyst_norm_nominal_name, 1)
      ggzz_djetsyst_norm_up_name = "ggzz_djetsyst_norm_up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggzz_djetsyst_norm_up = RooRealVar(ggzz_djetsyst_norm_up_name, ggzz_djetsyst_norm_up_name, 1.+self.theInputs['djetscale_ggzz'])
      ggzz_djetsyst_norm_dn_name = "ggzz_djetsyst_norm_dn_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggzz_djetsyst_norm_dn = RooRealVar(ggzz_djetsyst_norm_dn_name, ggzz_djetsyst_norm_dn_name, 1.-self.theInputs['djetscale_ggzz'])
      MorphNormList_ggZZ_Djet.add(ggzz_djetsyst_norm_nominal)
      MorphNormList_ggZZ_Djet.add(ggzz_djetsyst_norm_up)
      MorphNormList_ggZZ_Djet.add(ggzz_djetsyst_norm_dn)
      asympowname = "Asymquad_ggZZ_djet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_djet_ggZZ_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_ggZZ_Djet, morphVarListggZZ_djet, 1.0)

      morphVarListVBF_djet.add(CMS_zz4l_vbf_djet_syst)
      VBF_djetsyst_norm_nominal_name = "VBF_djetsyst_norm_nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBF_djetsyst_norm_nominal = RooRealVar(VBF_djetsyst_norm_nominal_name, VBF_djetsyst_norm_nominal_name, 1)
      VBF_djetsyst_norm_up_name = "VBF_djetsyst_norm_up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBF_djetsyst_norm_up = RooRealVar(VBF_djetsyst_norm_up_name, VBF_djetsyst_norm_up_name, 1.+self.theInputs['djetscale_vbf_offshell'])
      VBF_djetsyst_norm_dn_name = "VBF_djetsyst_norm_dn_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBF_djetsyst_norm_dn = RooRealVar(VBF_djetsyst_norm_dn_name, VBF_djetsyst_norm_dn_name, 1.-self.theInputs['djetscale_vbf_offshell'])
      MorphNormList_VBF_Djet.add(VBF_djetsyst_norm_nominal)
      MorphNormList_VBF_Djet.add(VBF_djetsyst_norm_up)
      MorphNormList_VBF_Djet.add(VBF_djetsyst_norm_dn)
      asympowname = "Asymquad_VBF_djet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_djet_VBF_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_VBF_Djet, morphVarListVBF_djet, 1.0)

      morphVarListqqZZ_djet.add(CMS_zz4l_qqzz_djet_syst)
      qqZZ_djetsyst_norm_nominal_name = "qqZZ_djetsyst_norm_nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      qqZZ_djetsyst_norm_nominal = RooRealVar(qqZZ_djetsyst_norm_nominal_name, qqZZ_djetsyst_norm_nominal_name, 1)
      qqZZ_djetsyst_norm_up_name = "qqZZ_djetsyst_norm_up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      qqZZ_djetsyst_norm_up = RooRealVar(qqZZ_djetsyst_norm_up_name, qqZZ_djetsyst_norm_up_name, 1.+self.theInputs['djetscale_bkg_qqzz'])
      qqZZ_djetsyst_norm_dn_name = "qqZZ_djetsyst_norm_dn_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      qqZZ_djetsyst_norm_dn = RooRealVar(qqZZ_djetsyst_norm_dn_name, qqZZ_djetsyst_norm_dn_name, 1.-self.theInputs['djetscale_bkg_qqzz'])
      MorphNormList_qqZZ_Djet.add(qqZZ_djetsyst_norm_nominal)
      MorphNormList_qqZZ_Djet.add(qqZZ_djetsyst_norm_up)
      MorphNormList_qqZZ_Djet.add(qqZZ_djetsyst_norm_dn)
      asympowname = "Asymquad_qqZZ_djet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_djet_qqZZ_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_qqZZ_Djet, morphVarListqqZZ_djet, 1.0)
      if DEBUG:
      print "theta qqZZ djet nominal: ",thetaSyst_djet_qqZZ_norm.getVal()
      CMS_zz4l_qqzz_djet_syst.setVal(1)
      print "theta qqZZ djet up: ",thetaSyst_djet_qqZZ_norm.getVal()
      CMS_zz4l_qqzz_djet_syst.setVal(-1)
      print "theta qqZZ djet dn: ",thetaSyst_djet_qqZZ_norm.getVal()
      CMS_zz4l_qqzz_djet_syst.setVal(0)

      morphVarListzjets_djet.add(CMS_zz4l_zjets_djet_syst)
      zjets_djetsyst_norm_nominal_name = "zjets_djetsyst_norm_nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      zjets_djetsyst_norm_nominal = RooRealVar(zjets_djetsyst_norm_nominal_name, zjets_djetsyst_norm_nominal_name, 1)
      zjets_djetsyst_norm_up_name = "zjets_djetsyst_norm_up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      zjets_djetsyst_norm_up = RooRealVar(zjets_djetsyst_norm_up_name, zjets_djetsyst_norm_up_name, 1.+self.theInputs['djetscale_bkg_zjets'])
      zjets_djetsyst_norm_dn_name = "zjets_djetsyst_norm_dn_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      zjets_djetsyst_norm_dn = RooRealVar(zjets_djetsyst_norm_dn_name, zjets_djetsyst_norm_dn_name, 1.-self.theInputs['djetscale_bkg_zjets'])
      MorphNormList_zjets_Djet.add(zjets_djetsyst_norm_nominal)
      MorphNormList_zjets_Djet.add(zjets_djetsyst_norm_up)
      MorphNormList_zjets_Djet.add(zjets_djetsyst_norm_dn)
      asympowname = "Asymquad_zjets_djet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_djet_zjets_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_zjets_Djet, morphVarListzjets_djet, 1.0)


      #--------------------- End Djet systematics ---------------------------#


      # Add if-then for Djet cut here
      #-------
      print '2D signal shapes for Width'

      templateSigNameMain = "HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine_"
      if(self.anomCoupl == 1):
      templateSigNameMain = "{0}{1}".format(templateSigNameMain,"fLQAdded_")
      templateSigNameMain = "{0}{1}".format(templateSigNameMain,"_GenLevelVBF")
      if(USELEGACY == 1):
      templateSigNameMain = "{0}{1}".format(templateSigNameMain,"_wResolution")
      templateSigNameMain = "{0}{1}".format(templateSigNameMain,"_D_Gamma_gg_r10")
      templateSigNameMain = "{0}/LHC_{1:.0f}TeV/{2}/{3}".format(self.templateDir, self.sqrts, self.appendNameAlt,templateSigNameMain)

      templateSigNameMain_Nominal = "{0}_Nominal".format(templateSigNameMain)
      templateSigNameMainUp_PDF = "{0}_SysUp_ggPDF".format(templateSigNameMain)
      templateSigNameMainDown_PDF = "{0}_SysDown_ggPDF".format(templateSigNameMain)
      templateSigNameMainUp_QCD = "{0}_SysUp_ggQCD".format(templateSigNameMain)
      templateSigNameMainDown_QCD = "{0}_SysDown_ggQCD".format(templateSigNameMain)

      templateSigNameMain_Nominal_OppositeDjet = templateSigNameMain_Nominal
      templateSigNameMainUp_PDF_OppositeDjet = templateSigNameMainUp_PDF
      templateSigNameMainDown_PDF_OppositeDjet = templateSigNameMainDown_PDF
      templateSigNameMainUp_QCD_OppositeDjet = templateSigNameMainUp_QCD
      templateSigNameMainDown_QCD_OppositeDjet = templateSigNameMainDown_QCD



      if(useDjet == 0):
      templateSigNameMain_Nominal = "{0}{1}".format(templateSigNameMain_Nominal,".root")
      templateSigNameMainUp_PDF = "{0}{1}".format(templateSigNameMainUp_PDF,".root")
      templateSigNameMainDown_PDF = "{0}{1}".format(templateSigNameMainDown_PDF,".root")
      templateSigNameMainUp_QCD = "{0}{1}".format(templateSigNameMainUp_QCD,".root")
      templateSigNameMainDown_QCD = "{0}{1}".format(templateSigNameMainDown_QCD,".root")
      if(useDjet == 1):
      templateSigNameMain_Nominal = "{0}{1}".format(templateSigNameMain_Nominal,"_nonDjet.root")
      templateSigNameMainUp_PDF = "{0}{1}".format(templateSigNameMainUp_PDF,"_nonDjet.root")
      templateSigNameMainDown_PDF = "{0}{1}".format(templateSigNameMainDown_PDF,"_nonDjet.root")
      templateSigNameMainUp_QCD = "{0}{1}".format(templateSigNameMainUp_QCD,"_nonDjet.root")
      templateSigNameMainDown_QCD = "{0}{1}".format(templateSigNameMainDown_QCD,"_nonDjet.root")
      if(useDjet == 2):
      templateSigNameMain_Nominal = "{0}{1}".format(templateSigNameMain_Nominal,"_Djet.root")
      templateSigNameMainUp_PDF = "{0}{1}".format(templateSigNameMainUp_PDF,"_Djet.root")
      templateSigNameMainDown_PDF = "{0}{1}".format(templateSigNameMainDown_PDF,"_Djet.root")
      templateSigNameMainUp_QCD = "{0}{1}".format(templateSigNameMainUp_QCD,"_Djet.root")
      templateSigNameMainDown_QCD = "{0}{1}".format(templateSigNameMainDown_QCD,"_Djet.root")

      if useDjet == 0:
      templateSigNameMain_Nominal_OppositeDjet = "{0}{1}".format(templateSigNameMain_Nominal_OppositeDjet,"_Djet.root")
      templateSigNameMainUp_PDF_OppositeDjet = "{0}{1}".format(templateSigNameMainUp_PDF_OppositeDjet,"_Djet.root")
      templateSigNameMainDown_PDF_OppositeDjet = "{0}{1}".format(templateSigNameMainDown_PDF_OppositeDjet,"_Djet.root")
      templateSigNameMainUp_QCD_OppositeDjet = "{0}{1}".format(templateSigNameMainUp_QCD_OppositeDjet,"_Djet.root")
      templateSigNameMainDown_QCD_OppositeDjet = "{0}{1}".format(templateSigNameMainDown_QCD_OppositeDjet,"_Djet.root")
      if useDjet == 1:
      templateSigNameMain_Nominal_OppositeDjet = "{0}{1}".format(templateSigNameMain_Nominal_OppositeDjet,"_Djet.root")
      templateSigNameMainUp_PDF_OppositeDjet = "{0}{1}".format(templateSigNameMainUp_PDF_OppositeDjet,"_Djet.root")
      templateSigNameMainDown_PDF_OppositeDjet = "{0}{1}".format(templateSigNameMainDown_PDF_OppositeDjet,"_Djet.root")
      templateSigNameMainUp_QCD_OppositeDjet = "{0}{1}".format(templateSigNameMainUp_QCD_OppositeDjet,"_Djet.root")
      templateSigNameMainDown_QCD_OppositeDjet = "{0}{1}".format(templateSigNameMainDown_QCD_OppositeDjet,"_Djet.root")
      if useDjet == 2:
      templateSigNameMain_Nominal_OppositeDjet = "{0}{1}".format(templateSigNameMain_Nominal_OppositeDjet,"_nonDjet.root")
      templateSigNameMainUp_PDF_OppositeDjet = "{0}{1}".format(templateSigNameMainUp_PDF_OppositeDjet,"_nonDjet.root")
      templateSigNameMainDown_PDF_OppositeDjet = "{0}{1}".format(templateSigNameMainDown_PDF_OppositeDjet,"_nonDjet.root")
      templateSigNameMainUp_QCD_OppositeDjet = "{0}{1}".format(templateSigNameMainUp_QCD_OppositeDjet,"_nonDjet.root")
      templateSigNameMainDown_QCD_OppositeDjet = "{0}{1}".format(templateSigNameMainDown_QCD_OppositeDjet,"_nonDjet.root")

      sigTempFileU = ROOT.TFile(templateSigNameMain_Nominal)
      sigTempFileUp_PDF = ROOT.TFile(templateSigNameMainUp_PDF)
      sigTempFileDown_PDF = ROOT.TFile(templateSigNameMainDown_PDF)
      sigTempFileUp_QCD = ROOT.TFile(templateSigNameMainUp_QCD)
      sigTempFileDown_QCD = ROOT.TFile(templateSigNameMainDown_QCD)

      sigTempFileU_OppositeDjet = ROOT.TFile(templateSigNameMain_Nominal_OppositeDjet)
      sigTempFileUp_PDF_OppositeDjet = ROOT.TFile(templateSigNameMainUp_PDF_OppositeDjet)
      sigTempFileDown_PDF_OppositeDjet = ROOT.TFile(templateSigNameMainDown_PDF_OppositeDjet)
      sigTempFileUp_QCD_OppositeDjet = ROOT.TFile(templateSigNameMainUp_QCD_OppositeDjet)
      sigTempFileDown_QCD_OppositeDjet = ROOT.TFile(templateSigNameMainDown_QCD_OppositeDjet)


      #---------- SIGNAL TEMPLATES -------------

      # Bare SM

      Sig_T_2 = sigTempFileU.Get("T_2D_1").Clone("T_2D_1_Nominal")
      Sig_T_2_Up_PDF = sigTempFileUp_PDF.Get("T_2D_1").Clone("T_2D_1_PDFUp")
      Sig_T_2_Up_QCD = sigTempFileUp_QCD.Get("T_2D_1").Clone("T_2D_1_QCDUp")
      Sig_T_2_Down_PDF = sigTempFileDown_PDF.Get("T_2D_1").Clone("T_2D_1_PDFDown")
      Sig_T_2_Down_QCD = sigTempFileDown_QCD.Get("T_2D_1").Clone("T_2D_1_QCDDown")

      Sig_T_1 = sigTempFileU.Get("T_2D_2").Clone("T_2D_2_Nominal")
      Sig_T_1_Up_PDF = sigTempFileUp_PDF.Get("T_2D_2").Clone("T_2D_2_PDFUp")
      Sig_T_1_Down_PDF = sigTempFileDown_PDF.Get("T_2D_2").Clone("T_2D_2_PDFDown")
      Sig_T_1_Up_QCD = sigTempFileUp_QCD.Get("T_2D_2").Clone("T_2D_2_QCDUp")
      Sig_T_1_Down_QCD = sigTempFileDown_QCD.Get("T_2D_2").Clone("T_2D_2_QCDDown")

      Sig_T_4 = sigTempFileU.Get("T_2D_4").Clone("T_2D_4_Nominal")
      Sig_T_4_Up_PDF = sigTempFileUp_PDF.Get("T_2D_4").Clone("T_2D_4_PDFUp")
      Sig_T_4_Up_QCD = sigTempFileUp_QCD.Get("T_2D_4").Clone("T_2D_4_QCDUp")
      Sig_T_4_Down_PDF = sigTempFileDown_PDF.Get("T_2D_4").Clone("T_2D_4_PDFDown")
      Sig_T_4_Down_QCD = sigTempFileDown_QCD.Get("T_2D_4").Clone("T_2D_4_QCDDown")

      VBF_T_2 = sigTempFileU.Get("T_2D_VBF_1").Clone("T_2D_VBF_1_Nominal")
      VBF_T_2_Up = sigTempFileUp_PDF.Get("T_2D_VBF_1").Clone("T_2D_VBF_1_PDFUp")
      VBF_T_2_Down = sigTempFileDown_PDF.Get("T_2D_VBF_1").Clone("T_2D_VBF_1_PDFDown")

      VBF_T_1 = sigTempFileU.Get("T_2D_VBF_2").Clone("T_2D_VBF_2_Nominal")
      VBF_T_1_Up = sigTempFileUp_PDF.Get("T_2D_VBF_2").Clone("T_2D_VBF_2_PDFUp")
      VBF_T_1_Down = sigTempFileDown_PDF.Get("T_2D_VBF_2").Clone("T_2D_VBF_2_PDFDown")

      VBF_T_4 = sigTempFileU.Get("T_2D_VBF_4").Clone("T_2D_VBF_4_Nominal")
      VBF_T_4_Up = sigTempFileUp_PDF.Get("T_2D_VBF_4").Clone("T_2D_VBF_4_PDFUp")
      VBF_T_4_Down = sigTempFileDown_PDF.Get("T_2D_VBF_4").Clone("T_2D_VBF_4_PDFDown")



      Sig_T_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_1").Clone("T_2D_1_Nominal_OppositeDjet")
      Sig_T_2_Up_PDF_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_1").Clone("T_2D_1_PDFUp_OppositeDjet")
      Sig_T_2_Up_QCD_OppositeDjet = sigTempFileUp_QCD_OppositeDjet.Get("T_2D_1").Clone("T_2D_1_QCDUp_OppositeDjet")
      Sig_T_2_Down_PDF_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_1").Clone("T_2D_1_PDFDown_OppositeDjet")
      Sig_T_2_Down_QCD_OppositeDjet = sigTempFileDown_QCD_OppositeDjet.Get("T_2D_1").Clone("T_2D_1_QCDDown_OppositeDjet")

      Sig_T_1_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_2").Clone("T_2D_2_Nominal_OppositeDjet")
      Sig_T_1_Up_PDF_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_2").Clone("T_2D_2_PDFUp_OppositeDjet")
      Sig_T_1_Down_PDF_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_2").Clone("T_2D_2_PDFDown_OppositeDjet")
      Sig_T_1_Up_QCD_OppositeDjet = sigTempFileUp_QCD_OppositeDjet.Get("T_2D_2").Clone("T_2D_2_QCDUp_OppositeDjet")
      Sig_T_1_Down_QCD_OppositeDjet = sigTempFileDown_QCD_OppositeDjet.Get("T_2D_2").Clone("T_2D_2_QCDDown_OppositeDjet")

      Sig_T_4_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_4").Clone("T_2D_4_Nominal_OppositeDjet")
      Sig_T_4_Up_PDF_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_4").Clone("T_2D_4_PDFUp_OppositeDjet")
      Sig_T_4_Up_QCD_OppositeDjet = sigTempFileUp_QCD_OppositeDjet.Get("T_2D_4").Clone("T_2D_4_QCDUp_OppositeDjet")
      Sig_T_4_Down_PDF_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_4").Clone("T_2D_4_PDFDown_OppositeDjet")
      Sig_T_4_Down_QCD_OppositeDjet = sigTempFileDown_QCD_OppositeDjet.Get("T_2D_4").Clone("T_2D_4_QCDDown_OppositeDjet")

      VBF_T_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_1").Clone("T_2D_VBF_1_Nominal_OppositeDjet")
      VBF_T_2_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_1").Clone("T_2D_VBF_1_PDFUp_OppositeDjet")
      VBF_T_2_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_1").Clone("T_2D_VBF_1_PDFDown_OppositeDjet")

      VBF_T_1_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_2").Clone("T_2D_VBF_2_Nominal_OppositeDjet")
      VBF_T_1_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_2").Clone("T_2D_VBF_2_PDFUp_OppositeDjet")
      VBF_T_1_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_2").Clone("T_2D_VBF_2_PDFDown_OppositeDjet")

      VBF_T_4_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_4").Clone("T_2D_VBF_4_Nominal_OppositeDjet")
      VBF_T_4_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_4").Clone("T_2D_VBF_4_PDFUp_OppositeDjet")
      VBF_T_4_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_4").Clone("T_2D_VBF_4_PDFDown_OppositeDjet")



      integral_Sig_T_1 = Sig_T_1.Integral("width")
      integral_Sig_T_1_Up_PDF = Sig_T_1_Up_PDF.Integral("width")
      integral_Sig_T_1_Down_PDF = Sig_T_1_Down_PDF.Integral("width")
      integral_Sig_T_1_Up_QCD = Sig_T_1_Up_QCD.Integral("width")
      integral_Sig_T_1_Down_QCD = Sig_T_1_Down_QCD.Integral("width")
      integral_Sig_T_2 = Sig_T_2.Integral("width")
      integral_Sig_T_2_Up_PDF = Sig_T_2_Up_PDF.Integral("width")
      integral_Sig_T_2_Down_PDF = Sig_T_2_Down_PDF.Integral("width")
      integral_Sig_T_2_Up_QCD = Sig_T_2_Up_QCD.Integral("width")
      integral_Sig_T_2_Down_QCD = Sig_T_2_Down_QCD.Integral("width")
      integral_Sig_T_4 = Sig_T_4.Integral("width")
      integral_Sig_T_4_Up_PDF = Sig_T_4_Up_PDF.Integral("width")
      integral_Sig_T_4_Down_PDF = Sig_T_4_Down_PDF.Integral("width")
      integral_Sig_T_4_Up_QCD = Sig_T_4_Up_QCD.Integral("width")
      integral_Sig_T_4_Down_QCD = Sig_T_4_Down_QCD.Integral("width")

      integral_VBF_T_1 = VBF_T_1.Integral("width")
      integral_VBF_T_1_Up = VBF_T_1_Up.Integral("width")
      integral_VBF_T_1_Down = VBF_T_1_Down.Integral("width")
      integral_VBF_T_2 = VBF_T_2.Integral("width")
      integral_VBF_T_2_Up = VBF_T_2_Up.Integral("width")
      integral_VBF_T_2_Down = VBF_T_2_Down.Integral("width")
      integral_VBF_T_4 = VBF_T_4.Integral("width")
      integral_VBF_T_4_Up = VBF_T_4_Up.Integral("width")
      integral_VBF_T_4_Down = VBF_T_4_Down.Integral("width")


      integral_Sig_T_1_OppositeDjet = Sig_T_1_OppositeDjet.Integral("width")
      integral_Sig_T_1_Up_PDF_OppositeDjet = Sig_T_1_Up_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_1_Down_PDF_OppositeDjet = Sig_T_1_Down_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_1_Up_QCD_OppositeDjet = Sig_T_1_Up_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_1_Down_QCD_OppositeDjet = Sig_T_1_Down_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_2_OppositeDjet = Sig_T_2_OppositeDjet.Integral("width")
      integral_Sig_T_2_Up_PDF_OppositeDjet = Sig_T_2_Up_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_2_Down_PDF_OppositeDjet = Sig_T_2_Down_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_2_Up_QCD_OppositeDjet = Sig_T_2_Up_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_2_Down_QCD_OppositeDjet = Sig_T_2_Down_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_4_OppositeDjet = Sig_T_4_OppositeDjet.Integral("width")
      integral_Sig_T_4_Up_PDF_OppositeDjet = Sig_T_4_Up_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_4_Down_PDF_OppositeDjet = Sig_T_4_Down_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_4_Up_QCD_OppositeDjet = Sig_T_4_Up_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_4_Down_QCD_OppositeDjet = Sig_T_4_Down_QCD_OppositeDjet.Integral("width")

      integral_VBF_T_1_OppositeDjet = VBF_T_1_OppositeDjet.Integral("width")
      integral_VBF_T_1_Up_OppositeDjet = VBF_T_1_Up_OppositeDjet.Integral("width")
      integral_VBF_T_1_Down_OppositeDjet = VBF_T_1_Down_OppositeDjet.Integral("width")
      integral_VBF_T_2_OppositeDjet = VBF_T_2_OppositeDjet.Integral("width")
      integral_VBF_T_2_Up_OppositeDjet = VBF_T_2_Up_OppositeDjet.Integral("width")
      integral_VBF_T_2_Down_OppositeDjet = VBF_T_2_Down_OppositeDjet.Integral("width")
      integral_VBF_T_4_OppositeDjet = VBF_T_4_OppositeDjet.Integral("width")
      integral_VBF_T_4_Up_OppositeDjet = VBF_T_4_Up_OppositeDjet.Integral("width")
      integral_VBF_T_4_Down_OppositeDjet = VBF_T_4_Down_OppositeDjet.Integral("width")



      # (mZZ/mH)**2 Terms

      Sig_T_mZZ2_1_2 = sigTempFileU.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_Nominal")
      Sig_T_mZZ2_1_2_Up_PDF = sigTempFileUp_PDF.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_PDFUp")
      Sig_T_mZZ2_1_2_Up_QCD = sigTempFileUp_QCD.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_QCDUp")
      Sig_T_mZZ2_1_2_Down_PDF = sigTempFileDown_PDF.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_PDFDown")
      Sig_T_mZZ2_1_2_Down_QCD = sigTempFileDown_QCD.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_QCDDown")

      Sig_T_mZZ2_1_4 = sigTempFileU.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_Nominal")
      Sig_T_mZZ2_1_4_Up_PDF = sigTempFileUp_PDF.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_PDFUp")
      Sig_T_mZZ2_1_4_Up_QCD = sigTempFileUp_QCD.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_QCDUp")
      Sig_T_mZZ2_1_4_Down_PDF = sigTempFileDown_PDF.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_PDFDown")
      Sig_T_mZZ2_1_4_Down_QCD = sigTempFileDown_QCD.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_QCDDown")

      VBF_T_mZZ2_1_2 = sigTempFileU.Get("T_2D_VBF_1_mZZ2_1").Clone("T_2D_VBF_1_mZZ2_1_Nominal")
      VBF_T_mZZ2_1_2_Up = sigTempFileUp_PDF.Get("T_2D_VBF_1_mZZ2_1").Clone("T_2D_VBF_1_mZZ2_1_PDFUp")
      VBF_T_mZZ2_1_2_Down = sigTempFileDown_PDF.Get("T_2D_VBF_1_mZZ2_1").Clone("T_2D_VBF_1_mZZ2_1_PDFDown")

      VBF_T_mZZ2_1_4 = sigTempFileU.Get("T_2D_VBF_4_mZZ2_1").Clone("T_2D_VBF_4_mZZ2_1_Nominal")
      VBF_T_mZZ2_1_4_Up = sigTempFileUp_PDF.Get("T_2D_VBF_4_mZZ2_1").Clone("T_2D_VBF_4_mZZ2_1_PDFUp")
      VBF_T_mZZ2_1_4_Down = sigTempFileDown_PDF.Get("T_2D_VBF_4_mZZ2_1").Clone("T_2D_VBF_4_mZZ2_1_PDFDown")


      Sig_T_mZZ2_1_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_Nominal_OppositeDjet")
      Sig_T_mZZ2_1_2_Up_PDF_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_PDFUp_OppositeDjet")
      Sig_T_mZZ2_1_2_Up_QCD_OppositeDjet = sigTempFileUp_QCD_OppositeDjet.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_QCDUp_OppositeDjet")
      Sig_T_mZZ2_1_2_Down_PDF_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_PDFDown_OppositeDjet")
      Sig_T_mZZ2_1_2_Down_QCD_OppositeDjet = sigTempFileDown_QCD_OppositeDjet.Get("T_2D_1_mZZ2_1").Clone("T_2D_1_mZZ2_1_QCDDown_OppositeDjet")

      Sig_T_mZZ2_1_4_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_Nominal_OppositeDjet")
      Sig_T_mZZ2_1_4_Up_PDF_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_PDFUp_OppositeDjet")
      Sig_T_mZZ2_1_4_Up_QCD_OppositeDjet = sigTempFileUp_QCD_OppositeDjet.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_QCDUp_OppositeDjet")
      Sig_T_mZZ2_1_4_Down_PDF_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_PDFDown_OppositeDjet")
      Sig_T_mZZ2_1_4_Down_QCD_OppositeDjet = sigTempFileDown_QCD_OppositeDjet.Get("T_2D_4_mZZ2_1").Clone("T_2D_4_mZZ2_1_QCDDown_OppositeDjet")

      VBF_T_mZZ2_1_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_1_mZZ2_1").Clone("T_2D_VBF_1_mZZ2_1_Nominal_OppositeDjet")
      VBF_T_mZZ2_1_2_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_1").Clone("T_2D_VBF_1_mZZ2_1_PDFUp_OppositeDjet")
      VBF_T_mZZ2_1_2_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_1").Clone("T_2D_VBF_1_mZZ2_1_PDFDown_OppositeDjet")

      VBF_T_mZZ2_1_4_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_4_mZZ2_1").Clone("T_2D_VBF_4_mZZ2_1_Nominal_OppositeDjet")
      VBF_T_mZZ2_1_4_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_4_mZZ2_1").Clone("T_2D_VBF_4_mZZ2_1_PDFUp_OppositeDjet")
      VBF_T_mZZ2_1_4_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_4_mZZ2_1").Clone("T_2D_VBF_4_mZZ2_1_PDFDown_OppositeDjet")


      integral_Sig_T_mZZ2_1_2 = Sig_T_mZZ2_1_2.Integral("width")
      integral_Sig_T_mZZ2_1_2_Up_PDF = Sig_T_mZZ2_1_2_Up_PDF.Integral("width")
      integral_Sig_T_mZZ2_1_2_Down_PDF = Sig_T_mZZ2_1_2_Down_PDF.Integral("width")
      integral_Sig_T_mZZ2_1_2_Up_QCD = Sig_T_mZZ2_1_2_Up_QCD.Integral("width")
      integral_Sig_T_mZZ2_1_2_Down_QCD = Sig_T_mZZ2_1_2_Down_QCD.Integral("width")

      integral_Sig_T_mZZ2_1_4 = Sig_T_mZZ2_1_4.Integral("width")
      integral_Sig_T_mZZ2_1_4_Up_PDF = Sig_T_mZZ2_1_4_Up_PDF.Integral("width")
      integral_Sig_T_mZZ2_1_4_Down_PDF = Sig_T_mZZ2_1_4_Down_PDF.Integral("width")
      integral_Sig_T_mZZ2_1_4_Up_QCD = Sig_T_mZZ2_1_4_Up_QCD.Integral("width")
      integral_Sig_T_mZZ2_1_4_Down_QCD = Sig_T_mZZ2_1_4_Down_QCD.Integral("width")

      integral_VBF_T_mZZ2_1_2 = VBF_T_mZZ2_1_2.Integral("width")
      integral_VBF_T_mZZ2_1_2_Up = VBF_T_mZZ2_1_2_Up.Integral("width")
      integral_VBF_T_mZZ2_1_2_Down = VBF_T_mZZ2_1_2_Down.Integral("width")

      integral_VBF_T_mZZ2_1_4 = VBF_T_mZZ2_1_4.Integral("width")
      integral_VBF_T_mZZ2_1_4_Up = VBF_T_mZZ2_1_4_Up.Integral("width")
      integral_VBF_T_mZZ2_1_4_Down = VBF_T_mZZ2_1_4_Down.Integral("width")


      integral_Sig_T_mZZ2_1_2_OppositeDjet = Sig_T_mZZ2_1_2_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_2_Up_PDF_OppositeDjet = Sig_T_mZZ2_1_2_Up_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_2_Down_PDF_OppositeDjet = Sig_T_mZZ2_1_2_Down_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_2_Up_QCD_OppositeDjet = Sig_T_mZZ2_1_2_Up_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_2_Down_QCD_OppositeDjet = Sig_T_mZZ2_1_2_Down_QCD_OppositeDjet.Integral("width")

      integral_Sig_T_mZZ2_1_4_OppositeDjet = Sig_T_mZZ2_1_4_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_4_Up_PDF_OppositeDjet = Sig_T_mZZ2_1_4_Up_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_4_Down_PDF_OppositeDjet = Sig_T_mZZ2_1_4_Down_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_4_Up_QCD_OppositeDjet = Sig_T_mZZ2_1_4_Up_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_1_4_Down_QCD_OppositeDjet = Sig_T_mZZ2_1_4_Down_QCD_OppositeDjet.Integral("width")

      integral_VBF_T_mZZ2_1_2_OppositeDjet = VBF_T_mZZ2_1_2_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_1_2_Up_OppositeDjet = VBF_T_mZZ2_1_2_Up_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_1_2_Down_OppositeDjet = VBF_T_mZZ2_1_2_Down_OppositeDjet.Integral("width")

      integral_VBF_T_mZZ2_1_4_OppositeDjet = VBF_T_mZZ2_1_4_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_1_4_Up_OppositeDjet = VBF_T_mZZ2_1_4_Up_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_1_4_Down_OppositeDjet = VBF_T_mZZ2_1_4_Down_OppositeDjet.Integral("width")


      # (mZZ/mH)**2**2 Terms

      Sig_T_mZZ2_2_2 = sigTempFileU.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_Nominal")
      Sig_T_mZZ2_2_2_Up_PDF = sigTempFileUp_PDF.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_PDFUp")
      Sig_T_mZZ2_2_2_Up_QCD = sigTempFileUp_QCD.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_QCDUp")
      Sig_T_mZZ2_2_2_Down_PDF = sigTempFileDown_PDF.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_PDFDown")
      Sig_T_mZZ2_2_2_Down_QCD = sigTempFileDown_QCD.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_QCDDown")

      VBF_T_mZZ2_2_2 = sigTempFileU.Get("T_2D_VBF_1_mZZ2_2").Clone("T_2D_VBF_1_mZZ2_2_Nominal")
      VBF_T_mZZ2_2_2_Up = sigTempFileUp_PDF.Get("T_2D_VBF_1_mZZ2_2").Clone("T_2D_VBF_1_mZZ2_2_PDFUp")
      VBF_T_mZZ2_2_2_Down = sigTempFileDown_PDF.Get("T_2D_VBF_1_mZZ2_2").Clone("T_2D_VBF_1_mZZ2_2_PDFDown")

      VBF_T_mZZ2_2_4 = sigTempFileU.Get("T_2D_VBF_4_mZZ2_2").Clone("T_2D_VBF_4_mZZ2_2_Nominal")
      VBF_T_mZZ2_2_4_Up = sigTempFileUp_PDF.Get("T_2D_VBF_4_mZZ2_2").Clone("T_2D_VBF_4_mZZ2_2_PDFUp")
      VBF_T_mZZ2_2_4_Down = sigTempFileDown_PDF.Get("T_2D_VBF_4_mZZ2_2").Clone("T_2D_VBF_4_mZZ2_2_PDFDown")


      Sig_T_mZZ2_2_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_Nominal_OppositeDjet")
      Sig_T_mZZ2_2_2_Up_PDF_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_PDFUp_OppositeDjet")
      Sig_T_mZZ2_2_2_Up_QCD_OppositeDjet = sigTempFileUp_QCD_OppositeDjet.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_QCDUp_OppositeDjet")
      Sig_T_mZZ2_2_2_Down_PDF_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_PDFDown_OppositeDjet")
      Sig_T_mZZ2_2_2_Down_QCD_OppositeDjet = sigTempFileDown_QCD_OppositeDjet.Get("T_2D_1_mZZ2_2").Clone("T_2D_1_mZZ2_2_QCDDown_OppositeDjet")

      VBF_T_mZZ2_2_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_1_mZZ2_2").Clone("T_2D_VBF_1_mZZ2_2_Nominal_OppositeDjet")
      VBF_T_mZZ2_2_2_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_2").Clone("T_2D_VBF_1_mZZ2_2_PDFUp_OppositeDjet")
      VBF_T_mZZ2_2_2_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_2").Clone("T_2D_VBF_1_mZZ2_2_PDFDown_OppositeDjet")

      VBF_T_mZZ2_2_4_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_4_mZZ2_2").Clone("T_2D_VBF_4_mZZ2_2_Nominal_OppositeDjet")
      VBF_T_mZZ2_2_4_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_4_mZZ2_2").Clone("T_2D_VBF_4_mZZ2_2_PDFUp_OppositeDjet")
      VBF_T_mZZ2_2_4_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_4_mZZ2_2").Clone("T_2D_VBF_4_mZZ2_2_PDFDown_OppositeDjet")



      integral_Sig_T_mZZ2_2_2 = Sig_T_mZZ2_2_2.Integral("width")
      integral_Sig_T_mZZ2_2_2_Up_PDF = Sig_T_mZZ2_2_2_Up_PDF.Integral("width")
      integral_Sig_T_mZZ2_2_2_Down_PDF = Sig_T_mZZ2_2_2_Down_PDF.Integral("width")
      integral_Sig_T_mZZ2_2_2_Up_QCD = Sig_T_mZZ2_2_2_Up_QCD.Integral("width")
      integral_Sig_T_mZZ2_2_2_Down_QCD = Sig_T_mZZ2_2_2_Down_QCD.Integral("width")

      integral_VBF_T_mZZ2_2_2 = VBF_T_mZZ2_2_2.Integral("width")
      integral_VBF_T_mZZ2_2_2_Up = VBF_T_mZZ2_2_2_Up.Integral("width")
      integral_VBF_T_mZZ2_2_2_Down = VBF_T_mZZ2_2_2_Down.Integral("width")

      integral_VBF_T_mZZ2_2_4 = VBF_T_mZZ2_2_4.Integral("width")
      integral_VBF_T_mZZ2_2_4_Up = VBF_T_mZZ2_2_4_Up.Integral("width")
      integral_VBF_T_mZZ2_2_4_Down = VBF_T_mZZ2_2_4_Down.Integral("width")


      integral_Sig_T_mZZ2_2_2_OppositeDjet = Sig_T_mZZ2_2_2_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_2_2_Up_PDF_OppositeDjet = Sig_T_mZZ2_2_2_Up_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_2_2_Down_PDF_OppositeDjet = Sig_T_mZZ2_2_2_Down_PDF_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_2_2_Up_QCD_OppositeDjet = Sig_T_mZZ2_2_2_Up_QCD_OppositeDjet.Integral("width")
      integral_Sig_T_mZZ2_2_2_Down_QCD_OppositeDjet = Sig_T_mZZ2_2_2_Down_QCD_OppositeDjet.Integral("width")

      integral_VBF_T_mZZ2_2_2_OppositeDjet = VBF_T_mZZ2_2_2_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_2_2_Up_OppositeDjet = VBF_T_mZZ2_2_2_Up_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_2_2_Down_OppositeDjet = VBF_T_mZZ2_2_2_Down_OppositeDjet.Integral("width")

      integral_VBF_T_mZZ2_2_4_OppositeDjet = VBF_T_mZZ2_2_4_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_2_4_Up_OppositeDjet = VBF_T_mZZ2_2_4_Up_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_2_4_Down_OppositeDjet = VBF_T_mZZ2_2_4_Down_OppositeDjet.Integral("width")


      # (mZZ/mH)**2**3 Terms


      VBF_T_mZZ2_3_2 = sigTempFileU.Get("T_2D_VBF_1_mZZ2_3").Clone("T_2D_VBF_1_mZZ2_3_Nominal")
      VBF_T_mZZ2_3_2_Up = sigTempFileUp_PDF.Get("T_2D_VBF_1_mZZ2_3").Clone("T_2D_VBF_1_mZZ2_3_PDFUp")
      VBF_T_mZZ2_3_2_Down = sigTempFileDown_PDF.Get("T_2D_VBF_1_mZZ2_3").Clone("T_2D_VBF_1_mZZ2_3_PDFDown")

      VBF_T_mZZ2_3_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_1_mZZ2_3").Clone("T_2D_VBF_1_mZZ2_3_Nominal_OppositeDjet")
      VBF_T_mZZ2_3_2_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_3").Clone("T_2D_VBF_1_mZZ2_3_PDFUp_OppositeDjet")
      VBF_T_mZZ2_3_2_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_3").Clone("T_2D_VBF_1_mZZ2_3_PDFDown_OppositeDjet")


      integral_VBF_T_mZZ2_3_2 = VBF_T_mZZ2_3_2.Integral("width")
      integral_VBF_T_mZZ2_3_2_Up = VBF_T_mZZ2_3_2_Up.Integral("width")
      integral_VBF_T_mZZ2_3_2_Down = VBF_T_mZZ2_3_2_Down.Integral("width")

      integral_VBF_T_mZZ2_3_2_OppositeDjet = VBF_T_mZZ2_3_2_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_3_2_Up_OppositeDjet = VBF_T_mZZ2_3_2_Up_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_3_2_Down_OppositeDjet = VBF_T_mZZ2_3_2_Down_OppositeDjet.Integral("width")


      # (mZZ/mH)**2**4 Terms


      VBF_T_mZZ2_4_2 = sigTempFileU.Get("T_2D_VBF_1_mZZ2_4").Clone("T_2D_VBF_1_mZZ2_4_Nominal")
      VBF_T_mZZ2_4_2_Up = sigTempFileUp_PDF.Get("T_2D_VBF_1_mZZ2_4").Clone("T_2D_VBF_1_mZZ2_4_PDFUp")
      VBF_T_mZZ2_4_2_Down = sigTempFileDown_PDF.Get("T_2D_VBF_1_mZZ2_4").Clone("T_2D_VBF_1_mZZ2_4_PDFDown")

      VBF_T_mZZ2_4_2_OppositeDjet = sigTempFileU_OppositeDjet.Get("T_2D_VBF_1_mZZ2_4").Clone("T_2D_VBF_1_mZZ2_4_Nominal_OppositeDjet")
      VBF_T_mZZ2_4_2_Up_OppositeDjet = sigTempFileUp_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_4").Clone("T_2D_VBF_1_mZZ2_4_PDFUp_OppositeDjet")
      VBF_T_mZZ2_4_2_Down_OppositeDjet = sigTempFileDown_PDF_OppositeDjet.Get("T_2D_VBF_1_mZZ2_4").Clone("T_2D_VBF_1_mZZ2_4_PDFDown_OppositeDjet")


      integral_VBF_T_mZZ2_4_2 = VBF_T_mZZ2_4_2.Integral("width")
      integral_VBF_T_mZZ2_4_2_Up = VBF_T_mZZ2_4_2_Up.Integral("width")
      integral_VBF_T_mZZ2_4_2_Down = VBF_T_mZZ2_4_2_Down.Integral("width")

      integral_VBF_T_mZZ2_4_2_OppositeDjet = VBF_T_mZZ2_4_2_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_4_2_Up_OppositeDjet = VBF_T_mZZ2_4_2_Up_OppositeDjet.Integral("width")
      integral_VBF_T_mZZ2_4_2_Down_OppositeDjet = VBF_T_mZZ2_4_2_Down_OppositeDjet.Integral("width")


      #-------- BACKGROUND TEMPLATES ------------------

      Bkg_T = sigTempFileU.Get("T_2D_qqZZ_UnConditional").Clone("Bkg_QCDContinuum")
      Bkg_ZX = sigTempFileU.Get("T_2D_ZX_UnConditional").Clone("Bkg_ZX_Nominal")
      Bkg_ZX_Up = sigTempFileUp_PDF.Get("T_2D_ZX_UnConditional").Clone("Bkg_ZX_Up")
      Bkg_ZX_Down = sigTempFileDown_PDF.Get("T_2D_ZX_UnConditional").Clone("Bkg_ZX_Down")

      #---------- RATES ---------------

      totalRate_ggzz = integral_Sig_T_1+integral_Sig_T_2+integral_Sig_T_4
      totalRate_ggzz_Shape = totalRate_ggzz * self.lumi
      rate_signal_ggzz_Shape = integral_Sig_T_2 * self.lumi  # *2.3
      rate_bkg_ggzz_Shape = integral_Sig_T_1 * self.lumi  # *2.3
      rate_interf_ggzz_Shape = integral_Sig_T_4 * self.lumi  # *2.3

      #----------- OWN-CATEGORY DJET ----------#

      sigRateQCDUpName = "signal_ggZZQCDUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateQCDUpName = "bkg_ggZZQCDUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDUpName = "interf_ggZZQCDUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp = ROOT.RooRealVar(sigRateQCDUpName, sigRateQCDUpName, integral_Sig_T_2_Up_QCD)
      bkgRates_QCDUp = ROOT.RooRealVar(bkgRateQCDUpName, bkgRateQCDUpName, integral_Sig_T_1_Up_QCD)
      interfRates_QCDUp = ROOT.RooRealVar(interfRateQCDUpName, interfRateQCDUpName, integral_Sig_T_4_Up_QCD)
      sigRates_QCDUp.setConstant(True)
      bkgRates_QCDUp.setConstant(True)
      interfRates_QCDUp.setConstant(True)

      sigRateQCDDownName = "signal_ggZZQCDDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateQCDDownName = "bkg_ggZZQCDDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDDownName = "interf_ggZZQCDDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown = ROOT.RooRealVar(sigRateQCDDownName, sigRateQCDDownName, integral_Sig_T_2_Down_QCD)
      bkgRates_QCDDown = ROOT.RooRealVar(bkgRateQCDDownName, bkgRateQCDDownName, integral_Sig_T_1_Down_QCD)
      interfRates_QCDDown = ROOT.RooRealVar(interfRateQCDDownName, interfRateQCDDownName, integral_Sig_T_4_Down_QCD)
      sigRates_QCDDown.setConstant(True)
      bkgRates_QCDDown.setConstant(True)
      interfRates_QCDDown.setConstant(True)

      sigRatePDFUpName = "signal_ggZZPDFUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRatePDFUpName = "bkg_ggZZPDFUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFUpName = "interf_ggZZPDFUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp = ROOT.RooRealVar(sigRatePDFUpName, sigRatePDFUpName, integral_Sig_T_2_Up_PDF)
      bkgRates_PDFUp = ROOT.RooRealVar(bkgRatePDFUpName, bkgRatePDFUpName, integral_Sig_T_1_Up_PDF)
      interfRates_PDFUp = ROOT.RooRealVar(interfRatePDFUpName, interfRatePDFUpName, integral_Sig_T_4_Up_PDF)
      sigRates_PDFUp.setConstant(True)
      bkgRates_PDFUp.setConstant(True)
      interfRates_PDFUp.setConstant(True)

      sigRatePDFDownName = "signal_ggZZPDFDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRatePDFDownName = "bkg_ggZZPDFDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFDownName = "interf_ggZZPDFDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown = ROOT.RooRealVar(sigRatePDFDownName, sigRatePDFDownName, integral_Sig_T_2_Down_PDF)
      bkgRates_PDFDown = ROOT.RooRealVar(bkgRatePDFDownName, bkgRatePDFDownName, integral_Sig_T_1_Down_PDF)
      interfRates_PDFDown = ROOT.RooRealVar(interfRatePDFDownName, interfRatePDFDownName, integral_Sig_T_4_Down_PDF)
      sigRates_PDFDown.setConstant(True)
      bkgRates_PDFDown.setConstant(True)
      interfRates_PDFDown.setConstant(True)

      sigRateNominalName = "signal_ggZZNominalrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateNominalName = "bkg_ggZZNominalrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateNominalName = "interf_ggZZNominalrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal = ROOT.RooRealVar(sigRateNominalName, sigRateNominalName, integral_Sig_T_2)
      bkgRates_Nominal = ROOT.RooRealVar(bkgRateNominalName, bkgRateNominalName, integral_Sig_T_1)
      interfRates_Nominal = ROOT.RooRealVar(interfRateNominalName, interfRateNominalName, integral_Sig_T_4)
      sigRates_Nominal.setConstant(True)
      bkgRates_Nominal.setConstant(True)
      interfRates_Nominal.setConstant(True)


      # ggH mZZ/mH**2

      sigRateQCDUp_mZZ2_1_Name = "signal_ggZZQCDUprate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDUp_mZZ2_1_Name = "interf_ggZZQCDUprate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_mZZ2_1 = ROOT.RooRealVar(sigRateQCDUp_mZZ2_1_Name, sigRateQCDUp_mZZ2_1_Name, integral_Sig_T_mZZ2_1_2_Up_QCD)
      interfRates_QCDUp_mZZ2_1 = ROOT.RooRealVar(interfRateQCDUp_mZZ2_1_Name, interfRateQCDUp_mZZ2_1_Name, integral_Sig_T_mZZ2_1_4_Up_QCD)
      sigRates_QCDUp_mZZ2_1.setConstant(True)
      interfRates_QCDUp_mZZ2_1.setConstant(True)

      sigRateQCDDown_mZZ2_1_Name = "signal_ggZZQCDDownrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDDown_mZZ2_1_Name = "interf_ggZZQCDDownrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_mZZ2_1 = ROOT.RooRealVar(sigRateQCDDown_mZZ2_1_Name, sigRateQCDDown_mZZ2_1_Name, integral_Sig_T_mZZ2_1_2_Down_QCD)
      interfRates_QCDDown_mZZ2_1 = ROOT.RooRealVar(interfRateQCDDown_mZZ2_1_Name, interfRateQCDDown_mZZ2_1_Name, integral_Sig_T_mZZ2_1_4_Down_QCD)
      sigRates_QCDDown_mZZ2_1.setConstant(True)
      interfRates_QCDDown_mZZ2_1.setConstant(True)

      sigRatePDFUp_mZZ2_1_Name = "signal_ggZZPDFUprate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFUp_mZZ2_1_Name = "interf_ggZZPDFUprate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_mZZ2_1 = ROOT.RooRealVar(sigRatePDFUp_mZZ2_1_Name, sigRatePDFUp_mZZ2_1_Name, integral_Sig_T_mZZ2_1_2_Up_PDF)
      interfRates_PDFUp_mZZ2_1 = ROOT.RooRealVar(interfRatePDFUp_mZZ2_1_Name, interfRatePDFUp_mZZ2_1_Name, integral_Sig_T_mZZ2_1_4_Up_PDF)
      sigRates_PDFUp_mZZ2_1.setConstant(True)
      interfRates_PDFUp_mZZ2_1.setConstant(True)

      sigRatePDFDown_mZZ2_1_Name = "signal_ggZZPDFDownrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFDown_mZZ2_1_Name = "interf_ggZZPDFDownrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_mZZ2_1 = ROOT.RooRealVar(sigRatePDFDown_mZZ2_1_Name, sigRatePDFDown_mZZ2_1_Name, integral_Sig_T_mZZ2_1_2_Down_PDF)
      interfRates_PDFDown_mZZ2_1 = ROOT.RooRealVar(interfRatePDFDown_mZZ2_1_Name, interfRatePDFDown_mZZ2_1_Name, integral_Sig_T_mZZ2_1_4_Down_PDF)
      sigRates_PDFDown_mZZ2_1.setConstant(True)
      interfRates_PDFDown_mZZ2_1.setConstant(True)

      sigRateNominal_mZZ2_1_Name = "signal_ggZZNominalrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateNominal_mZZ2_1_Name = "interf_ggZZNominalrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_mZZ2_1 = ROOT.RooRealVar(sigRateNominal_mZZ2_1_Name, sigRateNominal_mZZ2_1_Name, integral_Sig_T_mZZ2_1_2)
      interfRates_Nominal_mZZ2_1 = ROOT.RooRealVar(interfRateNominal_mZZ2_1_Name, interfRateNominal_mZZ2_1_Name, integral_Sig_T_mZZ2_1_4)
      sigRates_Nominal_mZZ2_1.setConstant(True)
      interfRates_Nominal_mZZ2_1.setConstant(True)


      # ggH mZZ/mH**2**2

      sigRateQCDUp_mZZ2_2_Name = "signal_ggZZQCDUprate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_mZZ2_2 = ROOT.RooRealVar(sigRateQCDUp_mZZ2_2_Name, sigRateQCDUp_mZZ2_2_Name, integral_Sig_T_mZZ2_2_2_Up_QCD)
      sigRates_QCDUp_mZZ2_2.setConstant(True)

      sigRateQCDDown_mZZ2_2_Name = "signal_ggZZQCDDownrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_mZZ2_2 = ROOT.RooRealVar(sigRateQCDDown_mZZ2_2_Name, sigRateQCDDown_mZZ2_2_Name, integral_Sig_T_mZZ2_2_2_Down_QCD)
      sigRates_QCDDown_mZZ2_2.setConstant(True)

      sigRatePDFUp_mZZ2_2_Name = "signal_ggZZPDFUprate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_mZZ2_2 = ROOT.RooRealVar(sigRatePDFUp_mZZ2_2_Name, sigRatePDFUp_mZZ2_2_Name, integral_Sig_T_mZZ2_2_2_Up_PDF)
      sigRates_PDFUp_mZZ2_2.setConstant(True)

      sigRatePDFDown_mZZ2_2_Name = "signal_ggZZPDFDownrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_mZZ2_2 = ROOT.RooRealVar(sigRatePDFDown_mZZ2_2_Name, sigRatePDFDown_mZZ2_2_Name, integral_Sig_T_mZZ2_2_2_Down_PDF)
      sigRates_PDFDown_mZZ2_2.setConstant(True)

      sigRateNominal_mZZ2_2_Name = "signal_ggZZNominalrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_mZZ2_2 = ROOT.RooRealVar(sigRateNominal_mZZ2_2_Name, sigRateNominal_mZZ2_2_Name, integral_Sig_T_mZZ2_2_2)
      sigRates_Nominal_mZZ2_2.setConstant(True)


      #------------ OPPOSITE DJET CATEGORY --------------#


      sigRateQCDUpName_OppositeDjet = "signal_ggZZQCDUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateQCDUpName_OppositeDjet = "bkg_ggZZQCDUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDUpName_OppositeDjet = "interf_ggZZQCDUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_OppositeDjet = ROOT.RooRealVar(sigRateQCDUpName_OppositeDjet, sigRateQCDUpName_OppositeDjet, integral_Sig_T_2_Up_QCD_OppositeDjet)
      bkgRates_QCDUp_OppositeDjet = ROOT.RooRealVar(bkgRateQCDUpName_OppositeDjet, bkgRateQCDUpName_OppositeDjet, integral_Sig_T_1_Up_QCD_OppositeDjet)
      interfRates_QCDUp_OppositeDjet = ROOT.RooRealVar(interfRateQCDUpName_OppositeDjet, interfRateQCDUpName_OppositeDjet, integral_Sig_T_4_Up_QCD_OppositeDjet)
      sigRates_QCDUp_OppositeDjet.setConstant(True)
      bkgRates_QCDUp_OppositeDjet.setConstant(True)
      interfRates_QCDUp_OppositeDjet.setConstant(True)

      sigRateQCDDownName_OppositeDjet = "signal_ggZZQCDDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateQCDDownName_OppositeDjet = "bkg_ggZZQCDDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDDownName_OppositeDjet = "interf_ggZZQCDDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_OppositeDjet = ROOT.RooRealVar(sigRateQCDDownName_OppositeDjet, sigRateQCDDownName_OppositeDjet, integral_Sig_T_2_Down_QCD_OppositeDjet)
      bkgRates_QCDDown_OppositeDjet = ROOT.RooRealVar(bkgRateQCDDownName_OppositeDjet, bkgRateQCDDownName_OppositeDjet, integral_Sig_T_1_Down_QCD_OppositeDjet)
      interfRates_QCDDown_OppositeDjet = ROOT.RooRealVar(interfRateQCDDownName_OppositeDjet, interfRateQCDDownName_OppositeDjet, integral_Sig_T_4_Down_QCD_OppositeDjet)
      sigRates_QCDDown_OppositeDjet.setConstant(True)
      bkgRates_QCDDown_OppositeDjet.setConstant(True)
      interfRates_QCDDown_OppositeDjet.setConstant(True)

      sigRatePDFUpName_OppositeDjet = "signal_ggZZPDFUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRatePDFUpName_OppositeDjet = "bkg_ggZZPDFUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFUpName_OppositeDjet = "interf_ggZZPDFUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_OppositeDjet = ROOT.RooRealVar(sigRatePDFUpName_OppositeDjet, sigRatePDFUpName_OppositeDjet, integral_Sig_T_2_Up_PDF_OppositeDjet)
      bkgRates_PDFUp_OppositeDjet = ROOT.RooRealVar(bkgRatePDFUpName_OppositeDjet, bkgRatePDFUpName_OppositeDjet, integral_Sig_T_1_Up_PDF_OppositeDjet)
      interfRates_PDFUp_OppositeDjet = ROOT.RooRealVar(interfRatePDFUpName_OppositeDjet, interfRatePDFUpName_OppositeDjet, integral_Sig_T_4_Up_PDF_OppositeDjet)
      sigRates_PDFUp_OppositeDjet.setConstant(True)
      bkgRates_PDFUp_OppositeDjet.setConstant(True)
      interfRates_PDFUp_OppositeDjet.setConstant(True)

      sigRatePDFDownName_OppositeDjet = "signal_ggZZPDFDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRatePDFDownName_OppositeDjet = "bkg_ggZZPDFDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFDownName_OppositeDjet = "interf_ggZZPDFDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_OppositeDjet = ROOT.RooRealVar(sigRatePDFDownName_OppositeDjet, sigRatePDFDownName_OppositeDjet, integral_Sig_T_2_Down_PDF_OppositeDjet)
      bkgRates_PDFDown_OppositeDjet = ROOT.RooRealVar(bkgRatePDFDownName_OppositeDjet, bkgRatePDFDownName_OppositeDjet, integral_Sig_T_1_Down_PDF_OppositeDjet)
      interfRates_PDFDown_OppositeDjet = ROOT.RooRealVar(interfRatePDFDownName_OppositeDjet, interfRatePDFDownName_OppositeDjet, integral_Sig_T_4_Down_PDF_OppositeDjet)
      sigRates_PDFDown_OppositeDjet.setConstant(True)
      bkgRates_PDFDown_OppositeDjet.setConstant(True)
      interfRates_PDFDown_OppositeDjet.setConstant(True)

      sigRateNominalName_OppositeDjet = "signal_ggZZNominal_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateNominalName_OppositeDjet = "bkg_ggZZNominal_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateNominalName_OppositeDjet = "interf_ggZZNominal_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_OppositeDjet = ROOT.RooRealVar(sigRateNominalName_OppositeDjet, sigRateNominalName_OppositeDjet, integral_Sig_T_2_OppositeDjet)
      bkgRates_Nominal_OppositeDjet = ROOT.RooRealVar(bkgRateNominalName_OppositeDjet, bkgRateNominalName_OppositeDjet, integral_Sig_T_1_OppositeDjet)
      interfRates_Nominal_OppositeDjet = ROOT.RooRealVar(interfRateNominalName_OppositeDjet, interfRateNominalName_OppositeDjet, integral_Sig_T_4_OppositeDjet)
      sigRates_Nominal_OppositeDjet.setConstant(True)
      bkgRates_Nominal_OppositeDjet.setConstant(True)
      interfRates_Nominal_OppositeDjet.setConstant(True)


      # ggH mZZ/mH**2

      sigRateQCDUp_mZZ2_1_Name_OppositeDjet = "signal_ggZZQCDUp_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDUp_mZZ2_1_Name_OppositeDjet = "interf_ggZZQCDUp_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_mZZ2_1_OppositeDjet = ROOT.RooRealVar(sigRateQCDUp_mZZ2_1_Name_OppositeDjet, sigRateQCDUp_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_2_Up_QCD_OppositeDjet)
      interfRates_QCDUp_mZZ2_1_OppositeDjet = ROOT.RooRealVar(interfRateQCDUp_mZZ2_1_Name_OppositeDjet, interfRateQCDUp_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_4_Up_QCD_OppositeDjet)
      sigRates_QCDUp_mZZ2_1_OppositeDjet.setConstant(True)
      interfRates_QCDUp_mZZ2_1_OppositeDjet.setConstant(True)

      sigRateQCDDown_mZZ2_1_Name_OppositeDjet = "signal_ggZZQCDDown_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateQCDDown_mZZ2_1_Name_OppositeDjet = "interf_ggZZQCDDown_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_mZZ2_1_OppositeDjet = ROOT.RooRealVar(sigRateQCDDown_mZZ2_1_Name_OppositeDjet, sigRateQCDDown_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_2_Down_QCD_OppositeDjet)
      interfRates_QCDDown_mZZ2_1_OppositeDjet = ROOT.RooRealVar(interfRateQCDDown_mZZ2_1_Name_OppositeDjet, interfRateQCDDown_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_4_Down_QCD_OppositeDjet)
      sigRates_QCDDown_mZZ2_1_OppositeDjet.setConstant(True)
      interfRates_QCDDown_mZZ2_1_OppositeDjet.setConstant(True)

      sigRatePDFUp_mZZ2_1_Name_OppositeDjet = "signal_ggZZPDFUp_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFUp_mZZ2_1_Name_OppositeDjet = "interf_ggZZPDFUp_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_mZZ2_1_OppositeDjet = ROOT.RooRealVar(sigRatePDFUp_mZZ2_1_Name_OppositeDjet, sigRatePDFUp_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_2_Up_PDF_OppositeDjet)
      interfRates_PDFUp_mZZ2_1_OppositeDjet = ROOT.RooRealVar(interfRatePDFUp_mZZ2_1_Name_OppositeDjet, interfRatePDFUp_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_4_Up_PDF_OppositeDjet)
      sigRates_PDFUp_mZZ2_1_OppositeDjet.setConstant(True)
      interfRates_PDFUp_mZZ2_1_OppositeDjet.setConstant(True)

      sigRatePDFDown_mZZ2_1_Name_OppositeDjet = "signal_ggZZPDFDown_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRatePDFDown_mZZ2_1_Name_OppositeDjet = "interf_ggZZPDFDown_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_mZZ2_1_OppositeDjet = ROOT.RooRealVar(sigRatePDFDown_mZZ2_1_Name_OppositeDjet, sigRatePDFDown_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_2_Down_PDF_OppositeDjet)
      interfRates_PDFDown_mZZ2_1_OppositeDjet = ROOT.RooRealVar(interfRatePDFDown_mZZ2_1_Name_OppositeDjet, interfRatePDFDown_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_4_Down_PDF_OppositeDjet)
      sigRates_PDFDown_mZZ2_1_OppositeDjet.setConstant(True)
      interfRates_PDFDown_mZZ2_1_OppositeDjet.setConstant(True)

      sigRateNominal_mZZ2_1_Name_OppositeDjet = "signal_ggZZNominal_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRateNominal_mZZ2_1_Name_OppositeDjet = "interf_ggZZNominal_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_mZZ2_1_OppositeDjet = ROOT.RooRealVar(sigRateNominal_mZZ2_1_Name_OppositeDjet, sigRateNominal_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_2_OppositeDjet)
      interfRates_Nominal_mZZ2_1_OppositeDjet = ROOT.RooRealVar(interfRateNominal_mZZ2_1_Name_OppositeDjet, interfRateNominal_mZZ2_1_Name_OppositeDjet, integral_Sig_T_mZZ2_1_4_OppositeDjet)
      sigRates_Nominal_mZZ2_1_OppositeDjet.setConstant(True)
      interfRates_Nominal_mZZ2_1_OppositeDjet.setConstant(True)


      # ggH mZZ/mH**2**2

      sigRateQCDUp_mZZ2_2_Name_OppositeDjet = "signal_ggZZQCDUp_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_mZZ2_2_OppositeDjet = ROOT.RooRealVar(sigRateQCDUp_mZZ2_2_Name_OppositeDjet, sigRateQCDUp_mZZ2_2_Name_OppositeDjet, integral_Sig_T_mZZ2_2_2_Up_QCD_OppositeDjet)
      sigRates_QCDUp_mZZ2_2_OppositeDjet.setConstant(True)

      sigRateQCDDown_mZZ2_2_Name_OppositeDjet = "signal_ggZZQCDDown_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_mZZ2_2_OppositeDjet = ROOT.RooRealVar(sigRateQCDDown_mZZ2_2_Name_OppositeDjet, sigRateQCDDown_mZZ2_2_Name_OppositeDjet, integral_Sig_T_mZZ2_2_2_Down_QCD_OppositeDjet)
      sigRates_QCDDown_mZZ2_2_OppositeDjet.setConstant(True)

      sigRatePDFUp_mZZ2_2_Name_OppositeDjet = "signal_ggZZPDFUp_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_mZZ2_2_OppositeDjet = ROOT.RooRealVar(sigRatePDFUp_mZZ2_2_Name_OppositeDjet, sigRatePDFUp_mZZ2_2_Name_OppositeDjet, integral_Sig_T_mZZ2_2_2_Up_PDF_OppositeDjet)
      sigRates_PDFUp_mZZ2_2_OppositeDjet.setConstant(True)

      sigRatePDFDown_mZZ2_2_Name_OppositeDjet = "signal_ggZZPDFDown_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_mZZ2_2_OppositeDjet = ROOT.RooRealVar(sigRatePDFDown_mZZ2_2_Name_OppositeDjet, sigRatePDFDown_mZZ2_2_Name_OppositeDjet, integral_Sig_T_mZZ2_2_2_Down_PDF_OppositeDjet)
      sigRates_PDFDown_mZZ2_2_OppositeDjet.setConstant(True)

      sigRateNominal_mZZ2_2_Name_OppositeDjet = "signal_ggZZNominal_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_mZZ2_2_OppositeDjet = ROOT.RooRealVar(sigRateNominal_mZZ2_2_Name_OppositeDjet, sigRateNominal_mZZ2_2_Name_OppositeDjet, integral_Sig_T_mZZ2_2_2_OppositeDjet)
      sigRates_Nominal_mZZ2_2_OppositeDjet.setConstant(True)

      #---------------------------------------------------#


      # ggH anomalous coupling rate parameterizations

      sigRates_Nominal_AnomCoupl_Name = "signal_ggZZNominalrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_AnomCoupl = ROOT.RooFormulaVar(
      sigRates_Nominal_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_Nominal, sigRates_Nominal_mZZ2_1, sigRates_Nominal_mZZ2_2, fai1)
      )
      interfRates_Nominal_AnomCoupl_Name = "interf_ggZZNominalrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_Nominal_AnomCoupl = ROOT.RooFormulaVar(
      interfRates_Nominal_AnomCoupl_Name,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_Nominal, interfRates_Nominal_mZZ2_1, fai1)
      )

      sigRates_QCDUp_AnomCoupl_Name = "signal_ggZZQCDUprate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_AnomCoupl = ROOT.RooFormulaVar(
      sigRates_QCDUp_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_QCDUp, sigRates_QCDUp_mZZ2_1, sigRates_QCDUp_mZZ2_2, fai1)
      )
      interfRates_QCDUp_AnomCoupl_Name = "interf_ggZZQCDUprate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_QCDUp_AnomCoupl = ROOT.RooFormulaVar(
      interfRates_QCDUp_AnomCoupl_Name,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_QCDUp, interfRates_QCDUp_mZZ2_1, fai1)
      )

      sigRates_QCDDown_AnomCoupl_Name = "signal_ggZZQCDDownrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_AnomCoupl = ROOT.RooFormulaVar(
      sigRates_QCDDown_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_QCDDown, sigRates_QCDDown_mZZ2_1, sigRates_QCDDown_mZZ2_2, fai1)
      )
      interfRates_QCDDown_AnomCoupl_Name = "interf_ggZZQCDDownrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_QCDDown_AnomCoupl = ROOT.RooFormulaVar(
      interfRates_QCDDown_AnomCoupl_Name,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_QCDDown, interfRates_QCDDown_mZZ2_1, fai1)
      )

      sigRates_PDFUp_AnomCoupl_Name = "signal_ggZZPDFUprate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_AnomCoupl = ROOT.RooFormulaVar(
      sigRates_PDFUp_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_PDFUp, sigRates_PDFUp_mZZ2_1, sigRates_PDFUp_mZZ2_2, fai1)
      )
      interfRates_PDFUp_AnomCoupl_Name = "interf_ggZZPDFUprate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_PDFUp_AnomCoupl = ROOT.RooFormulaVar(
      interfRates_PDFUp_AnomCoupl_Name,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_PDFUp, interfRates_PDFUp_mZZ2_1, fai1)
      )

      sigRates_PDFDown_AnomCoupl_Name = "signal_ggZZPDFDownrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_AnomCoupl = ROOT.RooFormulaVar(
      sigRates_PDFDown_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_PDFDown, sigRates_PDFDown_mZZ2_1, sigRates_PDFDown_mZZ2_2, fai1)
      )
      interfRates_PDFDown_AnomCoupl_Name = "interf_ggZZPDFDownrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_PDFDown_AnomCoupl = ROOT.RooFormulaVar(
      interfRates_PDFDown_AnomCoupl_Name,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_PDFDown, interfRates_PDFDown_mZZ2_1, fai1)
      )



      sigRates_Nominal_AnomCoupl_Name_OppositeDjet = "signal_ggZZNominalrate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      sigRates_Nominal_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_Nominal_OppositeDjet, sigRates_Nominal_mZZ2_1_OppositeDjet, sigRates_Nominal_mZZ2_2_OppositeDjet, fai1)
      )
      interfRates_Nominal_AnomCoupl_Name_OppositeDjet = "interf_ggZZNominalrate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_Nominal_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      interfRates_Nominal_AnomCoupl_Name_OppositeDjet,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_Nominal_OppositeDjet, interfRates_Nominal_mZZ2_1_OppositeDjet, fai1)
      )

      sigRates_QCDUp_AnomCoupl_Name_OppositeDjet = "signal_ggZZQCDUprate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      sigRates_QCDUp_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_QCDUp_OppositeDjet, sigRates_QCDUp_mZZ2_1_OppositeDjet, sigRates_QCDUp_mZZ2_2_OppositeDjet, fai1)
      )
      interfRates_QCDUp_AnomCoupl_Name_OppositeDjet = "interf_ggZZQCDUprate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_QCDUp_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      interfRates_QCDUp_AnomCoupl_Name_OppositeDjet,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_QCDUp_OppositeDjet, interfRates_QCDUp_mZZ2_1_OppositeDjet, fai1)
      )

      sigRates_QCDDown_AnomCoupl_Name_OppositeDjet = "signal_ggZZQCDDownrate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      sigRates_QCDDown_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_QCDDown_OppositeDjet, sigRates_QCDDown_mZZ2_1_OppositeDjet, sigRates_QCDDown_mZZ2_2_OppositeDjet, fai1)
      )
      interfRates_QCDDown_AnomCoupl_Name_OppositeDjet = "interf_ggZZQCDDownrate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_QCDDown_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      interfRates_QCDDown_AnomCoupl_Name_OppositeDjet,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_QCDDown_OppositeDjet, interfRates_QCDDown_mZZ2_1_OppositeDjet, fai1)
      )

      sigRates_PDFUp_AnomCoupl_Name_OppositeDjet = "signal_ggZZPDFUprate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      sigRates_PDFUp_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_PDFUp_OppositeDjet, sigRates_PDFUp_mZZ2_1_OppositeDjet, sigRates_PDFUp_mZZ2_2_OppositeDjet, fai1)
      )
      interfRates_PDFUp_AnomCoupl_Name_OppositeDjet = "interf_ggZZPDFUprate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_PDFUp_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      interfRates_PDFUp_AnomCoupl_Name_OppositeDjet,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_PDFUp_OppositeDjet, interfRates_PDFUp_mZZ2_1_OppositeDjet, fai1)
      )

      sigRates_PDFDown_AnomCoupl_Name_OppositeDjet = "signal_ggZZPDFDownrate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      sigRates_PDFDown_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(sigRates_PDFDown_OppositeDjet, sigRates_PDFDown_mZZ2_1_OppositeDjet, sigRates_PDFDown_mZZ2_2_OppositeDjet, fai1)
      )
      interfRates_PDFDown_AnomCoupl_Name_OppositeDjet = "interf_ggZZPDFDownrate_AnomCoupl_OppositeDjet_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_PDFDown_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      interfRates_PDFDown_AnomCoupl_Name_OppositeDjet,
      "( sqrt(1-abs(@2))*@0 + sign(@2)*sqrt(abs(@2))*@1 )",
      ROOT.RooArgList(interfRates_PDFDown_OppositeDjet, interfRates_PDFDown_mZZ2_1_OppositeDjet, fai1)
      )


      #------- Combined ggZZ Djet ----------------#

      sigRates_Nominal_AnomCoupl_CombinedJet_Name = "signal_ggZZNominal_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_Nominal_AnomCoupl_CombinedJet_Name = "interf_ggZZNominal_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateNominalName_CombinedDjet = "bkg_ggZZNominal_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(sigRates_Nominal_AnomCoupl)
      )
      interfRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(interfRates_Nominal_AnomCoupl)
      )
      bkgRates_Nominal_CombinedJet = ROOT.RooFormulaVar(
      bkgRateNominalName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(bkgRates_Nominal)
      )
      if useDjet == 1:
      sigRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(sigRates_Nominal_AnomCoupl, sigRates_Nominal_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(interfRates_Nominal_AnomCoupl, interfRates_Nominal_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_Nominal_CombinedJet = ROOT.RooFormulaVar(
      bkgRateNominalName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(bkgRates_Nominal, bkgRates_Nominal_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      if useDjet == 2:
      sigRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(sigRates_Nominal_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(interfRates_Nominal_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_Nominal_CombinedJet = ROOT.RooFormulaVar(
      bkgRateNominalName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(bkgRates_Nominal, thetaSyst_djet_ggZZ_norm)
      )

      sigRates_QCDUp_AnomCoupl_CombinedJet_Name = "signal_ggZZQCDUp_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_QCDUp_AnomCoupl_CombinedJet_Name = "interf_ggZZQCDUp_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateQCDUpName_CombinedDjet = "bkg_ggZZQCDUp_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_QCDUp_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(sigRates_QCDUp_AnomCoupl)
      )
      interfRates_QCDUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_QCDUp_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(interfRates_QCDUp_AnomCoupl)
      )
      bkgRates_QCDUp_CombinedJet = ROOT.RooFormulaVar(
      bkgRateQCDUpName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(bkgRates_QCDUp)
      )
      if useDjet == 1:
      sigRates_QCDUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_QCDUp_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(sigRates_QCDUp_AnomCoupl, sigRates_QCDUp_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_QCDUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_QCDUp_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(interfRates_QCDUp_AnomCoupl, interfRates_QCDUp_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_QCDUp_CombinedJet = ROOT.RooFormulaVar(
      bkgRateQCDUpName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(bkgRates_QCDUp, bkgRates_QCDUp_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      if useDjet == 2:
      sigRates_QCDUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_QCDUp_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(sigRates_QCDUp_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_QCDUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_QCDUp_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(interfRates_QCDUp_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_QCDUp_CombinedJet = ROOT.RooFormulaVar(
      bkgRateQCDUpName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(bkgRates_QCDUp, thetaSyst_djet_ggZZ_norm)
      )

      sigRates_QCDDown_AnomCoupl_CombinedJet_Name = "signal_ggZZQCDDown_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_QCDDown_AnomCoupl_CombinedJet_Name = "interf_ggZZQCDDown_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRateQCDDownName_CombinedDjet = "bkg_ggZZQCDDown_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_QCDDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_QCDDown_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(sigRates_QCDDown_AnomCoupl)
      )
      interfRates_QCDDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_QCDDown_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(interfRates_QCDDown_AnomCoupl)
      )
      bkgRates_QCDDown_CombinedJet = ROOT.RooFormulaVar(
      bkgRateQCDDownName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(bkgRates_QCDDown)
      )
      if useDjet == 1:
      sigRates_QCDDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_QCDDown_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(sigRates_QCDDown_AnomCoupl, sigRates_QCDDown_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_QCDDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_QCDDown_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(interfRates_QCDDown_AnomCoupl, interfRates_QCDDown_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_QCDDown_CombinedJet = ROOT.RooFormulaVar(
      bkgRateQCDDownName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(bkgRates_QCDDown, bkgRates_QCDDown_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      if useDjet == 2:
      sigRates_QCDDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_QCDDown_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(sigRates_QCDDown_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_QCDDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_QCDDown_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(interfRates_QCDDown_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_QCDDown_CombinedJet = ROOT.RooFormulaVar(
      bkgRateQCDDownName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(bkgRates_QCDDown, thetaSyst_djet_ggZZ_norm)
      )


      sigRates_PDFUp_AnomCoupl_CombinedJet_Name = "signal_ggZZPDFUp_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_PDFUp_AnomCoupl_CombinedJet_Name = "interf_ggZZPDFUp_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRatePDFUpName_CombinedDjet = "bkg_ggZZPDFUp_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_PDFUp_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(sigRates_PDFUp_AnomCoupl)
      )
      interfRates_PDFUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_PDFUp_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(interfRates_PDFUp_AnomCoupl)
      )
      bkgRates_PDFUp_CombinedJet = ROOT.RooFormulaVar(
      bkgRatePDFUpName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(bkgRates_PDFUp)
      )
      if useDjet == 1:
      sigRates_PDFUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_PDFUp_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(sigRates_PDFUp_AnomCoupl, sigRates_PDFUp_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_PDFUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_PDFUp_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(interfRates_PDFUp_AnomCoupl, interfRates_PDFUp_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_PDFUp_CombinedJet = ROOT.RooFormulaVar(
      bkgRatePDFUpName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(bkgRates_PDFUp, bkgRates_PDFUp_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      if useDjet == 2:
      sigRates_PDFUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_PDFUp_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(sigRates_PDFUp_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_PDFUp_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_PDFUp_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(interfRates_PDFUp_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_PDFUp_CombinedJet = ROOT.RooFormulaVar(
      bkgRatePDFUpName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(bkgRates_PDFUp, thetaSyst_djet_ggZZ_norm)
      )

      sigRates_PDFDown_AnomCoupl_CombinedJet_Name = "signal_ggZZPDFDown_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      interfRates_PDFDown_AnomCoupl_CombinedJet_Name = "interf_ggZZPDFDown_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkgRatePDFDownName_CombinedDjet = "bkg_ggZZPDFDown_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      sigRates_PDFDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_PDFDown_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(sigRates_PDFDown_AnomCoupl)
      )
      interfRates_PDFDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_PDFDown_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(interfRates_PDFDown_AnomCoupl)
      )
      bkgRates_PDFDown_CombinedJet = ROOT.RooFormulaVar(
      bkgRatePDFDownName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(bkgRates_PDFDown)
      )
      if useDjet == 1:
      sigRates_PDFDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_PDFDown_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(sigRates_PDFDown_AnomCoupl, sigRates_PDFDown_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_PDFDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_PDFDown_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(interfRates_PDFDown_AnomCoupl, interfRates_PDFDown_AnomCoupl_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_PDFDown_CombinedJet = ROOT.RooFormulaVar(
      bkgRatePDFDownName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(bkgRates_PDFDown, bkgRates_PDFDown_OppositeDjet, thetaSyst_djet_ggZZ_norm)
      )
      if useDjet == 2:
      sigRates_PDFDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      sigRates_PDFDown_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(sigRates_PDFDown_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      interfRates_PDFDown_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      interfRates_PDFDown_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(interfRates_PDFDown_AnomCoupl, thetaSyst_djet_ggZZ_norm)
      )
      bkgRates_PDFDown_CombinedJet = ROOT.RooFormulaVar(
      bkgRatePDFDownName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(bkgRates_PDFDown, thetaSyst_djet_ggZZ_norm)
      )



      ggZZVarNorm_Name = "ggZZVarNominalNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZNominal_norm = ROOT.RooFormulaVar(
      ggZZVarNorm_Name, "(@0*@3*@6*@4+@1*sqrt(@3*@6*@4)*sign(@5)*sqrt(abs(@5))+@2*abs(@5))",
      ROOT.RooArgList(sigRates_Nominal_AnomCoupl_CombinedJet, interfRates_Nominal_AnomCoupl_CombinedJet, bkgRates_Nominal_CombinedJet, x, mu, kbkg_gg, muF)
      )
      ggZZVarNormQCDUp_Name = "ggZZVarQCDUpNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZQCDUp_norm = ROOT.RooFormulaVar(
      ggZZVarNormQCDUp_Name, "(@0*@3*@6*@4+@1*sqrt(@3*@6*@4)*sign(@5)*sqrt(abs(@5))+@2*abs(@5))",
      ROOT.RooArgList(sigRates_QCDUp_AnomCoupl_CombinedJet, interfRates_QCDUp_AnomCoupl_CombinedJet, bkgRates_QCDUp_CombinedJet, x, mu, kbkg_gg, muF)
      )
      ggZZVarNormQCDDown_Name = "ggZZVarQCDDownNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZQCDDown_norm = ROOT.RooFormulaVar(
      ggZZVarNormQCDDown_Name, "(@0*@3*@6*@4+@1*sqrt(@3*@6*@4)*sign(@5)*sqrt(abs(@5))+@2*abs(@5))",
      ROOT.RooArgList(sigRates_QCDDown_AnomCoupl_CombinedJet, interfRates_QCDDown_AnomCoupl_CombinedJet, bkgRates_QCDDown_CombinedJet, x, mu, kbkg_gg, muF)
      )
      ggZZVarNormPDFUp_Name = "ggZZVarPDFUpNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZPDFUp_norm = ROOT.RooFormulaVar(
      ggZZVarNormPDFUp_Name, "(@0*@3*@6*@4+@1*sqrt(@3*@6*@4)*sign(@5)*sqrt(abs(@5))+@2*abs(@5))",
      ROOT.RooArgList(sigRates_PDFUp_AnomCoupl_CombinedJet, interfRates_PDFUp_AnomCoupl_CombinedJet, bkgRates_PDFUp_CombinedJet, x, mu, kbkg_gg, muF)
      )
      ggZZVarNormPDFDown_Name = "ggZZVarPDFDownNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZPDFDown_norm = ROOT.RooFormulaVar(
      ggZZVarNormPDFDown_Name, "(@0*@3*@6*@4+@1*sqrt(@3*@6*@4)*sign(@5)*sqrt(abs(@5))+@2*abs(@5))",
      ROOT.RooArgList(sigRates_PDFDown_AnomCoupl_CombinedJet, interfRates_PDFDown_AnomCoupl_CombinedJet, bkgRates_PDFDown_CombinedJet, x, mu, kbkg_gg, muF)
      )


      #--------------- VBF ------------

      # Set names, bounds, and values for rates of VBF
      totalRate_vbf = integral_VBF_T_1+integral_VBF_T_2+integral_VBF_T_4
      totalRate_vbf_Shape = totalRate_vbf * self.lumi
      rate_signal_vbf_Shape = integral_VBF_T_2 * self.lumi
      rate_bkg_vbf_Shape = integral_VBF_T_1 * self.lumi
      rate_interf_vbf_Shape = integral_VBF_T_4 * self.lumi

      VBFsigRateUpName = "signal_VBFUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateUpName = "bkg_VBFUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateUpName = "interf_VBFUprate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Up = ROOT.RooRealVar(VBFsigRateUpName, VBFsigRateUpName, integral_VBF_T_2_Up)
      VBFbkgRates_Up = ROOT.RooRealVar(VBFbkgRateUpName, VBFbkgRateUpName, integral_VBF_T_1_Up)
      VBFinterfRates_Up = ROOT.RooRealVar(VBFinterfRateUpName, VBFinterfRateUpName, integral_VBF_T_4_Up)
      VBFsigRates_Up.setConstant(True)
      VBFbkgRates_Up.setConstant(True)
      VBFinterfRates_Up.setConstant(True)

      VBFsigRateDownName = "signal_VBFDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateDownName = "bkg_VBFDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateDownName = "interf_VBFDownrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Down = ROOT.RooRealVar(VBFsigRateDownName, VBFsigRateDownName, integral_VBF_T_2_Down)
      VBFbkgRates_Down = ROOT.RooRealVar(VBFbkgRateDownName, VBFbkgRateDownName, integral_VBF_T_1_Down)
      VBFinterfRates_Down = ROOT.RooRealVar(VBFinterfRateDownName, VBFinterfRateDownName, integral_VBF_T_4_Down)
      VBFsigRates_Down.setConstant(True)
      VBFbkgRates_Down.setConstant(True)
      VBFinterfRates_Down.setConstant(True)

      VBFsigRateNominalName = "signal_VBFNominalrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateNominalName = "bkg_VBFNominalrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateNominalName = "interf_VBFNominalrate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Nominal = ROOT.RooRealVar(VBFsigRateNominalName, VBFsigRateNominalName, integral_VBF_T_2)
      VBFbkgRates_Nominal = ROOT.RooRealVar(VBFbkgRateNominalName, VBFbkgRateNominalName, integral_VBF_T_1)
      VBFinterfRates_Nominal = ROOT.RooRealVar(VBFinterfRateNominalName, VBFinterfRateNominalName, integral_VBF_T_4)
      VBFsigRates_Nominal.setConstant(True)
      VBFbkgRates_Nominal.setConstant(True)
      VBFinterfRates_Nominal.setConstant(True)


      # VBF mZZ/mH**2

      VBFsigRateUp_mZZ2_1_Name = "signal_VBFUprate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateUp_mZZ2_1_Name = "interf_VBFUprate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_1_Up = ROOT.RooRealVar(VBFsigRateUp_mZZ2_1_Name, VBFsigRateUp_mZZ2_1_Name, integral_VBF_T_mZZ2_1_2_Up)
      VBFinterfRates_mZZ2_1_Up = ROOT.RooRealVar(VBFinterfRateUp_mZZ2_1_Name, VBFinterfRateUp_mZZ2_1_Name, integral_VBF_T_mZZ2_1_4_Up)
      VBFsigRates_mZZ2_1_Up.setConstant(True)
      VBFinterfRates_mZZ2_1_Up.setConstant(True)

      VBFsigRateDown_mZZ2_1_Name = "signal_VBFDownrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateDown_mZZ2_1_Name = "interf_VBFDownrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_1_Down = ROOT.RooRealVar(VBFsigRateDown_mZZ2_1_Name, VBFsigRateDown_mZZ2_1_Name, integral_VBF_T_mZZ2_1_2_Down)
      VBFinterfRates_mZZ2_1_Down = ROOT.RooRealVar(VBFinterfRateDown_mZZ2_1_Name, VBFinterfRateDown_mZZ2_1_Name, integral_VBF_T_mZZ2_1_4_Down)
      VBFsigRates_mZZ2_1_Down.setConstant(True)
      VBFinterfRates_mZZ2_1_Down.setConstant(True)

      VBFsigRateNominal_mZZ2_1_Name = "signal_VBFNominalrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateNominal_mZZ2_1_Name = "interf_VBFNominalrate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_1_Nominal = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_1_Name, VBFsigRateNominal_mZZ2_1_Name, integral_VBF_T_mZZ2_1_2)
      VBFinterfRates_mZZ2_1_Nominal = ROOT.RooRealVar(VBFinterfRateNominal_mZZ2_1_Name, VBFinterfRateNominal_mZZ2_1_Name, integral_VBF_T_mZZ2_1_4)
      VBFsigRates_mZZ2_1_Nominal.setConstant(True)
      VBFinterfRates_mZZ2_1_Nominal.setConstant(True)

      # VBF mZZ/mH**2**2

      VBFsigRateUp_mZZ2_2_Name = "signal_VBFUprate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateUp_mZZ2_2_Name = "interf_VBFUprate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_2_Up = ROOT.RooRealVar(VBFsigRateUp_mZZ2_2_Name, VBFsigRateUp_mZZ2_2_Name, integral_VBF_T_mZZ2_2_2_Up)
      VBFinterfRates_mZZ2_2_Up = ROOT.RooRealVar(VBFinterfRateUp_mZZ2_2_Name, VBFinterfRateUp_mZZ2_2_Name, integral_VBF_T_mZZ2_2_4_Up)
      VBFsigRates_mZZ2_2_Up.setConstant(True)
      VBFinterfRates_mZZ2_2_Up.setConstant(True)

      VBFsigRateDown_mZZ2_2_Name = "signal_VBFDownrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateDown_mZZ2_2_Name = "interf_VBFDownrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_2_Down = ROOT.RooRealVar(VBFsigRateDown_mZZ2_2_Name, VBFsigRateDown_mZZ2_2_Name, integral_VBF_T_mZZ2_2_2_Down)
      VBFinterfRates_mZZ2_2_Down = ROOT.RooRealVar(VBFinterfRateDown_mZZ2_2_Name, VBFinterfRateDown_mZZ2_2_Name, integral_VBF_T_mZZ2_2_4_Down)
      VBFsigRates_mZZ2_2_Down.setConstant(True)
      VBFinterfRates_mZZ2_2_Down.setConstant(True)

      VBFsigRateNominal_mZZ2_2_Name = "signal_VBFNominalrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateNominal_mZZ2_2_Name = "interf_VBFNominalrate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_2_Nominal = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_2_Name, VBFsigRateNominal_mZZ2_2_Name, integral_VBF_T_mZZ2_2_2)
      VBFinterfRates_mZZ2_2_Nominal = ROOT.RooRealVar(VBFinterfRateNominal_mZZ2_2_Name, VBFinterfRateNominal_mZZ2_2_Name, integral_VBF_T_mZZ2_2_4)
      VBFsigRates_mZZ2_2_Nominal.setConstant(True)
      VBFinterfRates_mZZ2_2_Nominal.setConstant(True)

      # VBF mZZ/mH**2**3

      VBFsigRateUp_mZZ2_3_Name = "signal_VBFUprate_mZZ2_3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_3_Up = ROOT.RooRealVar(VBFsigRateUp_mZZ2_3_Name, VBFsigRateUp_mZZ2_3_Name, integral_VBF_T_mZZ2_3_2_Up)
      VBFsigRates_mZZ2_3_Up.setConstant(True)

      VBFsigRateDown_mZZ2_3_Name = "signal_VBFDownrate_mZZ2_3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_3_Down = ROOT.RooRealVar(VBFsigRateDown_mZZ2_3_Name, VBFsigRateDown_mZZ2_3_Name, integral_VBF_T_mZZ2_3_2_Down)
      VBFsigRates_mZZ2_3_Down.setConstant(True)

      VBFsigRateNominal_mZZ2_3_Name = "signal_VBFNominalrate_mZZ2_3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_3_Nominal = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_3_Name, VBFsigRateNominal_mZZ2_3_Name, integral_VBF_T_mZZ2_3_2)
      VBFsigRates_mZZ2_3_Nominal.setConstant(True)

      # VBF mZZ/mH**2**4

      VBFsigRateUp_mZZ2_4_Name = "signal_VBFUprate_mZZ2_4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_4_Up = ROOT.RooRealVar(VBFsigRateUp_mZZ2_4_Name, VBFsigRateUp_mZZ2_4_Name, integral_VBF_T_mZZ2_4_2_Up)
      VBFsigRates_mZZ2_4_Up.setConstant(True)

      VBFsigRateDown_mZZ2_4_Name = "signal_VBFDownrate_mZZ2_4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_4_Down = ROOT.RooRealVar(VBFsigRateDown_mZZ2_4_Name, VBFsigRateDown_mZZ2_4_Name, integral_VBF_T_mZZ2_4_2_Down)
      VBFsigRates_mZZ2_4_Down.setConstant(True)

      VBFsigRateNominal_mZZ2_4_Name = "signal_VBFNominalrate_mZZ2_4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_4_Nominal = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_4_Name, VBFsigRateNominal_mZZ2_4_Name, integral_VBF_T_mZZ2_4_2)
      VBFsigRates_mZZ2_4_Nominal.setConstant(True)


      #------------- VBF Opposite Djet --------------#

      VBFsigRateUpName_OppositeDjet = "signal_VBFUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateUpName_OppositeDjet = "bkg_VBFUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateUpName_OppositeDjet = "interf_VBFUp_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Up_OppositeDjet = ROOT.RooRealVar(VBFsigRateUpName_OppositeDjet, VBFsigRateUpName_OppositeDjet, integral_VBF_T_2_Up_OppositeDjet)
      VBFbkgRates_Up_OppositeDjet = ROOT.RooRealVar(VBFbkgRateUpName_OppositeDjet, VBFbkgRateUpName_OppositeDjet, integral_VBF_T_1_Up_OppositeDjet)
      VBFinterfRates_Up_OppositeDjet = ROOT.RooRealVar(VBFinterfRateUpName_OppositeDjet, VBFinterfRateUpName_OppositeDjet, integral_VBF_T_4_Up_OppositeDjet)
      VBFsigRates_Up_OppositeDjet.setConstant(True)
      VBFbkgRates_Up_OppositeDjet.setConstant(True)
      VBFinterfRates_Up_OppositeDjet.setConstant(True)

      VBFsigRateDownName_OppositeDjet = "signal_VBFDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateDownName_OppositeDjet = "bkg_VBFDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateDownName_OppositeDjet = "interf_VBFDown_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Down_OppositeDjet = ROOT.RooRealVar(VBFsigRateDownName_OppositeDjet, VBFsigRateDownName_OppositeDjet, integral_VBF_T_2_Down_OppositeDjet)
      VBFbkgRates_Down_OppositeDjet = ROOT.RooRealVar(VBFbkgRateDownName_OppositeDjet, VBFbkgRateDownName_OppositeDjet, integral_VBF_T_1_Down_OppositeDjet)
      VBFinterfRates_Down_OppositeDjet = ROOT.RooRealVar(VBFinterfRateDownName_OppositeDjet, VBFinterfRateDownName_OppositeDjet, integral_VBF_T_4_Down_OppositeDjet)
      VBFsigRates_Down_OppositeDjet.setConstant(True)
      VBFbkgRates_Down_OppositeDjet.setConstant(True)
      VBFinterfRates_Down_OppositeDjet.setConstant(True)

      VBFsigRateNominalName_OppositeDjet = "signal_VBFNominal_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateNominalName_OppositeDjet = "bkg_VBFNominal_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateNominalName_OppositeDjet = "interf_VBFNominal_OppositeDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Nominal_OppositeDjet = ROOT.RooRealVar(VBFsigRateNominalName_OppositeDjet, VBFsigRateNominalName_OppositeDjet, integral_VBF_T_2_OppositeDjet)
      VBFbkgRates_Nominal_OppositeDjet = ROOT.RooRealVar(VBFbkgRateNominalName_OppositeDjet, VBFbkgRateNominalName_OppositeDjet, integral_VBF_T_1_OppositeDjet)
      VBFinterfRates_Nominal_OppositeDjet = ROOT.RooRealVar(VBFinterfRateNominalName_OppositeDjet, VBFinterfRateNominalName_OppositeDjet, integral_VBF_T_4_OppositeDjet)
      VBFsigRates_Nominal_OppositeDjet.setConstant(True)
      VBFbkgRates_Nominal_OppositeDjet.setConstant(True)
      VBFinterfRates_Nominal_OppositeDjet.setConstant(True)


      # VBF mZZ/mH**2

      VBFsigRateUp_mZZ2_1_Name_OppositeDjet = "signal_VBFUp_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateUp_mZZ2_1_Name_OppositeDjet = "interf_VBFUp_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_1_Up_OppositeDjet = ROOT.RooRealVar(VBFsigRateUp_mZZ2_1_Name_OppositeDjet, VBFsigRateUp_mZZ2_1_Name_OppositeDjet, integral_VBF_T_mZZ2_1_2_Up_OppositeDjet)
      VBFinterfRates_mZZ2_1_Up_OppositeDjet = ROOT.RooRealVar(VBFinterfRateUp_mZZ2_1_Name_OppositeDjet, VBFinterfRateUp_mZZ2_1_Name_OppositeDjet, integral_VBF_T_mZZ2_1_4_Up_OppositeDjet)
      VBFsigRates_mZZ2_1_Up_OppositeDjet.setConstant(True)
      VBFinterfRates_mZZ2_1_Up_OppositeDjet.setConstant(True)

      VBFsigRateDown_mZZ2_1_Name_OppositeDjet = "signal_VBFDown_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateDown_mZZ2_1_Name_OppositeDjet = "interf_VBFDown_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_1_Down_OppositeDjet = ROOT.RooRealVar(VBFsigRateDown_mZZ2_1_Name_OppositeDjet, VBFsigRateDown_mZZ2_1_Name_OppositeDjet, integral_VBF_T_mZZ2_1_2_Down_OppositeDjet)
      VBFinterfRates_mZZ2_1_Down_OppositeDjet = ROOT.RooRealVar(VBFinterfRateDown_mZZ2_1_Name_OppositeDjet, VBFinterfRateDown_mZZ2_1_Name_OppositeDjet, integral_VBF_T_mZZ2_1_4_Down_OppositeDjet)
      VBFsigRates_mZZ2_1_Down_OppositeDjet.setConstant(True)
      VBFinterfRates_mZZ2_1_Down_OppositeDjet.setConstant(True)

      VBFsigRateNominal_mZZ2_1_Name_OppositeDjet = "signal_VBFNominal_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateNominal_mZZ2_1_Name_OppositeDjet = "interf_VBFNominal_OppositeDjet_rate_mZZ2_1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_1_Nominal_OppositeDjet = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_1_Name_OppositeDjet, VBFsigRateNominal_mZZ2_1_Name_OppositeDjet, integral_VBF_T_mZZ2_1_2_OppositeDjet)
      VBFinterfRates_mZZ2_1_Nominal_OppositeDjet = ROOT.RooRealVar(VBFinterfRateNominal_mZZ2_1_Name_OppositeDjet, VBFinterfRateNominal_mZZ2_1_Name_OppositeDjet, integral_VBF_T_mZZ2_1_4_OppositeDjet)
      VBFsigRates_mZZ2_1_Nominal_OppositeDjet.setConstant(True)
      VBFinterfRates_mZZ2_1_Nominal_OppositeDjet.setConstant(True)

      # VBF mZZ/mH**2**2

      VBFsigRateUp_mZZ2_2_Name_OppositeDjet = "signal_VBFUp_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateUp_mZZ2_2_Name_OppositeDjet = "interf_VBFUp_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_2_Up_OppositeDjet = ROOT.RooRealVar(VBFsigRateUp_mZZ2_2_Name_OppositeDjet, VBFsigRateUp_mZZ2_2_Name_OppositeDjet, integral_VBF_T_mZZ2_2_2_Up_OppositeDjet)
      VBFinterfRates_mZZ2_2_Up_OppositeDjet = ROOT.RooRealVar(VBFinterfRateUp_mZZ2_2_Name_OppositeDjet, VBFinterfRateUp_mZZ2_2_Name_OppositeDjet, integral_VBF_T_mZZ2_2_4_Up_OppositeDjet)
      VBFsigRates_mZZ2_2_Up_OppositeDjet.setConstant(True)
      VBFinterfRates_mZZ2_2_Up_OppositeDjet.setConstant(True)

      VBFsigRateDown_mZZ2_2_Name_OppositeDjet = "signal_VBFDown_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateDown_mZZ2_2_Name_OppositeDjet = "interf_VBFDown_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_2_Down_OppositeDjet = ROOT.RooRealVar(VBFsigRateDown_mZZ2_2_Name_OppositeDjet, VBFsigRateDown_mZZ2_2_Name_OppositeDjet, integral_VBF_T_mZZ2_2_2_Down_OppositeDjet)
      VBFinterfRates_mZZ2_2_Down_OppositeDjet = ROOT.RooRealVar(VBFinterfRateDown_mZZ2_2_Name_OppositeDjet, VBFinterfRateDown_mZZ2_2_Name_OppositeDjet, integral_VBF_T_mZZ2_2_4_Down_OppositeDjet)
      VBFsigRates_mZZ2_2_Down_OppositeDjet.setConstant(True)
      VBFinterfRates_mZZ2_2_Down_OppositeDjet.setConstant(True)

      VBFsigRateNominal_mZZ2_2_Name_OppositeDjet = "signal_VBFNominal_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRateNominal_mZZ2_2_Name_OppositeDjet = "interf_VBFNominal_OppositeDjet_rate_mZZ2_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_2_Nominal_OppositeDjet = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_2_Name_OppositeDjet, VBFsigRateNominal_mZZ2_2_Name_OppositeDjet, integral_VBF_T_mZZ2_2_2_OppositeDjet)
      VBFinterfRates_mZZ2_2_Nominal_OppositeDjet = ROOT.RooRealVar(VBFinterfRateNominal_mZZ2_2_Name_OppositeDjet, VBFinterfRateNominal_mZZ2_2_Name_OppositeDjet, integral_VBF_T_mZZ2_2_4_OppositeDjet)
      VBFsigRates_mZZ2_2_Nominal_OppositeDjet.setConstant(True)
      VBFinterfRates_mZZ2_2_Nominal_OppositeDjet.setConstant(True)

      # VBF mZZ/mH**2**3

      VBFsigRateUp_mZZ2_3_Name_OppositeDjet = "signal_VBFUp_OppositeDjet_rate_mZZ2_3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_3_Up_OppositeDjet = ROOT.RooRealVar(VBFsigRateUp_mZZ2_3_Name_OppositeDjet, VBFsigRateUp_mZZ2_3_Name_OppositeDjet, integral_VBF_T_mZZ2_3_2_Up_OppositeDjet)
      VBFsigRates_mZZ2_3_Up_OppositeDjet.setConstant(True)

      VBFsigRateDown_mZZ2_3_Name_OppositeDjet = "signal_VBFDown_OppositeDjet_rate_mZZ2_3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_3_Down_OppositeDjet = ROOT.RooRealVar(VBFsigRateDown_mZZ2_3_Name_OppositeDjet, VBFsigRateDown_mZZ2_3_Name_OppositeDjet, integral_VBF_T_mZZ2_3_2_Down_OppositeDjet)
      VBFsigRates_mZZ2_3_Down_OppositeDjet.setConstant(True)

      VBFsigRateNominal_mZZ2_3_Name_OppositeDjet = "signal_VBFNominal_OppositeDjet_rate_mZZ2_3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_3_Nominal_OppositeDjet = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_3_Name_OppositeDjet, VBFsigRateNominal_mZZ2_3_Name_OppositeDjet, integral_VBF_T_mZZ2_3_2_OppositeDjet)
      VBFsigRates_mZZ2_3_Nominal_OppositeDjet.setConstant(True)

      # VBF mZZ/mH**2**4

      VBFsigRateUp_mZZ2_4_Name_OppositeDjet = "signal_VBFUp_OppositeDjet_rate_mZZ2_4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_4_Up_OppositeDjet = ROOT.RooRealVar(VBFsigRateUp_mZZ2_4_Name_OppositeDjet, VBFsigRateUp_mZZ2_4_Name_OppositeDjet, integral_VBF_T_mZZ2_4_2_Up_OppositeDjet)
      VBFsigRates_mZZ2_4_Up_OppositeDjet.setConstant(True)

      VBFsigRateDown_mZZ2_4_Name_OppositeDjet = "signal_VBFDown_OppositeDjet_rate_mZZ2_4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_4_Down_OppositeDjet = ROOT.RooRealVar(VBFsigRateDown_mZZ2_4_Name_OppositeDjet, VBFsigRateDown_mZZ2_4_Name_OppositeDjet, integral_VBF_T_mZZ2_4_2_Down_OppositeDjet)
      VBFsigRates_mZZ2_4_Down_OppositeDjet.setConstant(True)

      VBFsigRateNominal_mZZ2_4_Name_OppositeDjet = "signal_VBFNominal_OppositeDjet_rate_mZZ2_4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_mZZ2_4_Nominal_OppositeDjet = ROOT.RooRealVar(VBFsigRateNominal_mZZ2_4_Name_OppositeDjet, VBFsigRateNominal_mZZ2_4_Name_OppositeDjet, integral_VBF_T_mZZ2_4_2_OppositeDjet)
      VBFsigRates_mZZ2_4_Nominal_OppositeDjet.setConstant(True)


      # VBF anomalous coupling rate parameterizations

      VBFsigRates_Nominal_AnomCoupl_Name = "signal_VBFNominalrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Nominal_AnomCoupl = ROOT.RooFormulaVar(
      VBFsigRates_Nominal_AnomCoupl_Name,
      "( pow((1-abs(@5)),2)*@0 + sign(@5)*sqrt(abs(@5))*pow(sqrt(1-abs(@5)),3)*@1 + abs(@5)*(1-abs(@5))*@2 + sign(@5)*pow(sqrt(abs(@5)),3)*sqrt(1-abs(@5))*@3 + pow(@5,2)*@4 )",
      ROOT.RooArgList(VBFsigRates_Nominal,VBFsigRates_mZZ2_1_Nominal,VBFsigRates_mZZ2_2_Nominal,VBFsigRates_mZZ2_3_Nominal,VBFsigRates_mZZ2_4_Nominal, fai1)
      )
      VBFinterfRates_Nominal_AnomCoupl_Name = "interf_VBFNominalrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Nominal_AnomCoupl = ROOT.RooFormulaVar(
      VBFinterfRates_Nominal_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(VBFinterfRates_Nominal,VBFinterfRates_mZZ2_1_Nominal,VBFinterfRates_mZZ2_2_Nominal, fai1)
      )

      VBFsigRates_Up_AnomCoupl_Name = "signal_VBFUprate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Up_AnomCoupl = ROOT.RooFormulaVar(
      VBFsigRates_Up_AnomCoupl_Name,
      "( pow((1-abs(@5)),2)*@0 + sign(@5)*sqrt(abs(@5))*pow(sqrt(1-abs(@5)),3)*@1 + abs(@5)*(1-abs(@5))*@2 + sign(@5)*pow(sqrt(abs(@5)),3)*sqrt(1-abs(@5))*@3 + pow(@5,2)*@4 )",
      ROOT.RooArgList(VBFsigRates_Up,VBFsigRates_mZZ2_1_Up,VBFsigRates_mZZ2_2_Up,VBFsigRates_mZZ2_3_Up,VBFsigRates_mZZ2_4_Up, fai1)
      )
      VBFinterfRates_Up_AnomCoupl_Name = "interf_VBFUprate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Up_AnomCoupl = ROOT.RooFormulaVar(
      VBFinterfRates_Up_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(VBFinterfRates_Up,VBFinterfRates_mZZ2_1_Up,VBFinterfRates_mZZ2_2_Up, fai1)
      )

      VBFsigRates_Down_AnomCoupl_Name = "signal_VBFDownrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Down_AnomCoupl = ROOT.RooFormulaVar(
      VBFsigRates_Down_AnomCoupl_Name,
      "( pow((1-abs(@5)),2)*@0 + sign(@5)*sqrt(abs(@5))*pow(sqrt(1-abs(@5)),3)*@1 + abs(@5)*(1-abs(@5))*@2 + sign(@5)*pow(sqrt(abs(@5)),3)*sqrt(1-abs(@5))*@3 + pow(@5,2)*@4 )",
      ROOT.RooArgList(VBFsigRates_Down,VBFsigRates_mZZ2_1_Down,VBFsigRates_mZZ2_2_Down,VBFsigRates_mZZ2_3_Down,VBFsigRates_mZZ2_4_Down, fai1)
      )
      VBFinterfRates_Down_AnomCoupl_Name = "interf_VBFDownrate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Down_AnomCoupl = ROOT.RooFormulaVar(
      VBFinterfRates_Down_AnomCoupl_Name,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(VBFinterfRates_Down,VBFinterfRates_mZZ2_1_Down,VBFinterfRates_mZZ2_2_Down, fai1)
      )


      VBFsigRates_Nominal_AnomCoupl_Name_OppositeDjet = "signal_VBFNominal_OppositeDjet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Nominal_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      VBFsigRates_Nominal_AnomCoupl_Name_OppositeDjet,
      "( pow((1-abs(@5)),2)*@0 + sign(@5)*sqrt(abs(@5))*pow(sqrt(1-abs(@5)),3)*@1 + abs(@5)*(1-abs(@5))*@2 + sign(@5)*pow(sqrt(abs(@5)),3)*sqrt(1-abs(@5))*@3 + pow(@5,2)*@4 )",
      ROOT.RooArgList(VBFsigRates_Nominal_OppositeDjet,VBFsigRates_mZZ2_1_Nominal_OppositeDjet,VBFsigRates_mZZ2_2_Nominal_OppositeDjet,VBFsigRates_mZZ2_3_Nominal_OppositeDjet,VBFsigRates_mZZ2_4_Nominal_OppositeDjet, fai1)
      )
      VBFinterfRates_Nominal_AnomCoupl_Name_OppositeDjet = "interf_VBFNominal_OppositeDjet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Nominal_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      VBFinterfRates_Nominal_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(VBFinterfRates_Nominal_OppositeDjet,VBFinterfRates_mZZ2_1_Nominal_OppositeDjet,VBFinterfRates_mZZ2_2_Nominal_OppositeDjet, fai1)
      )

      VBFsigRates_Up_AnomCoupl_Name_OppositeDjet = "signal_VBFUp_OppositeDjet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Up_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      VBFsigRates_Up_AnomCoupl_Name_OppositeDjet,
      "( pow((1-abs(@5)),2)*@0 + sign(@5)*sqrt(abs(@5))*pow(sqrt(1-abs(@5)),3)*@1 + abs(@5)*(1-abs(@5))*@2 + sign(@5)*pow(sqrt(abs(@5)),3)*sqrt(1-abs(@5))*@3 + pow(@5,2)*@4 )",
      ROOT.RooArgList(VBFsigRates_Up_OppositeDjet,VBFsigRates_mZZ2_1_Up_OppositeDjet,VBFsigRates_mZZ2_2_Up_OppositeDjet,VBFsigRates_mZZ2_3_Up_OppositeDjet,VBFsigRates_mZZ2_4_Up_OppositeDjet, fai1)
      )
      VBFinterfRates_Up_AnomCoupl_Name_OppositeDjet = "interf_VBFUp_OppositeDjet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Up_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      VBFinterfRates_Up_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(VBFinterfRates_Up_OppositeDjet,VBFinterfRates_mZZ2_1_Up_OppositeDjet,VBFinterfRates_mZZ2_2_Up_OppositeDjet, fai1)
      )

      VBFsigRates_Down_AnomCoupl_Name_OppositeDjet = "signal_VBFDown_OppositeDjet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Down_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      VBFsigRates_Down_AnomCoupl_Name_OppositeDjet,
      "( pow((1-abs(@5)),2)*@0 + sign(@5)*sqrt(abs(@5))*pow(sqrt(1-abs(@5)),3)*@1 + abs(@5)*(1-abs(@5))*@2 + sign(@5)*pow(sqrt(abs(@5)),3)*sqrt(1-abs(@5))*@3 + pow(@5,2)*@4 )",
      ROOT.RooArgList(VBFsigRates_Down_OppositeDjet,VBFsigRates_mZZ2_1_Down_OppositeDjet,VBFsigRates_mZZ2_2_Down_OppositeDjet,VBFsigRates_mZZ2_3_Down_OppositeDjet,VBFsigRates_mZZ2_4_Down_OppositeDjet, fai1)
      )
      VBFinterfRates_Down_AnomCoupl_Name_OppositeDjet = "interf_VBFDown_OppositeDjet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Down_AnomCoupl_OppositeDjet = ROOT.RooFormulaVar(
      VBFinterfRates_Down_AnomCoupl_Name_OppositeDjet,
      "( (1-abs(@3))*@0 + sign(@3)*sqrt(abs(@3)*(1-abs(@3)))*@1 + abs(@3)*@2 )",
      ROOT.RooArgList(VBFinterfRates_Down_OppositeDjet,VBFinterfRates_mZZ2_1_Down_OppositeDjet,VBFinterfRates_mZZ2_2_Down_OppositeDjet, fai1)
      )


      #------- Combined VBF Djet ----------------#

      VBFsigRates_Nominal_AnomCoupl_CombinedJet_Name = "signal_VBFNominal_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet_Name = "interf_VBFNominal_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateNominalName_CombinedDjet = "bkg_VBFNominal_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(VBFsigRates_Nominal_AnomCoupl)
      )
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(VBFinterfRates_Nominal_AnomCoupl)
      )
      VBFbkgRates_Nominal_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateNominalName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(VBFbkgRates_Nominal)
      )
      if useDjet == 1:
      VBFsigRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFsigRates_Nominal_AnomCoupl, VBFsigRates_Nominal_AnomCoupl_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFinterfRates_Nominal_AnomCoupl, VBFinterfRates_Nominal_AnomCoupl_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      VBFbkgRates_Nominal_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateNominalName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFbkgRates_Nominal, VBFbkgRates_Nominal_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      if useDjet == 2:
      VBFsigRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFsigRates_Nominal_AnomCoupl, thetaSyst_djet_VBF_norm)
      )
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Nominal_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFinterfRates_Nominal_AnomCoupl, thetaSyst_djet_VBF_norm)
      )
      VBFbkgRates_Nominal_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateNominalName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFbkgRates_Nominal, thetaSyst_djet_VBF_norm)
      )


      VBFsigRates_Up_AnomCoupl_CombinedJet_Name = "signal_VBFUp_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Up_AnomCoupl_CombinedJet_Name = "interf_VBFUp_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateUpName_CombinedDjet = "bkg_VBFUp_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Up_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Up_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(VBFsigRates_Up_AnomCoupl)
      )
      VBFinterfRates_Up_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Up_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(VBFinterfRates_Up_AnomCoupl)
      )
      VBFbkgRates_Up_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateUpName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(VBFbkgRates_Up)
      )
      if useDjet == 1:
      VBFsigRates_Up_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Up_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFsigRates_Up_AnomCoupl, VBFsigRates_Up_AnomCoupl_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      VBFinterfRates_Up_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Up_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFinterfRates_Up_AnomCoupl, VBFinterfRates_Up_AnomCoupl_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      VBFbkgRates_Up_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateUpName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFbkgRates_Up, VBFbkgRates_Up_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      if useDjet == 2:
      VBFsigRates_Up_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Up_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFsigRates_Up_AnomCoupl, thetaSyst_djet_VBF_norm)
      )
      VBFinterfRates_Up_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Up_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFinterfRates_Up_AnomCoupl, thetaSyst_djet_VBF_norm)
      )
      VBFbkgRates_Up_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateUpName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFbkgRates_Up, thetaSyst_djet_VBF_norm)
      )


      VBFsigRates_Down_AnomCoupl_CombinedJet_Name = "signal_VBFDown_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFinterfRates_Down_AnomCoupl_CombinedJet_Name = "interf_VBFDown_CombinedJet_rate_AnomCoupl_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFbkgRateDownName_CombinedDjet = "bkg_VBFDown_CombinedDjet_rate_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFsigRates_Down_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Down_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(VBFsigRates_Down_AnomCoupl)
      )
      VBFinterfRates_Down_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Down_AnomCoupl_CombinedJet_Name,
      "( @0 )",
      ROOT.RooArgList(VBFinterfRates_Down_AnomCoupl)
      )
      VBFbkgRates_Down_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateDownName_CombinedDjet,
      "( @0 )",
      ROOT.RooArgList(VBFbkgRates_Down)
      )
      if useDjet == 1:
      VBFsigRates_Down_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Down_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFsigRates_Down_AnomCoupl, VBFsigRates_Down_AnomCoupl_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      VBFinterfRates_Down_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Down_AnomCoupl_CombinedJet_Name,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFinterfRates_Down_AnomCoupl, VBFinterfRates_Down_AnomCoupl_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      VBFbkgRates_Down_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateDownName_CombinedDjet,
      "( @0 + ( 1-TMath::Max(@2,0) )*@1 )",
      ROOT.RooArgList(VBFbkgRates_Down, VBFbkgRates_Down_OppositeDjet, thetaSyst_djet_VBF_norm)
      )
      if useDjet == 2:
      VBFsigRates_Down_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFsigRates_Down_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFsigRates_Down_AnomCoupl, thetaSyst_djet_VBF_norm)
      )
      VBFinterfRates_Down_AnomCoupl_CombinedJet = ROOT.RooFormulaVar(
      VBFinterfRates_Down_AnomCoupl_CombinedJet_Name,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFinterfRates_Down_AnomCoupl, thetaSyst_djet_VBF_norm)
      )
      VBFbkgRates_Down_CombinedJet = ROOT.RooFormulaVar(
      VBFbkgRateDownName_CombinedDjet,
      "( @0*TMath::Max(@1,0) )",
      ROOT.RooArgList(VBFbkgRates_Down, thetaSyst_djet_VBF_norm)
      )


      VBFVarNormNominal_Name = "VBFVarNominalNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFNominal_norm = ROOT.RooFormulaVar(
      VBFVarNormNominal_Name, "(@0*@3*@5*@4+@1*sqrt(@3*@5*@4)+@2)",
      ROOT.RooArgList(VBFsigRates_Nominal_AnomCoupl_CombinedJet, VBFinterfRates_Nominal_AnomCoupl_CombinedJet, VBFbkgRates_Nominal_CombinedJet, x, mu, muV)
      )
      VBFVarNormUp_Name = "VBFVarUpNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFUp_norm = ROOT.RooFormulaVar(
      VBFVarNormUp_Name, "(@0*@3*@5*@4+@1*sqrt(@3*@5*@4)+@2)",
      ROOT.RooArgList(VBFsigRates_Up_AnomCoupl_CombinedJet, VBFinterfRates_Up_AnomCoupl_CombinedJet, VBFbkgRates_Up_CombinedJet, x, mu, muV)
      )
      VBFVarNormDown_Name = "VBFVarDownNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFDown_norm = ROOT.RooFormulaVar(
      VBFVarNormDown_Name, "(@0*@3*@5*@4+@1*sqrt(@3*@5*@4)+@2)",
      ROOT.RooArgList(VBFsigRates_Down_AnomCoupl_CombinedJet, VBFinterfRates_Down_AnomCoupl_CombinedJet, VBFbkgRates_Down_CombinedJet, x, mu, muV)
      )



      # normalization on background and protection against negative
      # fluctuations
      for ix in range(0, Bkg_T.GetXaxis().GetNbins() + 2):
      yNorm = Bkg_T.Integral(ix, ix, 1, Bkg_T.GetYaxis().GetNbins())
      yNorm_zx = Bkg_ZX.Integral(ix, ix, 1, Bkg_ZX.GetYaxis().GetNbins())
      yNorm_zx_Up = Bkg_ZX_Up.Integral(ix, ix, 1, Bkg_ZX_Up.GetYaxis().GetNbins())
      yNorm_zx_Down = Bkg_ZX_Down.Integral(ix, ix, 1, Bkg_ZX_Down.GetYaxis().GetNbins())
      for iy in range(0, Bkg_T.GetYaxis().GetNbins() + 2):
      if yNorm != 0:
      Bkg_T.SetBinContent(ix, iy, Bkg_T.GetBinContent(ix, iy) / yNorm)
      if yNorm_zx != 0:
      Bkg_ZX.SetBinContent(ix, iy, Bkg_ZX.GetBinContent(ix, iy) / yNorm_zx)
      if yNorm_zx_Up != 0:
      Bkg_ZX_Up.SetBinContent(ix, iy, Bkg_ZX_Up.GetBinContent(ix, iy) / yNorm_zx_Up)
      if yNorm_zx_Down != 0:
      Bkg_ZX_Down.SetBinContent(ix, iy, Bkg_ZX_Down.GetBinContent(ix, iy) / yNorm_zx_Down)

      dBinsX = Sig_T_1.GetXaxis().GetNbins()
      print "X bins: ", dBinsX

      dBinsY = Sig_T_1.GetYaxis().GetNbins()
      print "Y bins: ", dBinsY

      CMS_zz4l_widthMass.setBins(dBinsX)
      CMS_zz4l_widthKD.setBins(dBinsY)

      # -------------------------- gg2ZZ SHAPES ---------------------------------- ##
      sigRatesNormList=[]
      interfRatesNormList=[]
      bkgRatesNormList=[]

      bkg_ggZZ_RawHistList=[Sig_T_1]
      bkg_ggZZ_DataHistList=[]
      bkg_ggZZ_HistFuncList=[]

      bkg_ggZZ_RawHistPDFUpList=[Sig_T_1_Up_PDF]
      bkg_ggZZ_DataHistPDFUpList=[]
      bkg_ggZZ_HistFuncPDFUpList=[]

      bkg_ggZZ_RawHistPDFDownList=[Sig_T_1_Down_PDF]
      bkg_ggZZ_DataHistPDFDownList=[]
      bkg_ggZZ_HistFuncPDFDownList=[]

      bkg_ggZZ_RawHistQCDUpList=[Sig_T_1_Up_QCD]
      bkg_ggZZ_DataHistQCDUpList=[]
      bkg_ggZZ_HistFuncQCDUpList=[]

      bkg_ggZZ_RawHistQCDDownList=[Sig_T_1_Down_QCD]
      bkg_ggZZ_DataHistQCDDownList=[]
      bkg_ggZZ_HistFuncQCDDownList=[]


      signal_ggZZ_RawHistList=[Sig_T_2,Sig_T_mZZ2_1_2,Sig_T_mZZ2_2_2]
      signal_ggZZ_DataHistList=[]
      signal_ggZZ_HistFuncList=[]

      signal_ggZZ_RawHistPDFUpList=[Sig_T_2_Up_PDF,Sig_T_mZZ2_1_2_Up_PDF,Sig_T_mZZ2_2_2_Up_PDF]
      signal_ggZZ_DataHistPDFUpList=[]
      signal_ggZZ_HistFuncPDFUpList=[]

      signal_ggZZ_RawHistPDFDownList=[Sig_T_2_Down_PDF,Sig_T_mZZ2_1_2_Down_PDF,Sig_T_mZZ2_2_2_Down_PDF]
      signal_ggZZ_DataHistPDFDownList=[]
      signal_ggZZ_HistFuncPDFDownList=[]

      signal_ggZZ_RawHistQCDUpList=[Sig_T_2_Up_QCD,Sig_T_mZZ2_1_2_Up_QCD,Sig_T_mZZ2_2_2_Up_QCD]
      signal_ggZZ_DataHistQCDUpList=[]
      signal_ggZZ_HistFuncQCDUpList=[]

      signal_ggZZ_RawHistQCDDownList=[Sig_T_2_Down_QCD,Sig_T_mZZ2_1_2_Down_QCD,Sig_T_mZZ2_2_2_Down_QCD]
      signal_ggZZ_DataHistQCDDownList=[]
      signal_ggZZ_HistFuncQCDDownList=[]


      interf_ggZZ_RawHistList=[Sig_T_4,Sig_T_mZZ2_1_4]
      interf_ggZZ_DataHistList=[]
      interf_ggZZ_HistFuncList=[]

      interf_ggZZ_RawHistPDFUpList=[Sig_T_4_Up_PDF,Sig_T_mZZ2_1_4_Up_PDF]
      interf_ggZZ_DataHistPDFUpList=[]
      interf_ggZZ_HistFuncPDFUpList=[]

      interf_ggZZ_RawHistPDFDownList=[Sig_T_4_Down_PDF,Sig_T_mZZ2_1_4_Down_PDF]
      interf_ggZZ_DataHistPDFDownList=[]
      interf_ggZZ_HistFuncPDFDownList=[]

      interf_ggZZ_RawHistQCDUpList=[Sig_T_4_Up_QCD,Sig_T_mZZ2_1_4_Up_QCD]
      interf_ggZZ_DataHistQCDUpList=[]
      interf_ggZZ_HistFuncQCDUpList=[]

      interf_ggZZ_RawHistQCDDownList=[Sig_T_4_Down_QCD,Sig_T_mZZ2_1_4_Down_QCD]
      interf_ggZZ_DataHistQCDDownList=[]
      interf_ggZZ_HistFuncQCDDownList=[]


      bkgRateNameNorm = "bkgNorm_ggZZrate"
      bkgRatesNorm = ROOT.RooFormulaVar(bkgRateNameNorm, "@0", ROOT.RooArgList(kbkg_gg))

      sigRateNameWidthNorm = "signalWidthNorm_ggZZrate"
      interfRateNameWidthNorm = "interfWidthNorm_ggZZrate"

      sigRatesWidthNorm = ROOT.RooFormulaVar(sigRateNameWidthNorm, "@0*@1*@2", ROOT.RooArgList(x, mu, muF))
      interfRatesWidthNorm = ROOT.RooFormulaVar(interfRateNameWidthNorm, "sqrt(@0*@1*@2)*sign(@3)*sqrt(abs(@3))", ROOT.RooArgList(x, mu, muF, kbkg_gg))


      # ggZZ Bkg histfunc construction
      TemplateName = "ggZZbkg_TempDataHist_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "ggZZbkg_TempHistFunc_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_ggZZ_RawHistList[0])
      bkg_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_ggZZ_RawHistList[0].ProjectionX())
      bkg_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_ggZZ_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_ggZZ_RawHistList[0].ProjectionY())
      bkg_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncList.append(TempHistFunc)
      TemplateName = "ggZZbkg_TempDataHist_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "ggZZbkg_TempHistFunc_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_ggZZ_RawHistQCDDownList[0])
      bkg_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_ggZZ_RawHistQCDDownList[0].ProjectionX())
      bkg_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_ggZZ_RawHistQCDDownList[0].ProjectionY())
      bkg_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      TemplateName = "ggZZbkg_TempDataHist_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "ggZZbkg_TempHistFunc_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_ggZZ_RawHistQCDUpList[0])
      bkg_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_ggZZ_RawHistQCDUpList[0].ProjectionX())
      bkg_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_ggZZ_RawHistQCDUpList[0].ProjectionY())
      bkg_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      TemplateName = "ggZZbkg_TempDataHist_Down_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "ggZZbkg_TempHistFunc_Down_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_ggZZ_RawHistPDFDownList[0])
      bkg_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_ggZZ_RawHistPDFDownList[0].ProjectionX())
      bkg_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_ggZZ_RawHistPDFDownList[0].ProjectionY())
      bkg_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      TemplateName = "ggZZbkg_TempDataHist_Up_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "ggZZbkg_TempHistFunc_Up_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_ggZZ_RawHistPDFUpList[0])
      bkg_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_ggZZ_RawHistPDFUpList[0].ProjectionX())
      bkg_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_ggZZ_RawHistPDFUpList[0].ProjectionY())
      bkg_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_ggZZ_HistFuncPDFUpList.append(TempHistFunc)


      anomalousLoops = 1
      anomalousLoops_interf = 1
      if self.anomCoupl == 1:
      anomalousLoops = 3 # For ggZZ
      anomalousLoops_interf = 2

      for al in range(0,anomalousLoops) :
      TemplateName = "ggZZsignal_TempDataHist_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZsignal_TempHistFunc_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_ggZZ_RawHistList[al])
      signal_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_ggZZ_RawHistList[al].ProjectionX())
      signal_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_ggZZ_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_ggZZ_RawHistList[al].ProjectionY())
      signal_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncList.append(TempHistFunc)
      TemplateName = "ggZZsignal_TempDataHist_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZsignal_TempHistFunc_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_ggZZ_RawHistQCDUpList[al])
      signal_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_ggZZ_RawHistQCDUpList[al].ProjectionX())
      signal_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_ggZZ_RawHistQCDUpList[al].ProjectionY())
      signal_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      TemplateName = "ggZZsignal_TempDataHist_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZsignal_TempHistFunc_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_ggZZ_RawHistQCDDownList[al])
      signal_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_ggZZ_RawHistQCDDownList[al].ProjectionX())
      signal_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_ggZZ_RawHistQCDDownList[al].ProjectionY())
      signal_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      TemplateName = "ggZZsignal_TempDataHist_Up_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZsignal_TempHistFunc_Up_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_ggZZ_RawHistPDFUpList[al])
      signal_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_ggZZ_RawHistPDFUpList[al].ProjectionX())
      signal_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_ggZZ_RawHistPDFUpList[al].ProjectionY())
      signal_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      TemplateName = "ggZZsignal_TempDataHist_Down_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZsignal_TempHistFunc_Down_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_ggZZ_RawHistPDFDownList[al])
      signal_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_ggZZ_RawHistPDFDownList[al].ProjectionX())
      signal_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_ggZZ_RawHistPDFDownList[al].ProjectionY())
      signal_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_ggZZ_HistFuncPDFDownList.append(TempHistFunc)

      for al in range(0,anomalousLoops_interf) :
      TemplateName = "ggZZinterf_TempDataHist_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZinterf_TempHistFunc_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_ggZZ_RawHistList[al])
      interf_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_ggZZ_RawHistList[al].ProjectionX())
      interf_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_ggZZ_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_ggZZ_RawHistList[al].ProjectionY())
      interf_ggZZ_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncList.append(TempHistFunc)
      TemplateName = "ggZZinterf_TempDataHist_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZinterf_TempHistFunc_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_ggZZ_RawHistQCDUpList[al])
      interf_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_ggZZ_RawHistQCDUpList[al].ProjectionX())
      interf_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_ggZZ_RawHistQCDUpList[al].ProjectionY())
      interf_ggZZ_DataHistQCDUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncQCDUpList.append(TempHistFunc)
      TemplateName = "ggZZinterf_TempDataHist_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZinterf_TempHistFunc_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_ggZZ_RawHistQCDDownList[al])
      interf_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_ggZZ_RawHistQCDDownList[al].ProjectionX())
      interf_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_ggZZ_RawHistQCDDownList[al].ProjectionY())
      interf_ggZZ_DataHistQCDDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncQCDDownList.append(TempHistFunc)
      TemplateName = "ggZZinterf_TempDataHist_Up_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZinterf_TempHistFunc_Up_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_ggZZ_RawHistPDFUpList[al])
      interf_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_ggZZ_RawHistPDFUpList[al].ProjectionX())
      interf_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_ggZZ_RawHistPDFUpList[al].ProjectionY())
      interf_ggZZ_DataHistPDFUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncPDFUpList.append(TempHistFunc)
      TemplateName = "ggZZinterf_TempDataHist_Down_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "ggZZinterf_TempHistFunc_Down_pdf_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_ggZZ_RawHistPDFDownList[al])
      interf_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_ggZZ_RawHistPDFDownList[al].ProjectionX())
      interf_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_ggZZ_HistFuncPDFDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_ggZZ_RawHistPDFDownList[al].ProjectionY())
      interf_ggZZ_DataHistPDFDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_ggZZ_HistFuncPDFDownList.append(TempHistFunc)


      if self.anomCoupl == 0:
      sigRatesNormList.append(sigRatesWidthNorm)
      interfRatesNormList.append(interfRatesWidthNorm)
      elif self.anomCoupl == 1:
      sigRateNameACNorm = "signalNorm_AC_0_ggZZrate"
      sigRatesACNorm = ROOT.RooFormulaVar(sigRateNameACNorm, "(1.-abs(@0))*@1", ROOT.RooArgList(fai1,sigRatesWidthNorm))
      interfRateNameACNorm = "interfNorm_AC_0_ggZZrate"
      interfRatesACNorm = ROOT.RooFormulaVar(interfRateNameACNorm, "sqrt(1-abs(@0))*@1", ROOT.RooArgList(fai1,interfRatesWidthNorm))
      sigRatesNormList.append(sigRatesACNorm)
      interfRatesNormList.append(interfRatesACNorm)

      sigRateNameACNorm = "signalNorm_AC_1_ggZZrate"
      sigRatesACNorm = ROOT.RooFormulaVar(sigRateNameACNorm, "sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1", ROOT.RooArgList(fai1,sigRatesWidthNorm))
      interfRateNameACNorm = "interfNorm_AC_1_ggZZrate"
      interfRatesACNorm = ROOT.RooFormulaVar(interfRateNameACNorm, "sign(@0)*sqrt(abs(@0))*@1", ROOT.RooArgList(fai1,interfRatesWidthNorm))
      sigRatesNormList.append(sigRatesACNorm)
      interfRatesNormList.append(interfRatesACNorm)

      sigRateNameACNorm = "signalNorm_AC_2_ggZZrate"
      sigRatesACNorm = ROOT.RooFormulaVar(sigRateNameACNorm, "abs(@0)*@1", ROOT.RooArgList(fai1,sigRatesWidthNorm))
      sigRatesNormList.append(sigRatesACNorm)

      ggZZ_funcficients = ROOT.RooArgList()
      for al in range(0,len(sigRatesNormList)) :
      ggZZ_funcficients.add(sigRatesNormList[al])
      for al in range(0,len(interfRatesNormList)) :
      ggZZ_funcficients.add(interfRatesNormList[al])
      ggZZ_funcficients.add(bkgRatesNorm)

      ggZZ_Nominal_histfuncs = ROOT.RooArgList()
      ggZZ_QCDUp_histfuncs = ROOT.RooArgList()
      ggZZ_QCDDown_histfuncs = ROOT.RooArgList()
      ggZZ_PDFUp_histfuncs = ROOT.RooArgList()
      ggZZ_PDFDown_histfuncs = ROOT.RooArgList()
      for al in range(0,anomalousLoops) :
      ggZZ_Nominal_histfuncs.add(signal_ggZZ_HistFuncList[al])
      ggZZ_QCDUp_histfuncs.add(signal_ggZZ_HistFuncQCDUpList[al])
      ggZZ_QCDDown_histfuncs.add(signal_ggZZ_HistFuncQCDDownList[al])
      ggZZ_PDFUp_histfuncs.add(signal_ggZZ_HistFuncPDFUpList[al])
      ggZZ_PDFDown_histfuncs.add(signal_ggZZ_HistFuncPDFDownList[al])
      for al in range(0,anomalousLoops_interf) :
      ggZZ_Nominal_histfuncs.add(interf_ggZZ_HistFuncList[al])
      ggZZ_QCDUp_histfuncs.add(interf_ggZZ_HistFuncQCDUpList[al])
      ggZZ_QCDDown_histfuncs.add(interf_ggZZ_HistFuncQCDDownList[al])
      ggZZ_PDFUp_histfuncs.add(interf_ggZZ_HistFuncPDFUpList[al])
      ggZZ_PDFDown_histfuncs.add(interf_ggZZ_HistFuncPDFDownList[al])
      for al in range(0,1) :
      ggZZ_Nominal_histfuncs.add(bkg_ggZZ_HistFuncList[al])
      ggZZ_QCDUp_histfuncs.add(bkg_ggZZ_HistFuncQCDUpList[al])
      ggZZ_QCDDown_histfuncs.add(bkg_ggZZ_HistFuncQCDDownList[al])
      ggZZ_PDFUp_histfuncs.add(bkg_ggZZ_HistFuncPDFUpList[al])
      ggZZ_PDFDown_histfuncs.add(bkg_ggZZ_HistFuncPDFDownList[al])

      #        print "ggZZ_Nominal_histfuncs"
      #        ggZZ_Nominal_histfuncs.Print("v")

      #        print "ggZZ_funcficients"
      #        ggZZ_funcficients.Print("v")

      mixedList_ggZZ_Nominal = ROOT.RooArgList()
      mixedList_ggZZ_QCDUp = ROOT.RooArgList()
      mixedList_ggZZ_QCDDown = ROOT.RooArgList()
      mixedList_ggZZ_PDFUp = ROOT.RooArgList()
      mixedList_ggZZ_PDFDown = ROOT.RooArgList()
      totalsize = anomalousLoops + anomalousLoops_interf + 1
      genericpdf_ggZZ_formula = "TMath::Max( "
      for al in range(0,totalsize):
      mixedList_ggZZ_Nominal.add(ggZZ_Nominal_histfuncs[al])
      mixedList_ggZZ_Nominal.add(ggZZ_funcficients[al])
      mixedList_ggZZ_QCDUp.add(ggZZ_QCDUp_histfuncs[al])
      mixedList_ggZZ_QCDUp.add(ggZZ_funcficients[al])
      mixedList_ggZZ_QCDDown.add(ggZZ_QCDDown_histfuncs[al])
      mixedList_ggZZ_QCDDown.add(ggZZ_funcficients[al])
      mixedList_ggZZ_PDFUp.add(ggZZ_PDFUp_histfuncs[al])
      mixedList_ggZZ_PDFUp.add(ggZZ_funcficients[al])
      mixedList_ggZZ_PDFDown.add(ggZZ_PDFDown_histfuncs[al])
      mixedList_ggZZ_PDFDown.add(ggZZ_funcficients[al])
      indexA = 2*al
      indexB = 2*al+1
      addstring = "@{0}*@{1}".format(indexA,indexB)
      if al != (totalsize-1): addstring += "+"
      genericpdf_ggZZ_formula += addstring
      genericpdf_ggZZ_formula += " , 0.0000000001 )"
      print "ggZZ formula string: ",genericpdf_ggZZ_formula

      ggZZpdfName = "ggZZ_RooWidth_Nominal2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf2_Nominal = ROOT.RooRealSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_Nominal_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Up2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf2_Up = ROOT.RooRealSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_QCDUp_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Down2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf2_Down = ROOT.RooRealSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_QCDDown_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Up_pdf2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf2_Up_pdf = ROOT.RooRealSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_PDFUp_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Down_pdf2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf2_Down_pdf = ROOT.RooRealSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_PDFDown_histfuncs,ggZZ_funcficients
      )


      ggZZpdfName = "ggZZ_RooWidth_Nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf_Nominal = ROOT.RooRealFlooredSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_Nominal_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf_Up = ROOT.RooRealFlooredSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_QCDUp_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf_Down = ROOT.RooRealFlooredSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_QCDDown_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Up_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf_Up_pdf = ROOT.RooRealFlooredSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_PDFUp_histfuncs,ggZZ_funcficients
      )
      ggZZpdfName = "ggZZ_RooWidth_Down_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      ggZZpdf_Down_pdf = ROOT.RooRealFlooredSumPdf(
      ggZZpdfName, ggZZpdfName,
      ggZZ_PDFDown_histfuncs,ggZZ_funcficients
      )

      #        print "Original value: ",ggZZpdf2_Nominal.getVal()
      #        print "New value: ",ggZZpdf_Nominal.getVal()

      # Shape systematics variables for ggZZ

      CMS_zz4l_APscale_syst = w.factory("QCDscale_ggH[-7,7]")
      CMS_zz4l_pdf_gg_syst = w.factory("pdf_gg[-7,7]")
      morphVarListggZZ = ROOT.RooArgList()
      MorphList_ggZZ = ROOT.RooArgList()
      morphVarListggZZ.add(CMS_zz4l_APscale_syst)
      morphVarListggZZ.add(CMS_zz4l_pdf_gg_syst)
      MorphList_ggZZ.add(ggZZpdf_Nominal)
      MorphList_ggZZ.add(ggZZpdf_Up)
      MorphList_ggZZ.add(ggZZpdf_Down)
      MorphList_ggZZ.add(ggZZpdf_Up_pdf)
      MorphList_ggZZ.add(ggZZpdf_Down_pdf)
      ggZZpdf = ROOT.VerticalInterpPdf("ggzz", "ggzz", MorphList_ggZZ, morphVarListggZZ,1.0) # THIS IS THE ggZZ PDF!!!

      if DEBUG:
      ggzzsignalInt = signal_ggZZ_HistFuncList[0].analyticalIntegral(1000)
      ggzzinterfInt = interf_ggZZ_HistFuncList[0].analyticalIntegral(1000)
      ggzzbkgInt = bkg_ggZZ_HistFuncList[0].analyticalIntegral(1000)
      print "signal_ggZZ_HistFuncList[0]: ",ggzzsignalInt
      print "sigRates_Nominal_AnomCoupl: ",sigRates_Nominal_AnomCoupl.getVal()
      print "interf_ggZZ_HistFuncList[0]: ",ggzzinterfInt
      print "interfRates_Nominal_AnomCoupl: ",interfRates_Nominal_AnomCoupl.getVal()
      print "bkg_ggZZ_HistFuncList[0]: ",ggzzbkgInt
      print "bkgRates_Nominal: ",bkgRates_Nominal.getVal()
      print "ggZZNominal_norm (direct): ",ggZZNominal_norm.getVal()
      print "ggZZNominal_norm (calculated): {0:.12f}".format(sigRates_Nominal_AnomCoupl.getVal()+interfRates_Nominal_AnomCoupl.getVal()+bkgRates_Nominal.getVal())
      print "ggZZpdf_Nominal (calculated): {0:.12f}".format(ggzzsignalInt+ggzzinterfInt+ggzzbkgInt)
      print "ggZZpdf_Nominal PDF (direct): ",ggZZpdf_Nominal.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
      print "ggZZpdf_Nominal RooRealSum PDF (direct): ",ggZZpdf2_Nominal.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
      print "ggZZ Vertical interpolator PDF: ",ggZZpdf.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
      ggZZpdf.Print("v")
      ggZZpdf_Nominal.Print("v")
      signal_ggZZ_HistFuncList[0].Print("v")
      interf_ggZZ_HistFuncList[0].Print("v")
      bkg_ggZZ_HistFuncList[0].Print("v")
      ggZZNominal_norm.Print("v")
      sigRates_Nominal_AnomCoupl.Print("v")
      interfRates_Nominal_AnomCoupl.Print("v")
      bkgRates_Nominal.Print("v")


      asympowname = "kappalow_ggZZ_QCD_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      kappalow = ROOT.RooFormulaVar(
      asympowname, "@0/@1", ROOT.RooArgList(ggZZQCDDown_norm, ggZZNominal_norm)
      )
      asympowname = "kappahigh_ggZZ_QCD_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      kappahigh = ROOT.RooFormulaVar(
      asympowname, "@0/@1", ROOT.RooArgList(ggZZQCDUp_norm, ggZZNominal_norm)
      )
      asympowname = "kappalow_ggZZ_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      kappalow_pdf = ROOT.RooFormulaVar(
      asympowname, "@0/@1", ROOT.RooArgList(ggZZPDFDown_norm, ggZZNominal_norm)
      )
      asympowname = "kappahigh_ggZZ_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      kappahigh_pdf = ROOT.RooFormulaVar(
      asympowname, "@0/@1", ROOT.RooArgList(ggZZPDFUp_norm, ggZZNominal_norm)
      )

      MorphNormList_ggZZ = ROOT.RooArgList()
      MorphNormList_ggZZ.add(ggZZNominal_norm)
      MorphNormList_ggZZ.add(ggZZQCDUp_norm)
      MorphNormList_ggZZ.add(ggZZQCDDown_norm)
      MorphNormList_ggZZ.add(ggZZPDFUp_norm)
      MorphNormList_ggZZ.add(ggZZPDFDown_norm)


      #        asympowname = "Asympow_ggZZ_QCD_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      #        thetaSyst_ggZZ = AsymPow(
      #            asympowname, asympowname, kappalow, kappahigh, CMS_zz4l_APscale_syst
      #        )
      #        asympowname = "Asympow_ggZZ_pdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      #        thetaSyst_ggZZ_pdf = AsymPow(
      #            asympowname, asympowname, kappalow_pdf, kappahigh_pdf, CMS_zz4l_pdf_gg_syst
      #        )
      asympowname = "Asymquad_ggZZ_QCDPDF_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_ggZZ_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_ggZZ, morphVarListggZZ, 1.0) # THIS IS THE ggZZ PDF!!!



      # -------------------------- VBF offshell SHAPES ---------------------------------- ##
      VBFsigRatesNormList=[]
      VBFinterfRatesNormList=[]

      bkg_VBF_RawHistList=[VBF_T_1]
      bkg_VBF_DataHistList=[]
      bkg_VBF_HistFuncList=[]

      bkg_VBF_RawHistUpList=[VBF_T_1_Up]
      bkg_VBF_DataHistUpList=[]
      bkg_VBF_HistFuncUpList=[]

      bkg_VBF_RawHistDownList=[VBF_T_1_Down]
      bkg_VBF_DataHistDownList=[]
      bkg_VBF_HistFuncDownList=[]


      signal_VBF_RawHistList=[VBF_T_2,VBF_T_mZZ2_1_2,VBF_T_mZZ2_2_2,VBF_T_mZZ2_3_2,VBF_T_mZZ2_4_2]
      signal_VBF_DataHistList=[]
      signal_VBF_HistFuncList=[]

      signal_VBF_RawHistUpList=[VBF_T_2_Up,VBF_T_mZZ2_1_2_Up,VBF_T_mZZ2_2_2_Up,VBF_T_mZZ2_3_2_Up,VBF_T_mZZ2_4_2_Up]
      signal_VBF_DataHistUpList=[]
      signal_VBF_HistFuncUpList=[]

      signal_VBF_RawHistDownList=[VBF_T_2_Down,VBF_T_mZZ2_1_2_Down,VBF_T_mZZ2_2_2_Down,VBF_T_mZZ2_3_2_Down,VBF_T_mZZ2_4_2_Down]
      signal_VBF_DataHistDownList=[]
      signal_VBF_HistFuncDownList=[]


      interf_VBF_RawHistList=[VBF_T_4,VBF_T_mZZ2_1_4,VBF_T_mZZ2_2_4]
      interf_VBF_DataHistList=[]
      interf_VBF_HistFuncList=[]

      interf_VBF_RawHistUpList=[VBF_T_4_Up,VBF_T_mZZ2_1_4_Up,VBF_T_mZZ2_2_4_Up]
      interf_VBF_DataHistUpList=[]
      interf_VBF_HistFuncUpList=[]

      interf_VBF_RawHistDownList=[VBF_T_4_Down,VBF_T_mZZ2_1_4_Down,VBF_T_mZZ2_2_4_Down]
      interf_VBF_DataHistDownList=[]
      interf_VBF_HistFuncDownList=[]


      VBFbkgRateNameNorm = "bkgNorm_VBFrate"
      VBFbkgRatesNorm = ROOT.RooRealVar(VBFbkgRateNameNorm, VBFbkgRateNameNorm, 1.0)
      VBFbkgRatesNorm.setConstant(True)

      VBFsigRateNameWidthNorm = "signalWidthNorm_VBFrate"
      VBFinterfRateNameWidthNorm = "interfWidthNorm_VBFrate"

      VBFsigRatesWidthNorm = ROOT.RooFormulaVar(VBFsigRateNameWidthNorm, "@0*@1*@2", ROOT.RooArgList(x, mu, muV))
      VBFinterfRatesWidthNorm = ROOT.RooFormulaVar(VBFinterfRateNameWidthNorm, "sqrt(@0*@1*@2)", ROOT.RooArgList(x, mu, muV))


      # VBF Bkg histfunc construction
      TemplateName = "VBFbkg_TempDataHist_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "VBFbkg_TempHistFunc_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_VBF_RawHistList[0])
      bkg_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_VBF_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_VBF_RawHistList[0].ProjectionX())
      bkg_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_VBF_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_VBF_RawHistList[0].ProjectionY())
      bkg_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_VBF_HistFuncList.append(TempHistFunc)
      TemplateName = "VBFbkg_TempDataHist_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "VBFbkg_TempHistFunc_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_VBF_RawHistDownList[0])
      bkg_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_VBF_HistFuncDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_VBF_RawHistDownList[0].ProjectionX())
      bkg_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_VBF_HistFuncDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_VBF_RawHistDownList[0].ProjectionY())
      bkg_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_VBF_HistFuncDownList.append(TempHistFunc)
      TemplateName = "VBFbkg_TempDataHist_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      PdfName = "VBFbkg_TempHistFunc_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), bkg_VBF_RawHistUpList[0])
      bkg_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      bkg_VBF_HistFuncUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), bkg_VBF_RawHistUpList[0].ProjectionX())
      bkg_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      bkg_VBF_HistFuncUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), bkg_VBF_RawHistUpList[0].ProjectionY())
      bkg_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      bkg_VBF_HistFuncUpList.append(TempHistFunc)


      anomalousLoops = 1
      anomalousLoops_interf = 1
      if self.anomCoupl == 1:
      anomalousLoops = 5 # For VBF
      anomalousLoops_interf = 3

      for al in range(0,anomalousLoops) :
      TemplateName = "VBFsignal_TempDataHist_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "VBFsignal_TempHistFunc_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_VBF_RawHistList[al])
      signal_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_VBF_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_VBF_RawHistList[al].ProjectionX())
      signal_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_VBF_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_VBF_RawHistList[al].ProjectionY())
      signal_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_VBF_HistFuncList.append(TempHistFunc)
      TemplateName = "VBFsignal_TempDataHist_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "VBFsignal_TempHistFunc_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_VBF_RawHistUpList[al])
      signal_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_VBF_HistFuncUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_VBF_RawHistUpList[al].ProjectionX())
      signal_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_VBF_HistFuncUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_VBF_RawHistUpList[al].ProjectionY())
      signal_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_VBF_HistFuncUpList.append(TempHistFunc)
      TemplateName = "VBFsignal_TempDataHist_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "VBFsignal_TempHistFunc_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), signal_VBF_RawHistDownList[al])
      signal_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      signal_VBF_HistFuncDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), signal_VBF_RawHistDownList[al].ProjectionX())
      signal_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      signal_VBF_HistFuncDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), signal_VBF_RawHistDownList[al].ProjectionY())
      signal_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      signal_VBF_HistFuncDownList.append(TempHistFunc)

      for al in range(0,anomalousLoops_interf) :
      TemplateName = "VBFinterf_TempDataHist_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "VBFinterf_TempHistFunc_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_VBF_RawHistList[al])
      interf_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_VBF_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_VBF_RawHistList[al].ProjectionX())
      interf_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_VBF_HistFuncList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_VBF_RawHistList[al].ProjectionY())
      interf_VBF_DataHistList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_VBF_HistFuncList.append(TempHistFunc)
      TemplateName = "VBFinterf_TempDataHist_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "VBFinterf_TempHistFunc_Up_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_VBF_RawHistUpList[al])
      interf_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_VBF_HistFuncUpList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_VBF_RawHistUpList[al].ProjectionX())
      interf_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_VBF_HistFuncUpList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_VBF_RawHistUpList[al].ProjectionY())
      interf_VBF_DataHistUpList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_VBF_HistFuncUpList.append(TempHistFunc)
      TemplateName = "VBFinterf_TempDataHist_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      PdfName = "VBFinterf_TempHistFunc_Down_AC{3:.0f}_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet, al)
      if self.dimensions > 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), interf_VBF_RawHistDownList[al])
      interf_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), TempDataHist)
      interf_VBF_HistFuncDownList.append(TempHistFunc)
      elif self.dimensions == 1:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass), interf_VBF_RawHistDownList[al].ProjectionX())
      interf_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass), TempDataHist)
      interf_VBF_HistFuncDownList.append(TempHistFunc)
      elif self.dimensions == 0:
      TempDataHist = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthKD), interf_VBF_RawHistDownList[al].ProjectionY())
      interf_VBF_DataHistDownList.append(TempDataHist)
      TempHistFunc = ROOT.RooHistFunc(PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), TempDataHist)
      interf_VBF_HistFuncDownList.append(TempHistFunc)


      if self.anomCoupl == 0:
      VBFsigRatesNormList.append(VBFsigRatesWidthNorm)
      VBFinterfRatesNormList.append(VBFinterfRatesWidthNorm)
      elif self.anomCoupl == 1:
      VBFsigRateNameACNorm = "signalNorm_AC_0_VBFrate"
      VBFsigRatesACNorm = ROOT.RooFormulaVar(VBFsigRateNameACNorm, "pow((1.-abs(@0)),2)*@1", ROOT.RooArgList(fai1,VBFsigRatesWidthNorm))
      VBFinterfRateNameACNorm = "interfNorm_AC_0_VBFrate"
      VBFinterfRatesACNorm = ROOT.RooFormulaVar(VBFinterfRateNameACNorm, "(1-abs(@0))*@1", ROOT.RooArgList(fai1,VBFinterfRatesWidthNorm))
      VBFsigRatesNormList.append(VBFsigRatesACNorm)
      VBFinterfRatesNormList.append(VBFinterfRatesACNorm)

      VBFsigRateNameACNorm = "signalNorm_AC_1_VBFrate"
      VBFsigRatesACNorm = ROOT.RooFormulaVar(VBFsigRateNameACNorm, "sign(@0)*sqrt(abs(@0)*pow(sqrt(1-abs(@0)),3))*@1", ROOT.RooArgList(fai1,VBFsigRatesWidthNorm))
      VBFinterfRateNameACNorm = "interfNorm_AC_1_VBFrate"
      VBFinterfRatesACNorm = ROOT.RooFormulaVar(VBFinterfRateNameACNorm, "sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*@1", ROOT.RooArgList(fai1,VBFinterfRatesWidthNorm))
      VBFsigRatesNormList.append(VBFsigRatesACNorm)
      VBFinterfRatesNormList.append(VBFinterfRatesACNorm)

      VBFsigRateNameACNorm = "signalNorm_AC_2_VBFrate"
      VBFsigRatesACNorm = ROOT.RooFormulaVar(VBFsigRateNameACNorm, "abs(@0)*(1-abs(@0))*@1", ROOT.RooArgList(fai1,VBFsigRatesWidthNorm))
      VBFinterfRateNameACNorm = "interfNorm_AC_2_VBFrate"
      VBFinterfRatesACNorm = ROOT.RooFormulaVar(VBFinterfRateNameACNorm, "abs(@0)*@1", ROOT.RooArgList(fai1,VBFinterfRatesWidthNorm))
      VBFsigRatesNormList.append(VBFsigRatesACNorm)
      VBFinterfRatesNormList.append(VBFinterfRatesACNorm)

      VBFsigRateNameACNorm = "signalNorm_AC_3_VBFrate"
      VBFsigRatesACNorm = ROOT.RooFormulaVar(VBFsigRateNameACNorm, "sign(@0)*pow(sqrt(abs(@0)),3)*sqrt(1-abs(@0))*@1", ROOT.RooArgList(fai1,VBFsigRatesWidthNorm))
      VBFsigRatesNormList.append(VBFsigRatesACNorm)

      VBFsigRateNameACNorm = "signalNorm_AC_4_VBFrate"
      VBFsigRatesACNorm = ROOT.RooFormulaVar(VBFsigRateNameACNorm, "pow(@0,2)*@1", ROOT.RooArgList(fai1,VBFsigRatesWidthNorm))
      VBFsigRatesNormList.append(VBFsigRatesACNorm)

      VBF_funcficients = ROOT.RooArgList()
      for al in range(0,len(VBFsigRatesNormList)) :
      VBF_funcficients.add(VBFsigRatesNormList[al])
      for al in range(0,len(VBFinterfRatesNormList)) :
      VBF_funcficients.add(VBFinterfRatesNormList[al])
      VBF_funcficients.add(VBFbkgRatesNorm)


      VBF_Nominal_histfuncs = ROOT.RooArgList()
      VBF_Up_histfuncs = ROOT.RooArgList()
      VBF_Down_histfuncs = ROOT.RooArgList()
      for al in range(0,anomalousLoops) :
      VBF_Nominal_histfuncs.add(signal_VBF_HistFuncList[al])
      VBF_Up_histfuncs.add(signal_VBF_HistFuncUpList[al])
      VBF_Down_histfuncs.add(signal_VBF_HistFuncDownList[al])
      for al in range(0,anomalousLoops_interf) :
      VBF_Nominal_histfuncs.add(interf_VBF_HistFuncList[al])
      VBF_Up_histfuncs.add(interf_VBF_HistFuncUpList[al])
      VBF_Down_histfuncs.add(interf_VBF_HistFuncDownList[al])
      for al in range(0,1) :
      VBF_Nominal_histfuncs.add(bkg_VBF_HistFuncList[al])
      VBF_Up_histfuncs.add(bkg_VBF_HistFuncUpList[al])
      VBF_Down_histfuncs.add(bkg_VBF_HistFuncDownList[al])

      VBFpdfName = "VBF_RooWidth_Nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFpdf_Nominal = ROOT.RooRealFlooredSumPdf(
      VBFpdfName, VBFpdfName,
      VBF_Nominal_histfuncs,VBF_funcficients
      )
      VBFpdfName = "VBF_RooWidth_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFpdf_Up = ROOT.RooRealFlooredSumPdf(
      VBFpdfName, VBFpdfName,
      VBF_Up_histfuncs,VBF_funcficients
      )
      VBFpdfName = "VBF_RooWidth_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      VBFpdf_Down = ROOT.RooRealFlooredSumPdf(
      VBFpdfName, VBFpdfName,
      VBF_Down_histfuncs,VBF_funcficients
      )


      CMS_zz4l_VBFscale_syst = ROOT.RooRealVar("CMS_zz4l_VBFscale_syst", "CMS_zz4l_VBFscale_syst", 0.0, -1, 1)
      CMS_zz4l_VBFscale_syst.setConstant(True) #NOTE THIS IS KILLING THIS SYSTEMATIC
      morphVarListVBF = ROOT.RooArgList()
      morphVarListVBF.add(CMS_zz4l_VBFscale_syst)
      MorphList_VBF = ROOT.RooArgList()
      MorphList_VBF.add(VBFpdf_Nominal)
      MorphList_VBF.add(VBFpdf_Up)
      MorphList_VBF.add(VBFpdf_Down)

      VBFpdf = ROOT.VerticalInterpPdf("VBFpdf", "VBFpdf", MorphList_VBF, morphVarListVBF,1.0) # THIS IS THE VBF PDF!!!

      if DEBUG:
      vbfsignalInt = signal_VBF_HistFuncList[0].analyticalIntegral(1000)
      vbfinterfInt = interf_VBF_HistFuncList[0].analyticalIntegral(1000)
      vbfbkgInt = bkg_VBF_HistFuncList[0].analyticalIntegral(1000)
      print "signal_VBF_HistFuncList[0]: ",vbfsignalInt
      print "VBFsigRates_Nominal_AnomCoupl: ",VBFsigRates_Nominal_AnomCoupl.getVal()
      print "interf_VBF_HistFuncList[0]: ",vbfinterfInt
      print "VBFinterfRates_Nominal_AnomCoupl: ",VBFinterfRates_Nominal_AnomCoupl.getVal()
      print "bkg_VBF_HistFuncList[0]: ",vbfbkgInt
      print "VBFbkgRates_Nominal: ",VBFbkgRates_Nominal.getVal()
      print "VBFNominal_norm (direct): ",VBFNominal_norm.getVal()
      print "VBFNominal_norm (calculated): {0:.12f}".format(VBFsigRates_Nominal_AnomCoupl.getVal()+VBFinterfRates_Nominal_AnomCoupl.getVal()+VBFbkgRates_Nominal.getVal())
      print "VBFpdf_Nominal PDF (calculated): {0:.12f}".format(vbfsignalInt+vbfinterfInt+vbfbkgInt)
      print "VBFpdf_Nominal PDF (direct): ",VBFpdf_Nominal.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
      print "VBF Vertical interpolator PDF: ",VBFpdf.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
      VBFpdf.Print("v")
      VBFpdf_Nominal.Print("v")
      signal_VBF_HistFuncList[0].Print("v")
      interf_VBF_HistFuncList[0].Print("v")
      bkg_VBF_HistFuncList[0].Print("v")
      VBFNominal_norm.Print("v")
      VBFsigRates_Nominal_AnomCoupl.Print("v")
      VBFinterfRates_Nominal_AnomCoupl.Print("v")
      VBFbkgRates_Nominal.Print("v")

      MorphNormList_VBF = ROOT.RooArgList()
      MorphNormList_VBF.add(VBFNominal_norm)
      MorphNormList_VBF.add(VBFUp_norm)
      MorphNormList_VBF.add(VBFDown_norm)

      #        asympowname = "kappalow_VBF_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      #        kappalowVBF = ROOT.RooFormulaVar(asympowname, "@0/@1", ROOT.RooArgList(VBFDown_norm, VBFNominal_norm))
      #        asympowname = "kappahigh_VBF_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      #        kappahighVBF = ROOT.RooFormulaVar(asympowname, "@0/@1", ROOT.RooArgList(VBFUp_norm, VBFNominal_norm))
      #        asympowname = "Asympow_VBF_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      #        thetaSyst_VBF = AsymPow(asympowname, asympowname, kappalowVBF, kappahighVBF, CMS_zz4l_VBFscale_syst)
      asympowname = "Asymquad_VBF_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_VBF_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_VBF, morphVarListVBF, 1.0) # THIS IS THE ggZZ PDF!!!


      #------------------ SIGNAL PDF NORM VARIABLES -------------------------

      ggZZpdfNormName = "ggZZ_RooWidth_{0:.0f}_{1:.0f}_{2:.0f}_norm".format(self.channel, self.sqrts, useDjet)
      #        ggZZpdf_norm = ROOT.RooFormulaVar(ggZZpdfNormName, "@0*@1*@2*@3", ROOT.RooArgList(ggZZNominal_norm, self.LUMI, thetaSyst_ggZZ, thetaSyst_ggZZ_pdf))
      ggZZpdf_norm = ROOT.RooFormulaVar(ggZZpdfNormName, "TMath::Max(@0*@1,1e-10)", ROOT.RooArgList(thetaSyst_ggZZ_norm, self.LUMI))
      ggZZpdf_norm.SetNameTitle("ggzz_norm", "ggzz_norm")
      VBFpdfNormName = "VBF_RooWidth_{0:.0f}_{1:.0f}_{2:.0f}_norm".format(self.channel, self.sqrts, useDjet)
      #        VBFpdf_norm = ROOT.RooFormulaVar(VBFpdfNormName, "@0*@1*@2", ROOT.RooArgList(VBFNominal_norm, self.LUMI, thetaSyst_VBF))
      VBFpdf_norm = ROOT.RooFormulaVar(VBFpdfNormName, "TMath::Max(@0*@1,1e-10)", ROOT.RooArgList(thetaSyst_VBF_norm, self.LUMI))
      VBFpdf_norm.SetNameTitle("vbf_offshell_norm", "vbf_offshell_norm")


      # -------------------------- OTHER BACKGROUND SHAPES ---------------------------------- ##

      # rates per lumi for scaling
      bkgRate_qqzz = self.theInputs['qqZZ_rate'] / self.theInputs['qqZZ_lumi']  # *1.8
      bkgRate_zjets = self.theInputs['zjets_rate'] / self.theInputs['zjets_lumi']
      bkgRate_qqzz_Shape = bkgRate_qqzz * self.lumi
      bkgRate_zjets_Shape = bkgRate_zjets * self.lumi

      # qqZZ contribution
      name = "CMS_widthqqzzbkg_a0_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a0 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a0", self.theInputs['qqZZshape_a0'])
      name = "CMS_widthqqzzbkg_a1_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a1 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a1", self.theInputs['qqZZshape_a1'])
      name = "CMS_widthqqzzbkg_a2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a2 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a2", self.theInputs['qqZZshape_a2'])
      name = "CMS_widthqqzzbkg_a3_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a3 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a3", self.theInputs['qqZZshape_a3'])
      name = "CMS_widthqqzzbkg_a4_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a4 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a4", self.theInputs['qqZZshape_a4'])
      name = "CMS_widthqqzzbkg_a5_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a5 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a5", self.theInputs['qqZZshape_a5'])
      name = "CMS_widthqqzzbkg_a6_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a6 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a6", self.theInputs['qqZZshape_a6'])
      name = "CMS_widthqqzzbkg_a7_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a7 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a7", self.theInputs['qqZZshape_a7'])
      name = "CMS_widthqqzzbkg_a8_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a8 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a8", self.theInputs['qqZZshape_a8'])
      name = "CMS_widthqqzzbkg_a9_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a9 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a9", self.theInputs['qqZZshape_a9'])
      name = "CMS_widthqqzzbkg_a10_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a10 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a10", self.theInputs['qqZZshape_a10'])
      name = "CMS_widthqqzzbkg_a11_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a11 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a11", self.theInputs['qqZZshape_a11'])
      name = "CMS_widthqqzzbkg_a12_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a12 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a12", self.theInputs['qqZZshape_a12'])
      name = "CMS_widthqqzzbkg_a13_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      CMS_qqzzbkg_a13 = ROOT.RooRealVar(name, "CMS_widthqqzzbkg_a13", self.theInputs['qqZZshape_a13'])

      CMS_qqzzbkg_a0.setConstant(True)
      CMS_qqzzbkg_a1.setConstant(True)
      CMS_qqzzbkg_a2.setConstant(True)
      CMS_qqzzbkg_a3.setConstant(True)
      CMS_qqzzbkg_a4.setConstant(True)
      CMS_qqzzbkg_a5.setConstant(True)
      CMS_qqzzbkg_a6.setConstant(True)
      CMS_qqzzbkg_a7.setConstant(True)
      CMS_qqzzbkg_a8.setConstant(True)
      CMS_qqzzbkg_a9.setConstant(True)
      CMS_qqzzbkg_a10.setConstant(True)
      CMS_qqzzbkg_a11.setConstant(True)
      CMS_qqzzbkg_a12.setConstant(True)
      CMS_qqzzbkg_a13.setConstant(True)

      # TO BE CLEANED UP ->this part should be moved in inputs
      CMS_qqzzbkg_p0 = ROOT.RooRealVar("CMS_widthqqzzbkg_p0", "CMS_widthqqzzbkg_p0", 1.04012)
      CMS_qqzzbkg_p1 = ROOT.RooRealVar("CMS_widthqqzzbkg_p1", "CMS_widthqqzzbkg_p1", -0.000125088)
      CMS_qqzzbkg_p2 = ROOT.RooRealVar("CMS_widthqqzzbkg_p2", "CMS_widthqqzzbkg_p2", 2.39404e-07)
      CMS_qqzzbkg_p3 = ROOT.RooRealVar("CMS_widthqqzzbkg_p3", "CMS_widthqqzzbkg_p3", 1. - 0.034)
      CMS_qqzzbkg_p4 = ROOT.RooRealVar("CMS_widthqqzzbkg_p4", "CMS_widthqqzzbkg_p4", 1. + 0.027)
      CMS_qqzzbkg_p0.setConstant(True)
      CMS_qqzzbkg_p1.setConstant(True)
      CMS_qqzzbkg_p2.setConstant(True)
      CMS_qqzzbkg_p3.setConstant(True)
      CMS_qqzzbkg_p4.setConstant(True)

      # TO BE CLEANED UP ->this part should be moved in inputs
      CMS_qqzzbkg_EWK_p0 = ROOT.RooRealVar("CMS_widthqqzzbkg_EWK_p0", "CMS_widthqqzzbkg_EWK_p0", 0.953385)
      CMS_qqzzbkg_EWK_p1 = ROOT.RooRealVar("CMS_widthqqzzbkg_EWK_p1", "CMS_widthqqzzbkg_EWK_p1", 0.000412406)
      CMS_qqzzbkg_EWK_p2 = ROOT.RooRealVar("CMS_widthqqzzbkg_EWK_p2", "CMS_widthqqzzbkg_EWK_p2", -5.45148e-07)
      CMS_qqzzbkg_EWK_p3 = ROOT.RooRealVar("CMS_widthqqzzbkg_EWK_p3", "CMS_widthqqzzbkg_EWK_p3", 2.63944e-10)
      CMS_qqzzbkg_EWK_p4 = ROOT.RooRealVar("CMS_widthqqzzbkg_EWK_p4", "CMS_widthqqzzbkg_EWK_p4", 1. - 0.029)
      CMS_qqzzbkg_EWK_p5 = ROOT.RooRealVar("CMS_widthqqzzbkg_EWK_p5", "CMS_widthqqzzbkg_EWK_p5", 1. + 0.029)
      CMS_qqzzbkg_EWK_p0.setConstant(True)
      CMS_qqzzbkg_EWK_p1.setConstant(True)
      CMS_qqzzbkg_EWK_p2.setConstant(True)
      CMS_qqzzbkg_EWK_p3.setConstant(True)
      CMS_qqzzbkg_EWK_p4.setConstant(True)
      CMS_qqzzbkg_EWK_p5.setConstant(True)

      bkg_qqzz_mass_temp = ROOT.RooqqZZPdf_v2(
      "bkg_qqzz_widthmass_temp_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet),
      "bkg_qqzz_widthmass_temp", CMS_zz4l_widthMass,
      CMS_qqzzbkg_a0, CMS_qqzzbkg_a1, CMS_qqzzbkg_a2, CMS_qqzzbkg_a3,
      CMS_qqzzbkg_a4, CMS_qqzzbkg_a5, CMS_qqzzbkg_a6, CMS_qqzzbkg_a7,
      CMS_qqzzbkg_a8, CMS_qqzzbkg_a9, CMS_qqzzbkg_a10, CMS_qqzzbkg_a11,
      CMS_qqzzbkg_a12, CMS_qqzzbkg_a13
      )

      name = "qqzz_djetcut_p0_{0:.0f}_{1:.0f}".format(self.sqrts, useDjet)
      qqzz_djetcut_p0 = ROOT.RooRealVar(name, name, 1,-1e5,1e5)
      name = "qqzz_djetcut_p1_{0:.0f}_{1:.0f}".format(self.sqrts, useDjet)
      qqzz_djetcut_p1 = ROOT.RooRealVar(name, name, 0,-1e5,1e5)
      name = "qqzz_djetcut_p2_{0:.0f}_{1:.0f}".format(self.sqrts, useDjet)
      qqzz_djetcut_p2 = ROOT.RooRealVar(name, name, 0,-1e5,1e5)
      name = "qqzz_djetcut_p3_{0:.0f}_{1:.0f}".format(self.sqrts, useDjet)
      qqzz_djetcut_p3 = ROOT.RooRealVar(name, name, 1,-1e5,1e5)

      if useDjet != 0:
      # code for analytic form of Djet > 0.5 cut, useDjet == 2; useDjet == 1 is 1-(useDjet == 2)
      if self.sqrts == 7 :
      qqzz_djetcut_p0.setVal(3.98058e-03)
      qqzz_djetcut_p1.setVal(1.04364e-03)
      qqzz_djetcut_p2.setVal(2.00000e+02)
      qqzz_djetcut_p3.setVal(1.79004e+02)
      if self.sqrts == 8 :
      qqzz_djetcut_p0.setVal(5.64595e-03)
      qqzz_djetcut_p1.setVal(8.67370e-04)
      qqzz_djetcut_p2.setVal(2.73223e+02)
      qqzz_djetcut_p3.setVal(4.00070e+01)

      qqzz_djetcut_p0.setConstant(True)
      qqzz_djetcut_p1.setConstant(True)
      qqzz_djetcut_p2.setConstant(True)
      qqzz_djetcut_p3.setConstant(True)
      djetcutarglist = ROOT.RooArgList(
      CMS_zz4l_widthMass,
      qqzz_djetcut_p0, qqzz_djetcut_p1, qqzz_djetcut_p2, qqzz_djetcut_p3
      )
      bkg_qqzz_mass_Djet_shape_Name = "bkg_qqzz_widthmass_Djet_shape_{0:.0f}_{1:.0f}".format(self.sqrts, useDjet)
      bkg_qqzz_mass_Djet_shape_FormulaCore = "(@1+@2*TMath::Erf((@0-@3)/@4))"
      bkg_qqzz_mass_Djet_shape_Formula = "1"
      if useDjet!=0:
      bkg_qqzz_mass_Djet_shape_Formula = bkg_qqzz_mass_Djet_shape_FormulaCore
      bkg_qqzz_mass_Djet_shape = ROOT.RooFormulaVar(
      bkg_qqzz_mass_Djet_shape_Name, bkg_qqzz_mass_Djet_shape_Name,
      bkg_qqzz_mass_Djet_shape_Formula,
      djetcutarglist
      )

      qqZZ_Scale_Syst = w.factory("QCDscale_VV[-7,7]")
      qqZZ_EWK_Syst = w.factory("EWKcorr_VV[-7,7]")

      qqzzarglist = ROOT.RooArgList(
      qqZZ_EWK_Syst,
      CMS_qqzzbkg_EWK_p0,
      CMS_qqzzbkg_EWK_p1,
      CMS_qqzzbkg_EWK_p2,
      CMS_qqzzbkg_EWK_p3,

      CMS_zz4l_widthMass,

      qqZZ_Scale_Syst,
      CMS_qqzzbkg_p0,
      CMS_qqzzbkg_p1
      )
      qqzzarglist.add(CMS_qqzzbkg_p2)
      #        bkg_qqzz_mass_shape = ROOT.RooFormulaVar(
      bkg_qqzz_mass_shape = ROOT.RooGenericPdf(
      "bkg_qqzz_widthmass_shape_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_qqzz_mass_shape",
      "TMath::Max( 1 + @0*( @1-1 +@2*@5 +@3*@5*@5 +@4*@5*@5*@5 ) + @6*( @7-1 +@8*@5 +@9*@5*@5 ) , 1e-15 )",
      qqzzarglist
      )
      if useDjet != 0:
      qqzzarglist.add(bkg_qqzz_mass_Djet_shape)
      qqzzarglist.add(thetaSyst_djet_qqZZ_norm)
      if useDjet == 2:
      #                bkg_qqzz_mass_shape = ROOT.RooFormulaVar(
      bkg_qqzz_mass_shape = ROOT.RooGenericPdf(
      "bkg_qqzz_widthmass_shape_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_qqzz_mass_shape",
      "TMath::Max( ( TMath::Max(@11,0) + @0*( @1-1 +@2*@5 +@3*@5*@5 +@4*@5*@5*@5 ) + @6*( @7-1 +@8*@5 +@9*@5*@5 ) )*@10, 1e-15 )",
      qqzzarglist
      )
      elif useDjet == 1:
      #                bkg_qqzz_mass_shape = ROOT.RooFormulaVar(
      bkg_qqzz_mass_shape = ROOT.RooGenericPdf(
      "bkg_qqzz_widthmass_shape_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_qqzz_mass_shape",
      "TMath::Max( ( 1 + @0*( @1-1 +@2*@5 +@3*@5*@5 +@4*@5*@5*@5 ) + @6*( @7-1 +@8*@5 +@9*@5*@5 ) )*(1.-@10) + @10*(1-TMath::Max(@11,0)), 1e-15 )",
      qqzzarglist
      )

      bkg_qqzz_mass = ROOT.RooProdPdf(
      "bkg_qqzz_widthmass_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_qqzz_mass",
      ROOT.RooArgList(
      bkg_qqzz_mass_temp,
      bkg_qqzz_mass_shape
      )
      )
      #        bkg_qqzz_mass = ROOT.RooGenericPdf(
      #            "bkg_qqzz_widthmass_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_qqzz_mass","@0*@1",
      #            ROOT.RooArgList(
      #                bkg_qqzz_mass_temp,
      #                bkg_qqzz_mass_shape
      #            )
      #        )

      bkg_qqzz_mass.forceNumInt(True)
      if DEBUG:
      bkg_qqzz_mass_temp.Print("v")
      bkg_qqzz_mass_shape.Print("v")
      bkg_qqzz_mass.Print("v")

      bkg_qqzz_Nominal_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      bkg_qqzz_Nominal_Norm_Name = "bkg_qqzz_widthmass_NominalNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_Nominal_Norm = ROOT.RooRealVar(bkg_qqzz_Nominal_Norm_Name, bkg_qqzz_Nominal_Norm_Name, 1.)

      bkg_qqzz_DjetUp_NormVal = bkg_qqzz_Nominal_NormVal
      bkg_qqzz_DjetDown_NormVal = bkg_qqzz_Nominal_NormVal
      if useDjet!=0:
      CMS_zz4l_qqzz_djet_syst.setVal(1)
      bkg_qqzz_DjetUp_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      CMS_zz4l_qqzz_djet_syst.setVal(-1)
      bkg_qqzz_DjetDown_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      CMS_zz4l_qqzz_djet_syst.setVal(0)
      bkg_qqzz_DjetUp_NormVal = bkg_qqzz_DjetUp_NormVal/bkg_qqzz_Nominal_NormVal
      bkg_qqzz_DjetDown_NormVal = bkg_qqzz_DjetDown_NormVal/bkg_qqzz_Nominal_NormVal
      bkg_qqzz_DjetUp_Norm_Name = "bkg_qqzz_widthmass_DjetUpNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_DjetUp_Norm = ROOT.RooRealVar(bkg_qqzz_DjetUp_Norm_Name, bkg_qqzz_DjetUp_Norm_Name, bkg_qqzz_DjetUp_NormVal)
      bkg_qqzz_DjetDown_Norm_Name = "bkg_qqzz_widthmass_DjetDownNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_DjetDown_Norm = ROOT.RooRealVar(bkg_qqzz_DjetDown_Norm_Name, bkg_qqzz_DjetDown_Norm_Name, bkg_qqzz_DjetDown_NormVal)

      qqZZ_Scale_Syst.setVal(1)
      bkg_qqzz_QCDScaleUp_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      qqZZ_Scale_Syst.setVal(-1)
      bkg_qqzz_QCDScaleDown_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      qqZZ_Scale_Syst.setVal(0)
      bkg_qqzz_QCDScaleUp_NormVal = bkg_qqzz_QCDScaleUp_NormVal/bkg_qqzz_Nominal_NormVal
      bkg_qqzz_QCDScaleDown_NormVal = bkg_qqzz_QCDScaleDown_NormVal/bkg_qqzz_Nominal_NormVal
      bkg_qqzz_QCDScaleUp_Norm_Name = "bkg_qqzz_widthmass_QCDScaleUpNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_QCDScaleUp_Norm = ROOT.RooRealVar(bkg_qqzz_QCDScaleUp_Norm_Name, bkg_qqzz_QCDScaleUp_Norm_Name, bkg_qqzz_QCDScaleUp_NormVal)
      bkg_qqzz_QCDScaleDown_Norm_Name = "bkg_qqzz_widthmass_QCDScaleDownNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_QCDScaleDown_Norm = ROOT.RooRealVar(bkg_qqzz_QCDScaleDown_Norm_Name, bkg_qqzz_QCDScaleDown_Norm_Name, bkg_qqzz_QCDScaleDown_NormVal)

      qqZZ_EWK_Syst.setVal(1)
      bkg_qqzz_EWKcorrUp_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      qqZZ_EWK_Syst.setVal(-1)
      bkg_qqzz_EWKcorrDown_NormVal = bkg_qqzz_mass.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass)).getVal()
      qqZZ_EWK_Syst.setVal(0)
      bkg_qqzz_EWKcorrUp_NormVal = bkg_qqzz_EWKcorrUp_NormVal/bkg_qqzz_Nominal_NormVal
      bkg_qqzz_EWKcorrDown_NormVal = bkg_qqzz_EWKcorrDown_NormVal/bkg_qqzz_Nominal_NormVal
      bkg_qqzz_EWKcorrUp_Norm_Name = "bkg_qqzz_widthmass_EWKcorrUpNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_EWKcorrUp_Norm = ROOT.RooRealVar(bkg_qqzz_EWKcorrUp_Norm_Name, bkg_qqzz_EWKcorrUp_Norm_Name, bkg_qqzz_EWKcorrUp_NormVal)
      bkg_qqzz_EWKcorrDown_Norm_Name = "bkg_qqzz_widthmass_EWKcorrDownNorm_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      bkg_qqzz_EWKcorrDown_Norm = ROOT.RooRealVar(bkg_qqzz_EWKcorrDown_Norm_Name, bkg_qqzz_EWKcorrDown_Norm_Name, bkg_qqzz_EWKcorrDown_NormVal)

      morphVarListqqZZ_norm = ROOT.RooArgList()
      MorphNormList_qqZZ_norm = ROOT.RooArgList()
      morphVarListqqZZ_norm.add(qqZZ_Scale_Syst)
      morphVarListqqZZ_norm.add(qqZZ_EWK_Syst)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_Nominal_Norm)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_QCDScaleUp_Norm)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_QCDScaleDown_Norm)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_EWKcorrUp_Norm)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_EWKcorrDown_Norm)
      asympowname = "Asymquad_qqZZ_QCD_EWK_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      if useDjet!=0:
      morphVarListqqZZ_norm.add(CMS_zz4l_qqzz_djet_syst)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_DjetUp_Norm)
      MorphNormList_qqZZ_norm.add(bkg_qqzz_DjetDown_Norm)
      asympowname = "Asymquad_qqZZ_Djet_QCD_EWK_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      thetaSyst_combined_qqZZ_norm = ROOT.AsymQuad(asympowname, asympowname, MorphNormList_qqZZ_norm, morphVarListqqZZ_norm, 1.0)
      bkg_qqzz_norm = ROOT.RooFormulaVar("bkg_qqzz_norm", "TMath::Max(@0,1e-15)", ROOT.RooArgList(thetaSyst_combined_qqZZ_norm))

      if DEBUG:
      print "bkg_qqZZ_Nominal:",bkg_qqzz_Nominal_NormVal
      print "bkg_qqZZ_QCD Up Ratio:",bkg_qqzz_QCDScaleUp_Norm.getVal()
      print "bkg_qqZZ_QCD Dn Ratio:",bkg_qqzz_QCDScaleDown_Norm.getVal()
      print "bkg_qqZZ_EWK Up Ratio:",bkg_qqzz_EWKcorrUp_Norm.getVal()
      print "bkg_qqZZ_EWK Dn Ratio:",bkg_qqzz_EWKcorrDown_Norm.getVal()
      print "bkg_qqZZ_Djet Up Ratio:",bkg_qqzz_DjetUp_Norm.getVal()
      print "bkg_qqZZ_Djet Dn Ratio:",bkg_qqzz_DjetDown_Norm.getVal()



      TemplateName = "qqzz_TempDataHist_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      qqzz_TempDataHist = ROOT.RooDataHist(
      TemplateName, TemplateName,
      ROOT.RooArgList(
      CMS_zz4l_widthMass,
      CMS_zz4l_widthKD
      ),
      Bkg_T
      )
      PdfName = "qqzz_TemplatePdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      qqzz_TemplatePdf = ROOT.RooHistPdf(
      PdfName, PdfName,
      ROOT.RooArgSet(
      CMS_zz4l_widthMass,
      CMS_zz4l_widthKD
      ),
      qqzz_TempDataHist
      )
      bkg_qqzz = ROOT.RooProdPdf(
      "bkg_qqzz", "bkg_qqzz",
      ROOT.RooArgSet( bkg_qqzz_mass ),
      ROOT.RooFit.Conditional(
      ROOT.RooArgSet( qqzz_TemplatePdf ),
      ROOT.RooArgSet( CMS_zz4l_widthKD )
      )
      )
      if self.dimensions == 1:
      qqzz_TempDataHist1 = ROOT.RooDataHist(
      "{0}_proj".format(TemplateName), "{0}_proj".format(TemplateName),
      ROOT.RooArgList( CMS_zz4l_widthMass ),
      bkg_qqzz.createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD").ProjectionX()
      )
      bkg_qqzz = ROOT.RooHistPdf(
      "bkg_qqzz", "bkg_qqzz",
      ROOT.RooArgSet( CMS_zz4l_widthMass ),
      qqzz_TempDataHist1
      )
      elif self.dimensions == 0:
      qqzz_TempDataHist1 = ROOT.RooDataHist(
      "{0}_proj".format(TemplateName), "{0}_proj".format(TemplateName),
      ROOT.RooArgList( CMS_zz4l_widthKD ),
      bkg_qqzz.createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD").ProjectionY()
      )
      bkg_qqzz = ROOT.RooHistPdf(
      "bkg_qqzz", "bkg_qqzz",
      ROOT.RooArgSet( CMS_zz4l_widthKD ),
      qqzz_TempDataHist1
      )

      bkg_qqzz.SetNameTitle("bkg_qqzz", "bkg_qqzz")


      #------------ Z+X PDFs -----------------#

      val_meanL_3P1F = float(self.theInputs['zjetsShape_mean_3P1F'])
      val_sigmaL_3P1F = float(self.theInputs['zjetsShape_sigma_3P1F'])
      val_normL_3P1F = float(self.theInputs['zjetsShape_norm_3P1F'])

      val_meanL_2P2F = float(self.theInputs['zjetsShape_mean_2P2F'])
      val_sigmaL_2P2F = float(self.theInputs['zjetsShape_sigma_2P2F'])
      val_normL_2P2F = float(self.theInputs['zjetsShape_norm_2P2F'])
      val_pol0_2P2F = float(self.theInputs['zjetsShape_pol0_2P2F'])
      val_pol1_2P2F = float(self.theInputs['zjetsShape_pol1_2P2F'])

      val_meanL_2P2F_2 = float(self.theInputs['zjetsShape_mean_2P2F_2e2mu'])
      val_sigmaL_2P2F_2 = float(self.theInputs['zjetsShape_sigma_2P2F_2e2mu'])
      val_normL_2P2F_2 = float(self.theInputs['zjetsShape_norm_2P2F_2e2mu'])

      TemplateName = "zjet_TempDataHist_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHist = ROOT.RooDataHist(
      TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), Bkg_ZX)
      PdfName = "zjet_TemplatePdf_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_TemplatePdfNominal = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), zjet_TempDataHist)

      TemplateName = "zjet_TempDataHist_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHistUp = ROOT.RooDataHist(
      TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), Bkg_ZX_Up)
      PdfName = "zjet_TemplatePdf_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_TemplatePdfUp = ROOT.RooHistPdf(PdfName, PdfName, ROOT.RooArgSet(
      CMS_zz4l_widthMass, CMS_zz4l_widthKD), zjet_TempDataHistUp)

      TemplateName = "zjet_TempDataHist_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHistDown = ROOT.RooDataHist(
      TemplateName, TemplateName, ROOT.RooArgList(CMS_zz4l_widthMass, CMS_zz4l_widthKD), Bkg_ZX_Down)
      PdfName = "zjet_TemplatePdf_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_TemplatePdfDown = ROOT.RooHistPdf(PdfName, PdfName, ROOT.RooArgSet(
      CMS_zz4l_widthMass, CMS_zz4l_widthKD), zjet_TempDataHistDown)

      # Add a RooProdPDF after each iteration of bkg_zjets_mass with
      # bkg_zjets_mass = Prod(bkg_zjets_mass,analytic form)

      bkg_zjets_mass_Djet_ratio_Name = "bkg_zjets_mass_Djet_ratio_{0:.0f}_{1:.0f}_{2:.0f}".format(self.sqrts, self.channel, useDjet)
      bkg_zjets_mass_Djet_shape_Name = "bkg_zjets_mass_Djet_shape_{0:.0f}_{1:.0f}_{2:.0f}".format(self.sqrts, self.channel, useDjet)
      bkg_zjets_mass_Djet_shape_FormulaCore = "@0"
      bkg_zjets_mass_Djet_shape_Formula = "1"
      bkg_zjets_mass_Djet_ratio = ROOT.RooRealVar(bkg_zjets_mass_Djet_ratio_Name, "bkg_zjets_mass_Djet_ratio", 1,0,1)
      if useDjet != 0:
      # code for analytic form of Djet > 0.5 cut
      if(self.sqrts == 7):
      bkg_zjets_mass_Djet_ratio.setVal(9.94527e-03)
      if(self.sqrts == 8):
      bkg_zjets_mass_Djet_ratio.setVal(1.00038e-02)
      if useDjet==1:
      bkg_zjets_mass_Djet_shape_FormulaCore = "(1 + (@1/(1.-@1))*(1-TMath::Max(@0,0)))"
      bkg_zjets_mass_Djet_shape_Formula = bkg_zjets_mass_Djet_shape_FormulaCore
      bkg_zjets_mass_Djet_ratio.setConstant(True)
      bkg_zjets_djet_arglist = ROOT.RooArgList(
      thetaSyst_djet_zjets_norm
      )
      if useDjet==0:
      bkg_zjets_djet_arglist.removeAll()
      elif useDjet==1:
      bkg_zjets_djet_arglist.add(bkg_zjets_mass_Djet_ratio)
      bkg_zjets_norm = ROOT.RooFormulaVar("bkg_zjets_norm", bkg_zjets_mass_Djet_shape_Formula, bkg_zjets_djet_arglist)




      if (self.channel == self.ID_4mu):
      name = "mlZjet_width_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      mlZjet = ROOT.RooRealVar(name, "mean landau Zjet", val_meanL_2P2F)
      name = "slZjet_width_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      slZjet = ROOT.RooRealVar(
      name, "sigma landau Zjet", val_sigmaL_2P2F)
      print "mean 4mu: ", mlZjet.getVal()
      print "sigma 4mu: ", slZjet.getVal()

      bkg_zjets_mass = ROOT.RooLandau(
      "bkg_zjetsTmp_width_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp", CMS_zz4l_widthMass, mlZjet, slZjet)

      bkg_zjets_Nominal = ROOT.RooProdPdf("bkg_zjets_Nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Nominal", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfNominal), ROOT.RooArgSet(CMS_zz4l_widthKD)))
      bkg_zjets_Up = ROOT.RooProdPdf("bkg_zjets_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Up", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfUp), ROOT.RooArgSet(CMS_zz4l_widthKD)))
      bkg_zjets_Down = ROOT.RooProdPdf("bkg_zjets_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Down", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfDown), ROOT.RooArgSet(CMS_zz4l_widthKD)))

      elif (self.channel == self.ID_4e):

      summa = val_normL_2P2F + val_normL_3P1F

      name = "mlZjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      mlZjet_2p2f = ROOT.RooRealVar(
      name, "mean landau Zjet 2p2f", val_meanL_2P2F)
      name = "slZjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      slZjet_2p2f = ROOT.RooRealVar(
      name, "sigma landau Zjet 2p2f", val_sigmaL_2P2F)
      name = "nlZjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      nlZjet_2p2f = ROOT.RooRealVar(
      name, "norm landau Zjet 2p2f", val_normL_2P2F / summa)
      name = "p0Zjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      p0Zjet_2p2f = ROOT.RooRealVar(name, "p0 Zjet 2p2f", val_pol0_2P2F)
      name = "p1Zjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      p1Zjet_2p2f = ROOT.RooRealVar(name, "p1 Zjet 2p2f", val_pol1_2P2F)
      print "mean 2p2f 4e: ", mlZjet_2p2f.getVal()
      print "sigma 2p2f 4e: ", slZjet_2p2f.getVal()
      print "norm 2p2f 4e: ", nlZjet_2p2f.getVal()
      print "pol0 2p2f 4e: ", p0Zjet_2p2f.getVal()
      print "pol1 2p2f 4e: ", p1Zjet_2p2f.getVal()
      bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp_2p2f", "(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))", RooArgList(
      CMS_zz4l_widthMass, mlZjet_2p2f, slZjet_2p2f, p0Zjet_2p2f, p1Zjet_2p2f))

      name = "mlZjet_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      mlZjet_3p1f = ROOT.RooRealVar(
      name, "mean landau Zjet 3p1f", val_meanL_3P1F)
      name = "slZjet_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      slZjet_3p1f = ROOT.RooRealVar(
      name, "sigma landau Zjet 3p1f", val_sigmaL_3P1F)
      name = "nlZjet_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      nlZjet_3p1f = ROOT.RooRealVar(
      name, "norm landau Zjet 3p1f", val_normL_3P1F / summa)
      print "mean 3p1f 4e: ", mlZjet_3p1f.getVal()
      print "sigma 3p1f 4e: ", slZjet_3p1f.getVal()
      print "norm 3p1f 4e: ", nlZjet_3p1f.getVal()

      bkg_zjets_3p1f = ROOT.RooLandau(
      "bkg_zjetsTmp_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp_3p1f", CMS_zz4l_widthMass, mlZjet_3p1f, slZjet_3p1f)

      bkg_zjets_mass = ROOT.RooAddPdf("bkg_zjetsTmp_width_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp", ROOT.RooArgList(
      bkg_zjets_2p2f, bkg_zjets_3p1f), ROOT.RooArgList(nlZjet_2p2f, nlZjet_3p1f))

      bkg_zjets_Nominal = ROOT.RooProdPdf("bkg_zjets_Nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Nominal", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfNominal), ROOT.RooArgSet(CMS_zz4l_widthKD)))
      bkg_zjets_Up = ROOT.RooProdPdf("bkg_zjets_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Up", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfUp), ROOT.RooArgSet(CMS_zz4l_widthKD)))
      bkg_zjets_Down = ROOT.RooProdPdf("bkg_zjets_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Down", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfDown), ROOT.RooArgSet(CMS_zz4l_widthKD)))

      elif (self.channel == self.ID_2e2mu):

      summa = val_normL_2P2F + val_normL_2P2F_2 + val_normL_3P1F

      name = "mlZjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      mlZjet_2p2f = ROOT.RooRealVar(
      name, "mean landau Zjet 2p2f", val_meanL_2P2F)
      name = "slZjet_width_2p2f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      slZjet_2p2f = ROOT.RooRealVar(
      name, "sigma landau Zjet 2p2f", val_sigmaL_2P2F)
      name = "nlZjet_width_2p2f_{0:.0f}_{1:.0f}".format(
      self.channel, self.sqrts, useDjet)

      nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F/summa)
      name = "p0Zjet_width_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
      p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",val_pol0_2P2F)
      name = "p1Zjet_width_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
      p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",val_pol1_2P2F)
      print "mean 2p2f 2mu2e: ",mlZjet_2p2f.getVal()
      print "sigma 2p2f 2mu2e: ",slZjet_2p2f.getVal()
      print "norm 2p2f 2mu2e: ",nlZjet_2p2f.getVal()
      print "pol0 2p2f 2mu2e: ",p0Zjet_2p2f.getVal()
      print "pol1 2p2f 2mu2e: ",p1Zjet_2p2f.getVal()
      bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_width_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",RooArgList(CMS_zz4l_widthMass,mlZjet_2p2f,slZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))

      name = "mlZjet_width_2p2f_2_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      mlZjet_2p2f_2 = ROOT.RooRealVar(
      name, "mean landau Zjet 2p2f 2e2mu", val_meanL_2P2F_2)
      name = "slZjet_width_2p2f_2_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      slZjet_2p2f_2 = ROOT.RooRealVar(
      name, "sigma landau Zjet 2p2f 2e2mu", val_sigmaL_2P2F_2)
      name = "nlZjet_width_2p2f_2_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      nlZjet_2p2f_2 = ROOT.RooRealVar(
      name, "norm landau Zjet 2p2f 2e2mu", val_normL_2P2F_2 / summa)
      print "mean 2p2f 2e2mu: ", mlZjet_2p2f_2.getVal()
      print "sigma 2p2f 2e2mu: ", slZjet_2p2f_2.getVal()
      print "norm 2p2f 2e2mu: ", nlZjet_2p2f_2.getVal()
      bkg_zjets_2p2f_2 = ROOT.RooLandau(
      "bkg_zjetsTmp_width_2p2f_2_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp_2p2f_2", CMS_zz4l_widthMass, mlZjet_2p2f_2, slZjet_2p2f_2)

      name = "mlZjet_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      mlZjet_3p1f = ROOT.RooRealVar(
      name, "mean landau Zjet 3p1f", val_meanL_3P1F)
      name = "slZjet_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      slZjet_3p1f = ROOT.RooRealVar(
      name, "sigma landau Zjet 3p1f", val_sigmaL_3P1F)
      name = "nlZjet_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      nlZjet_3p1f = ROOT.RooRealVar(
      name, "norm landau Zjet 3p1f", val_normL_3P1F / summa)
      print "mean 3p1f 2mu2e: ", mlZjet_3p1f.getVal()
      print "sigma 3p1f 2mu2e: ", slZjet_3p1f.getVal()
      print "norm 3p1f 2mu2e: ", nlZjet_3p1f.getVal()
      bkg_zjets_3p1f = ROOT.RooLandau(
      "bkg_zjetsTmp_width_3p1f_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp_3p1f", CMS_zz4l_widthMass, mlZjet_3p1f, slZjet_3p1f)

      bkg_zjets_mass = ROOT.RooAddPdf("bkg_zjetsTmp_width_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjetsTmp", ROOT.RooArgList(
      bkg_zjets_2p2f, bkg_zjets_3p1f, bkg_zjets_2p2f_2), ROOT.RooArgList(nlZjet_2p2f, nlZjet_3p1f, nlZjet_2p2f_2))

      bkg_zjets_Nominal = ROOT.RooProdPdf("bkg_zjets_Nominal_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Nominal", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfNominal), ROOT.RooArgSet(CMS_zz4l_widthKD)))
      bkg_zjets_Up = ROOT.RooProdPdf("bkg_zjets_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Up", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfUp), ROOT.RooArgSet(CMS_zz4l_widthKD)))
      bkg_zjets_Down = ROOT.RooProdPdf("bkg_zjets_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet), "bkg_zjets_Down", ROOT.RooArgSet(
      bkg_zjets_mass), ROOT.RooFit.Conditional(ROOT.RooArgSet(zjet_TemplatePdfDown), ROOT.RooArgSet(CMS_zz4l_widthKD)))


      DataName = "ZX_FullDataHist_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      PdfName = "ZX_FullPdf_{0:.0f}_{1:.0f}_{2:.0f}".format(self.channel, self.sqrts, useDjet)
      zjet_DataHistNominal = RooDataHist(
      DataName, DataName, ROOT.RooArgList(
      CMS_zz4l_widthMass, CMS_zz4l_widthKD),
      bkg_zjets_Nominal.createHistogram("{0},{1}".format(CMS_zz4l_widthMass.GetName(), CMS_zz4l_widthKD.GetName())))
      zjet_HistPdfNominaltemp = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), zjet_DataHistNominal)

      DataName = "ZX_FullDataHist_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      PdfName = "ZX_FullPdf_Up_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_DataHistUp = RooDataHist(
      DataName, DataName, ROOT.RooArgList(
      CMS_zz4l_widthMass, CMS_zz4l_widthKD),
      bkg_zjets_Up.createHistogram("{0},{1}".format(CMS_zz4l_widthMass.GetName(), CMS_zz4l_widthKD.GetName())))
      zjet_HistPdfUptemp = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), zjet_DataHistUp)

      DataName = "ZX_FullDataHist_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      PdfName = "ZX_FullPdf_Down_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      zjet_DataHistDown = RooDataHist(
      DataName, DataName, ROOT.RooArgList(
      CMS_zz4l_widthMass, CMS_zz4l_widthKD),
      bkg_zjets_Down.createHistogram("{0},{1}".format(CMS_zz4l_widthMass.GetName(), CMS_zz4l_widthKD.GetName())))
      zjet_HistPdfDowntemp = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD), zjet_DataHistDown)

      if self.dimensions == 0:
      PdfName = "ZX_FullPdf_Nominal1_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      TemplateName = "zjet_TempDataHist1_{0:.0f}_{1:.0f}_{2:.0f}_Nominal".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHist1_Nominal = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(
      CMS_zz4l_widthKD), bkg_zjets_Nominal.createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD").ProjectionY())
      zjet_HistPdfNominal = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), zjet_TempDataHist1_Nominal)

      PdfName = "ZX_FullPdf_Up1_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      TemplateName = "zjet_TempDataHist1_{0:.0f}_{1:.0f}_{2:.0f}_Up".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHist1_Up = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(
      CMS_zz4l_widthKD), bkg_zjets_Up.createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD").ProjectionY())
      zjet_HistPdfUp = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), zjet_TempDataHist1_Up)

      PdfName = "ZX_FullPdf_Down1_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      TemplateName = "zjet_TempDataHist1_{0:.0f}_{1:.0f}_{2:.0f}_Down".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHist1_Down = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(
      CMS_zz4l_widthKD), bkg_zjets_Down.createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD").ProjectionY())
      zjet_HistPdfDown = ROOT.RooHistPdf(
      PdfName, PdfName, ROOT.RooArgSet(CMS_zz4l_widthKD), zjet_TempDataHist1_Down)

      if self.dimensions == 2:
      zjet_HistPdfNominal = zjet_HistPdfNominaltemp
      zjet_HistPdfUp = zjet_HistPdfUptemp
      zjet_HistPdfDown = zjet_HistPdfDowntemp
      if self.dimensions == 1:
      PdfName = "ZX_FullPdf_Nominal1_{0:.0f}_{1:.0f}_{2:.0f}".format(
      self.channel, self.sqrts, useDjet)
      TemplateName = "zjet_TempDataHist1_{0:.0f}_{1:.0f}_{2:.0f}_Nominal".format(
      self.channel, self.sqrts, useDjet)
      zjet_TempDataHist1_Nominal = ROOT.RooDataHist(TemplateName, TemplateName, ROOT.RooArgList(
      CMS_zz4l_widthMass), bkg_zjets_Nominal.createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD").ProjectionX())
      bkg_zjets = ROOT.RooHistPdf(
      "bkg_zjets", "bkg_zjets", ROOT.RooArgSet(CMS_zz4l_widthMass), zjet_TempDataHist1_Nominal)
      else:
      CMS_zz4l_ZXshape_syst = ROOT.RooRealVar(
      "CMS_zz4l_ZXshape_syst", "CMS_zz4l_ZXshape_syst", 0.0, -1, 1)
      morphVarListZX = ROOT.RooArgList()
      morphVarListZX.add(CMS_zz4l_ZXshape_syst)
      MorphList_ZX = ROOT.RooArgList()
      MorphList_ZX.add(zjet_HistPdfNominal)
      MorphList_ZX.add(zjet_HistPdfUp)
      MorphList_ZX.add(zjet_HistPdfDown)

      bkg_zjets = ROOT.VerticalInterpPdf("bkg_zjets", "bkg_zjets", MorphList_ZX, morphVarListZX)


      ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##

      if DEBUG:
      testhandle = ROOT.TFile("{0}/figs/xcheckPlots_{2}_{1}_{3}.root".format(self.theOutputDir,self.appendName, self.sqrts, useDjet),"recreate")

      print "Starting to plot signal for sanity checks"
      print "Plot: ggZZ Norm vs GGsm"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_norm.GetName(),x.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = ggZZpdf_norm.createHistogram("htemp",x)
      histo.SetLineWidth(2)
      histo.Draw()
      intpdftmp = ggZZpdf_Nominal.createIntegral(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
      histo2 = intpdftmp.createHistogram("mytemp",x)
      histo2.SetLineColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()
      print "Plot: ggZZ Norm vs fai1"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_norm.GetName(),fai1.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = ggZZpdf_norm.createHistogram("htemp",fai1)
      histo.Draw()
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()
      print "Plot: VBF Norm vs GGsm"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,VBFpdf_norm.GetName(),x.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = VBFpdf_norm.createHistogram("htemp",x)
      histo.Draw()
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()
      print "Plot: VBF Norm vs fai1"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,VBFpdf_norm.GetName(),fai1.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = VBFpdf_norm.createHistogram("htemp",fai1)
      histo.Draw()
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: ggZZ PDF on mZZ vs KD"
      canvasname = "c_{3}_vs_{4}_{5}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet, ggZZpdf.GetName(),CMS_zz4l_widthMass.GetName(),CMS_zz4l_widthKD.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = ggZZpdf.createHistogram("htemp",CMS_zz4l_widthMass,ROOT.RooFit.YVar(CMS_zz4l_widthKD),ROOT.RooFit.ZVar(x))
      for binx in range(1,histo.GetNbinsX()+1) :
      for biny in range(1,histo.GetNbinsY()+1) :
      for binz in range(1,histo.GetNbinsZ()+1) :
      bincontent = histo.GetBinContent(binx,biny,binz)
      if bincontent <= 0: print "Bin content {0:.5f} at mass bin = {1} KD bin = {2} GGsm val = {3} is not positive.".format(bincontent,binx,biny,histo.GetZaxis().GetBinCenter(binz))
      defaultx = x.getVal()
      x.setVal(30)
      histo = ggZZpdf.createHistogram("htemp",CMS_zz4l_widthMass,ROOT.RooFit.YVar(CMS_zz4l_widthKD))
      histo.Draw("colz")
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()
      x.setVal(defaultx)


      print "Plot: qqZZ Norm vs EWK"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_norm.GetName(),qqZZ_EWK_Syst.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_norm.createHistogram("htemp",qqZZ_EWK_Syst)
      histo.Draw()
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()
      print "Plot: qqZZ Norm vs QCD"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_norm.GetName(),qqZZ_Scale_Syst.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_norm.createHistogram("htemp",qqZZ_Scale_Syst)
      histo.Draw()
      #            canvasfullName = "{0}/figs/{1}.root".format(self.theOutputDir, canvasname)
      #            ctest.SaveAs(canvasfullName)
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: ggZZ Norm vs CMS_zz4l_APscale_syst"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_norm.GetName(),CMS_zz4l_APscale_syst.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = ggZZpdf_norm.createHistogram("htemp",CMS_zz4l_APscale_syst)
      histo.SetLineWidth(2)
      histo2 = histo.Clone("mytemp")
      histo2.SetLineColor(kRed)
      histo2.SetLineStyle(7)
      histo3 = histo.Clone("mytemp")
      histo3.SetLineColor(kBlue)
      histo3.SetLineStyle(2)
      for bin in range(1,histo2.GetNbinsX()+1):
      CMS_zz4l_APscale_syst.setVal(histo2.GetBinCenter(bin))
      histo2.SetBinContent(bin,ggZZpdf.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD)))
      ceff = CMS_zz4l_APscale_syst.getVal()
      ceff2 = ceff*ceff
      termnull = integral_Sig_T_1+integral_Sig_T_2+integral_Sig_T_4
      termup = integral_Sig_T_1_Up_QCD+integral_Sig_T_2_Up_QCD+integral_Sig_T_4_Up_QCD
      termdn = integral_Sig_T_1_Down_QCD+integral_Sig_T_2_Down_QCD+integral_Sig_T_4_Down_QCD
      val = termnull+ceff2*((termdn+termup)/2.-termnull)+ceff*(termup-termdn)/2.
      histo3.SetBinContent(bin,val)
      # BLACK AND RED SHOULD BE ON TOP OF EACH OTHER, ALWAYS!!!
      CMS_zz4l_APscale_syst.setVal(0)
      histo.Scale(histo3.GetBinContent(histo3.GetXaxis().FindBin(0))/histo.GetBinContent(histo.GetXaxis().FindBin(0)))
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: ggZZ Norm vs CMS_zz4l_pdf_gg_syst"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_norm.GetName(),CMS_zz4l_pdf_gg_syst.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = ggZZpdf_norm.createHistogram("htemp",CMS_zz4l_pdf_gg_syst)
      histo.SetLineWidth(2)
      histo2 = histo.Clone("mytemp")
      histo2.SetLineColor(kRed)
      histo2.SetLineStyle(7)
      histo3 = histo.Clone("mytemp")
      histo3.SetLineColor(kBlue)
      histo3.SetLineStyle(2)
      for bin in range(1,histo2.GetNbinsX()+1):
      CMS_zz4l_pdf_gg_syst.setVal(histo2.GetBinCenter(bin))
      histo2.SetBinContent(bin,ggZZpdf.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD)))
      ceff = CMS_zz4l_pdf_gg_syst.getVal()
      ceff2 = ceff*ceff
      termnull = integral_Sig_T_1+integral_Sig_T_2+integral_Sig_T_4
      termup = integral_Sig_T_1_Up_PDF+integral_Sig_T_2_Up_PDF+integral_Sig_T_4_Up_PDF
      termdn = integral_Sig_T_1_Down_PDF+integral_Sig_T_2_Down_PDF+integral_Sig_T_4_Down_PDF
      val = termnull+ceff2*((termdn+termup)/2.-termnull)+ceff*(termup-termdn)/2.
      histo3.SetBinContent(bin,val)
      # BLACK AND RED SHOULD BE ON TOP OF EACH OTHER, ALWAYS!!!
      CMS_zz4l_pdf_gg_syst.setVal(0)
      histo.Scale(histo3.GetBinContent(histo3.GetXaxis().FindBin(0))/histo.GetBinContent(histo.GetXaxis().FindBin(0)))
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: VBF Norm vs CMS_zz4l_VBFscale_syst"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,VBFpdf_norm.GetName(),CMS_zz4l_VBFscale_syst.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = VBFpdf_norm.createHistogram("htemp",CMS_zz4l_VBFscale_syst)
      histo.SetLineWidth(2)
      histo2 = histo.Clone("mytemp")
      histo2.SetLineColor(kRed)
      histo2.SetLineStyle(7)
      histo3 = histo.Clone("mytemp")
      histo3.SetLineColor(kBlue)
      histo3.SetLineStyle(2)
      for bin in range(1,histo2.GetNbinsX()+1):
      CMS_zz4l_VBFscale_syst.setVal(histo2.GetBinCenter(bin))
      histo2.SetBinContent(bin,VBFpdf.getNorm(ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD)))
      ceff = CMS_zz4l_VBFscale_syst.getVal()
      ceff2 = ceff*ceff
      termnull = integral_VBF_T_1+integral_VBF_T_2+integral_VBF_T_4
      termup = integral_VBF_T_1_Up+integral_VBF_T_2_Up+integral_VBF_T_4_Up
      termdn = integral_VBF_T_1_Down+integral_VBF_T_2_Down+integral_VBF_T_4_Down
      val = termnull+ceff2*((termdn+termup)/2.-termnull)+ceff*(termup-termdn)/2.
      histo3.SetBinContent(bin,val)
      # BLACK AND RED SHOULD BE ON TOP OF EACH OTHER, ALWAYS!!!
      CMS_zz4l_VBFscale_syst.setVal(0)
      histo.Scale(histo3.GetBinContent(histo3.GetXaxis().FindBin(0))/histo.GetBinContent(histo.GetXaxis().FindBin(0)))
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: ggZZ pdf proj (QCD Up)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_Up.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = ggZZpdf_Up.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = Sig_T_1_Up_QCD.Clone("mytemp1")
      histo2.Add(Sig_T_2_Up_QCD.Clone("mytemp2"),10)
      histo2.Add(Sig_T_4_Up_QCD.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      CMS_zz4l_APscale_syst.setVal(1)
      histo3 = ggZZpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      CMS_zz4l_APscale_syst.setVal(0)
      x.setVal(1)

      print "Plot: ggZZ pdf proj (QCD Down)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_Down.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = ggZZpdf_Down.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = Sig_T_1_Down_QCD.Clone("mytemp1")
      histo2.Add(Sig_T_2_Down_QCD.Clone("mytemp2"),10)
      histo2.Add(Sig_T_4_Down_QCD.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      CMS_zz4l_APscale_syst.setVal(-1)
      histo3 = ggZZpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      #            histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      CMS_zz4l_APscale_syst.setVal(0)
      x.setVal(1)

      print "Plot: ggZZ pdf proj (PDF Up)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_Up_pdf.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = ggZZpdf_Up_pdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = Sig_T_1_Up_PDF.Clone("mytemp1")
      histo2.Add(Sig_T_2_Up_PDF.Clone("mytemp2"),10)
      histo2.Add(Sig_T_4_Up_PDF.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      CMS_zz4l_pdf_gg_syst.setVal(1)
      histo3 = ggZZpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      CMS_zz4l_pdf_gg_syst.setVal(0)
      x.setVal(1)

      print "Plot: ggZZ pdf proj (PDF Down)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_Down_pdf.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = ggZZpdf_Down_pdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = Sig_T_1_Down_PDF.Clone("mytemp1")
      histo2.Add(Sig_T_2_Down_PDF.Clone("mytemp2"),10)
      histo2.Add(Sig_T_4_Down_PDF.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      CMS_zz4l_pdf_gg_syst.setVal(-1)
      histo3 = ggZZpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      #            histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      CMS_zz4l_pdf_gg_syst.setVal(0)
      x.setVal(1)

      print "Plot: ggZZ pdf proj (Nominal)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,ggZZpdf_Nominal.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = ggZZpdf_Nominal.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = Sig_T_1.Clone("mytemp1")
      histo2.Add(Sig_T_2.Clone("mytemp2"),10)
      histo2.Add(Sig_T_4.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      histo3 = ggZZpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      x.setVal(1)

      print "Plot: VBF pdf proj (Up)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,VBFpdf_Up.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = VBFpdf_Up.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = VBF_T_1_Up.Clone("mytemp1")
      histo2.Add(VBF_T_2_Up.Clone("mytemp2"),10)
      histo2.Add(VBF_T_4_Up.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      CMS_zz4l_VBFscale_syst.setVal(1)
      histo3 = VBFpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      CMS_zz4l_VBFscale_syst.setVal(0)
      x.setVal(1)

      print "Plot: VBF pdf proj (Down)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,VBFpdf_Down.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = VBFpdf_Down.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = VBF_T_1_Down.Clone("mytemp1")
      histo2.Add(VBF_T_2_Down.Clone("mytemp2"),10)
      histo2.Add(VBF_T_4_Down.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      CMS_zz4l_VBFscale_syst.setVal(-1)
      histo3 = VBFpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      #            histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      CMS_zz4l_VBFscale_syst.setVal(0)
      x.setVal(1)

      print "Plot: VBF pdf proj (Nominal)"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,VBFpdf_Nominal.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      x.setVal(10)
      histo = VBFpdf_Nominal.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo.Draw()
      histo2 = VBF_T_1.Clone("mytemp1")
      histo2.Add(VBF_T_2.Clone("mytemp2"),10)
      histo2.Add(VBF_T_4.Clone("mytemp3"),sqrt(10))
      histo2 = histo2.ProjectionX()
      histo2.SetLineColor(kRed)
      histo2.SetMarkerColor(kRed)
      histo2.Draw("same")
      histo2.Scale(histo.Integral()/histo2.Integral())
      histo3 = VBFpdf.createProjection(ROOT.RooArgSet(CMS_zz4l_widthKD)).createHistogram("htempfull",CMS_zz4l_widthMass)
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      histo3.Scale(histo.Integral()/histo3.Integral())
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      x.setVal(1)

      print "Plot: qqZZ additional wrt QCDScale"
      canvasname = "c_{3}_vs_{4}_wrtQCDScale_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_mass_shape.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_mass_shape.createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo2 = histo.Clone("htemp_Up")
      histo2.SetLineWidth(2)
      histo2.SetLineStyle(7)
      histo2.SetLineColor(kRed)
      histo3 = histo.Clone("htemp_Dn")
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      for bin in range(1,histo.GetNbinsX()+1):
      qqZZ_Scale_Syst.setVal(0)
      CMS_zz4l_widthMass.setVal(histo.GetBinCenter(bin))
      histo.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      qqZZ_Scale_Syst.setVal(1)
      histo2.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      qqZZ_Scale_Syst.setVal(-1)
      histo3.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      qqZZ_Scale_Syst.setVal(0)
      histo.GetYaxis().SetRangeUser(0,2)
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      qqZZ_Scale_Syst.setVal(0)


      print "Plot: qqZZ additional wrt EWKcorr"
      canvasname = "c_{3}_vs_{4}_wrtEWKcorr_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_mass_shape.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_mass_shape.createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo2 = histo.Clone("htemp_Up")
      histo2.SetLineWidth(2)
      histo2.SetLineStyle(7)
      histo2.SetLineColor(kRed)
      histo3 = histo.Clone("htemp_Dn")
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      for bin in range(1,histo.GetNbinsX()+1):
      qqZZ_EWK_Syst.setVal(0)
      CMS_zz4l_widthMass.setVal(histo.GetBinCenter(bin))
      histo.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      qqZZ_EWK_Syst.setVal(1)
      histo2.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      qqZZ_EWK_Syst.setVal(-1)
      histo3.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      qqZZ_EWK_Syst.setVal(0)
      histo.GetYaxis().SetRangeUser(0,2)
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      qqZZ_EWK_Syst.setVal(0)

      print "Plot: qqZZ additional wrt Djet"
      canvasname = "c_{3}_vs_{4}_wrtDjet_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_mass_shape.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_mass_shape.createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo2 = histo.Clone("htemp_Up")
      histo2.SetLineWidth(2)
      histo2.SetLineStyle(7)
      histo2.SetLineColor(kRed)
      histo3 = histo.Clone("htemp_Dn")
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(2)
      histo3.SetLineColor(kBlue)
      for bin in range(1,histo.GetNbinsX()+1):
      CMS_zz4l_qqzz_djet_syst.setVal(0)
      bincenter = histo.GetBinCenter(bin)
      CMS_zz4l_widthMass.setVal(bincenter)
      histo.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      CMS_zz4l_qqzz_djet_syst.setVal(1)
      histo2.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      CMS_zz4l_qqzz_djet_syst.setVal(-1)
      histo3.SetBinContent(bin,bkg_qqzz_mass_shape.getVal())
      CMS_zz4l_qqzz_djet_syst.setVal(0)
      histo.GetYaxis().SetRangeUser(0,2)
      histo2.Divide(histo)
      histo3.Divide(histo)
      histo.Divide(histo)
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      testhandle.WriteTObject(ctest)
      ctest.Close()
      qqZZ_Scale_Syst.setVal(0)


      print "Plot: qqZZ PDF wrt EWKcorr"
      canvasname = "c_{3}_vs_{4}_wrtEWKcorr_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_mass.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_mass.createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo2 = histo.Clone("htemp_Up")
      histo2.SetLineWidth(2)
      histo2.SetLineStyle(7)
      histo2.SetLineColor(kRed)
      histo3 = histo.Clone("htemp_Dn")
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(7)
      histo3.SetLineColor(kBlue)
      for bin in range(1,histo.GetNbinsX()+1):
      qqZZ_EWK_Syst.setVal(0)
      CMS_zz4l_widthMass.setVal(histo.GetBinCenter(bin))
      histo.SetBinContent(bin,bkg_qqzz_mass.getVal())
      qqZZ_EWK_Syst.setVal(1)
      histo2.SetBinContent(bin,bkg_qqzz_mass.getVal())
      qqZZ_EWK_Syst.setVal(-1)
      histo3.SetBinContent(bin,bkg_qqzz_mass.getVal())
      qqZZ_EWK_Syst.setVal(0)
      histon = bkg_qqzz_mass_temp.createHistogram("htempn",CMS_zz4l_widthMass)
      histon.SetLineWidth(2)
      histon.SetLineStyle(2)
      histon2 = histon.Clone("htempn_Up")
      histon2.SetLineWidth(2)
      histon2.SetLineStyle(2)
      histon2.SetLineColor(kRed)
      histon3 = histon.Clone("htempn_Dn")
      histon3.SetLineWidth(2)
      histon3.SetLineStyle(2)
      histon3.SetLineColor(kBlue)
      for bin in range(1,histon.GetNbinsX()+1):
      #                CMS_qqzzbkg_EWK_p0
      #                CMS_qqzzbkg_EWK_p1
      #                CMS_qqzzbkg_EWK_p2
      #                CMS_qqzzbkg_EWK_p3
      #                CMS_qqzzbkg_p0
      #                CMS_qqzzbkg_p1
      bincenter = histon.GetXaxis().GetBinCenter(bin)
      bincontent = histo.GetBinContent(bin)
      systcorr = CMS_qqzzbkg_EWK_p0.getVal() + CMS_qqzzbkg_EWK_p1.getVal()*bincenter + CMS_qqzzbkg_EWK_p2.getVal()*bincenter*bincenter + CMS_qqzzbkg_EWK_p3.getVal()*bincenter*bincenter*bincenter
      histon.SetBinContent(bin,bincontent)
      histon2.SetBinContent(bin,bincontent*systcorr)
      histon3.SetBinContent(bin,bincontent*(2-systcorr))
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      histon.Draw("same")
      histon2.Draw("same")
      histon3.Draw("same")
      print "EWK up / EWK nominal from histo: ",(histo2.Integral()/histo.Integral())
      print "EWK dn / EWK nominal from histo: ",(histo3.Integral()/histo.Integral())
      print "EWK up / EWK nominal from histon: ",(histon2.Integral()/histon.Integral())
      print "EWK dn / EWK nominal from histon: ",(histon3.Integral()/histon.Integral())
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: qqZZ PDF wrt QCDScale"
      canvasname = "c_{3}_vs_{4}_wrtQCDScale_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_qqzz_mass.GetName(),CMS_zz4l_widthMass.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_qqzz_mass.createHistogram("htemp",CMS_zz4l_widthMass)
      histo.SetLineWidth(2)
      histo.SetLineStyle(7)
      histo2 = histo.Clone("htemp_Up")
      histo2.SetLineWidth(2)
      histo2.SetLineStyle(7)
      histo2.SetLineColor(kRed)
      histo3 = histo.Clone("htemp_Dn")
      histo3.SetLineWidth(2)
      histo3.SetLineStyle(7)
      histo3.SetLineColor(kBlue)
      for bin in range(1,histo.GetNbinsX()+1):
      qqZZ_Scale_Syst.setVal(0)
      CMS_zz4l_widthMass.setVal(histo.GetBinCenter(bin))
      histo.SetBinContent(bin,bkg_qqzz_mass.getVal())
      qqZZ_Scale_Syst.setVal(1)
      histo2.SetBinContent(bin,bkg_qqzz_mass.getVal())
      qqZZ_Scale_Syst.setVal(-1)
      histo3.SetBinContent(bin,bkg_qqzz_mass.getVal())
      qqZZ_Scale_Syst.setVal(0)
      histon = bkg_qqzz_mass_temp.createHistogram("htempn",CMS_zz4l_widthMass)
      histon.SetLineWidth(2)
      histon.SetLineStyle(2)
      histon2 = histon.Clone("htempn_Up")
      histon2.SetLineWidth(2)
      histon2.SetLineStyle(2)
      histon2.SetLineColor(kRed)
      histon3 = histon.Clone("htempn_Dn")
      histon3.SetLineWidth(2)
      histon3.SetLineStyle(2)
      histon3.SetLineColor(kBlue)
      for bin in range(1,histon.GetNbinsX()+1):
      bincenter = histon.GetXaxis().GetBinCenter(bin)
      bincontent = histo.GetBinContent(bin)
      systcorr = CMS_qqzzbkg_p0.getVal() + CMS_qqzzbkg_p1.getVal()*bincenter + CMS_qqzzbkg_p2.getVal()*bincenter*bincenter
      histon.SetBinContent(bin,bincontent)
      histon2.SetBinContent(bin,bincontent*systcorr)
      histon3.SetBinContent(bin,bincontent*(2-systcorr))
      histo.Draw()
      histo2.Draw("same")
      histo3.Draw("same")
      histon.Draw("same")
      histon2.Draw("same")
      histon3.Draw("same")
      print "QCD up / QCD nominal from histo: ",(histo2.Integral()/histo.Integral())
      print "QCD dn / QCD nominal from histo: ",(histo3.Integral()/histo.Integral())
      print "QCD up / QCD nominal from histon: ",(histon2.Integral()/histon.Integral())
      print "QCD dn / QCD nominal from histon: ",(histon3.Integral()/histon.Integral())
      testhandle.WriteTObject(ctest)
      ctest.Close()

      print "Plot: Z+X Norm vs CMS_zz4l_zjets_djet_syst"
      canvasname = "c_{3}_vs_{4}_{1:.0f}TeV_{0}_djet{2}".format(self.appendName, self.sqrts, useDjet,bkg_zjets_norm.GetName(),CMS_zz4l_zjets_djet_syst.GetName())
      ctest = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      ctest.cd()
      histo = bkg_zjets_norm.createHistogram("htemp",CMS_zz4l_zjets_djet_syst)
      histo.SetLineWidth(2)
      histo.Scale(1/0.02)
      histo.Draw()
      testhandle.WriteTObject(ctest)
      ctest.Close()


      testhandle.Close()


      # --------------------------- DATASET --------------------------- ##

      if(USELEGACY == 0):
      dataFileDir = "CMSdata"
      if(USELEGACY == 1):
      dataFileDir = "CMSdata_Legacy"
      if (self.dataAppendDir != ''):
      dataFileDir = "{0}_{1}".format(dataFileDir,self.dataAppendDir)
      dataTreeName = "data_obs"
      if(useDjet == 0):
      dataFileName = "{0}/hzz{1}_{2}.root".format(
      dataFileDir, self.appendName, self.inputlumi)
      if(useDjet == 1):
      dataFileName = "{0}/hzz{1}_{2}_0.root".format(
      dataFileDir, self.appendName, self.inputlumi)
      if(useDjet == 2):
      dataFileName = "{0}/hzz{1}_{2}_1.root".format(
      dataFileDir, self.appendName, self.inputlumi)
      # if (DEBUG):
      print dataFileName, " ", dataTreeName
      data_obs_file = ROOT.TFile(dataFileName)

      print data_obs_file.Get(dataTreeName)

      if not (data_obs_file.Get(dataTreeName)):
      print "File, \"", dataFileName, "\", or tree, \"", dataTreeName, "\", not found"
      print "Exiting..."
      sys.exit()

      data_obs_tree = data_obs_file.Get(dataTreeName)
      data_obs = ROOT.RooDataSet()
      datasetName = "data_obs_{0}".format(self.appendName)

      data_obs = ROOT.RooDataSet(
      datasetName, datasetName, data_obs_tree, ROOT.RooArgSet(CMS_zz4l_widthMass, CMS_zz4l_widthKD))
      data_obs_red = data_obs.reduce(
      "CMS_zz4l_widthMass > {0}".format(self.low_M))

      # --------------------------- WORKSPACE -------------------------- ##
      endsInP5 = False
      tmpMH = self.low_M
      if (math.fabs(math.floor(tmpMH) - self.low_M) > 0.001):
      endsInP5 = True
      if (DEBUG):
      print "ENDS IN P5  ", endsInP5

      name_Shape = ""
      name_ShapeWS = ""
      name_ShapeWS2 = ""

      if (endsInP5):
      if(useDjet == 0):
      name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 1):
      name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_0.txt".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 2):
      name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_1.txt".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      else:
      if(useDjet == 0):
      name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 1):
      name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_0.txt".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 2):
      name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_1.txt".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)

      if (endsInP5):
      if(useDjet == 0):
      name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 1):
      name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_0.input.root".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 2):
      name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_1.input.root".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      else:
      if(useDjet == 0):
      name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 1):
      name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_0.input.root".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)
      if(useDjet == 2):
      name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_1.input.root".format(
      self.theOutputDir, self.low_M, self.appendName, self.sqrts)

      if(useDjet == 0):
      name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(
      self.appendName, self.sqrts)
      if(useDjet == 1):
      name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV_0.input.root".format(
      self.appendName, self.sqrts)
      if(useDjet == 2):
      name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV_1.input.root".format(
      self.appendName, self.sqrts)

      if(DEBUG):
      print name_Shape, "  ", name_ShapeWS2

      w.importClassCode(AsymPow.Class(),True)
      w.importClassCode(AsymQuad.Class(),True)
      w.importClassCode(RooqqZZPdf_v2.Class(), True)
      w.importClassCode(RooFormulaVar.Class(), True)
      w.importClassCode(RooRealFlooredSumPdf.Class(),True)
      w.importClassCode(VerticalInterpPdf.Class(),True)

      getattr(w, 'import')(data_obs_red, ROOT.RooFit.Rename("data_obs"))

      # Import dependencies explicitly first since text2workspace gives irrelevant error messages
      getattr(w, 'import')(self.LUMI, ROOT.RooFit.RecycleConflictNodes())
      getattr(w, 'import')(x, ROOT.RooFit.RecycleConflictNodes())
      getattr(w, 'import')(kbkg_gg, ROOT.RooFit.RecycleConflictNodes())
      getattr(w, 'import')(mu, ROOT.RooFit.RecycleConflictNodes())
      getattr(w, 'import')(muV, ROOT.RooFit.RecycleConflictNodes())
      getattr(w, 'import')(muF, ROOT.RooFit.RecycleConflictNodes())

      ggZZpdf.SetNameTitle("ggzz", "ggzz")
      getattr(w, 'import')(ggZZpdf, ROOT.RooFit.RecycleConflictNodes())
      VBFpdf.SetNameTitle("vbf_offshell", "vbf_offshell")
      getattr(w, 'import')(VBFpdf, ROOT.RooFit.RecycleConflictNodes())
      bkg_zjets.SetNameTitle("bkg_zjets", "bkg_zjets")
      getattr(w, 'import')(bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
      bkg_qqzz.SetNameTitle("bkg_qqzz", "bkg_qqzz")
      getattr(w, 'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())

      getattr(w, 'import')(ggZZpdf_norm, ROOT.RooFit.RecycleConflictNodes())
      getattr(w, 'import')(VBFpdf_norm, ROOT.RooFit.RecycleConflictNodes())
      bkg_qqzz_norm.SetNameTitle("bkg_qqzz_norm", "bkg_qqzz_norm")
      getattr(w, 'import')(bkg_qqzz_norm, ROOT.RooFit.RecycleConflictNodes())
      if useDjet!=0:
      bkg_zjets_norm.SetNameTitle("bkg_zjets_norm", "bkg_zjets_norm")
      getattr(w, 'import')(bkg_zjets_norm, ROOT.RooFit.RecycleConflictNodes())

      w.writeToFile(name_ShapeWS)

      # --------------------------- DATACARDS -------------------------- ##

      systematics.setSystematics(
      bkgRate_qqzz_Shape, totalRate_ggzz_Shape, bkgRate_zjets_Shape)
      systematics_forXSxBR.setSystematics(
      bkgRate_qqzz_Shape, totalRate_ggzz_Shape, bkgRate_zjets_Shape)

      # If the channel is not declared in inputs, set rate = 0
      if not self.ggH_chan:
      sigRate_ggH_Shape = 0
      if not self.qqH_chan:
      sigRate_VBF_Shape = 0
      if not self.WH_chan:
      sigRate_WH_Shape = 0
      if not self.ZH_chan:
      sigRate_ZH_Shape = 0
      if not self.ttH_chan:
      sigRate_ttH_Shape = 0

      if not self.qqZZ_chan:
      bkgRate_qqzz_Shape = 0
      if not self.ggZZ_chan:
      totalRate_ggzz_Shape = 0
      if not self.zjets_chan:
      bkgRate_zjets_Shape = 0
      if not self.VBF_offshell_chan:
      totalRate_vbf_Shape = 0

      rates = {}
      rates['ggH'] = sigRate_ggH_Shape
      rates['qqH'] = sigRate_VBF_Shape
      rates['WH'] = sigRate_WH_Shape
      rates['ZH'] = sigRate_ZH_Shape
      rates['ttH'] = sigRate_ttH_Shape

      rates['qqZZ'] = bkgRate_qqzz_Shape
      rates['ggZZ'] = 1
      rates['VBF_offshell'] = 1
      rates['ggZZ_signal'] = rate_signal_ggzz_Shape
      rates['ggZZbkg'] = rate_bkg_ggzz_Shape
      rates['ggZZ_interf'] = rate_interf_ggzz_Shape
      rates['VBF_offshell_signal'] = rate_signal_vbf_Shape
      rates['VBF_offshell_bkg'] = rate_bkg_vbf_Shape
      rates['VBF_offshell_interf'] = rate_interf_vbf_Shape
      rates['zjets'] = bkgRate_zjets_Shape
      rates['ggZZ_tot'] = totalRate_ggzz_Shape
      rates['VBF_offshell_tot'] = totalRate_vbf_Shape
      rates['ttbar'] = 0
      rates['zbb'] = 0

      # Write Datacards
      fo = open(name_Shape, "wb")
      self.WriteDatacard(
      fo, self.theInputs, name_ShapeWS2, rates, data_obs_red.numEntries())

      systematics.WriteSystematics(fo, self.theInputs)
      systematics.WriteShapeSystematics(fo, self.theInputs)

      fo.close()

   def WriteDatacard(self, file, self.theInputs, nameWS, theRates, obsEvents):

      numberSig = self.numberOfSigChan(self.theInputs)
      numberBg = self.numberOfBgChan(self.theInputs)

      file.write("imax 1\n")
      file.write("jmax {0}\n".format(numberSig + numberBg - 1))
      file.write("kmax *\n")

      file.write("------------\n")
      file.write(
      "shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
      file.write("------------\n")

      file.write("bin a{0} \n".format(self.channel))
      file.write("observation {0} \n".format(obsEvents))

      file.write("------------\n")
      file.write(
      "## mass window [{0},{1}] \n".format(self.low_M, self.high_M))
      file.write("## signal,bkg,interf,tot rates [{0:.4f}, {1:.4f}, {2:.4f}, {3:.4f}] \n".format(
      theRates["ggZZ_signal"], theRates["ggZZbkg"], theRates["ggZZ_interf"], theRates["ggZZ_tot"]))
      file.write("## vbfsig,vbfbkg,vbfinterf,vbftot rates [{0:.4f}, {1:.4f}, {2:.4f}, {3:.4f}] \n".format(
      theRates["VBF_offshell_signal"], theRates["VBF_offshell_bkg"], theRates["VBF_offshell_interf"], theRates["VBF_offshell_tot"]))
      file.write("bin ")

      # channelList=['ggZZ_signal','ggZZ_interf','ggZZbkg','qqZZ','zjets']
      channelList = ['ggZZ', 'VBF_offshell', 'qqZZ', 'zjets']

      # channelName=['ggsignalzz','gginterfzz','ggZZbkg','bkg_qqzz','bkg_zjets']
      channelName = ['ggzz', 'vbf_offshell', 'bkg_qqzz', 'bkg_zjets']

      for chan in channelList:
      if self.theInputs[chan]:
      file.write("a{0} ".format(self.channel))
      file.write("\n")

      file.write("process ")

      i = 0

      for chan in channelList:
      if self.theInputs[chan]:
      file.write("{0} ".format(channelName[i]))
      i += 1

      file.write("\n")

      processLine = "process "

      for x in range(-numberSig + 1, 1):
      processLine += "{0} ".format(x)

      for y in range(1, numberBg + 1):
      processLine += "{0} ".format(y)

      file.write(processLine)
      file.write("\n")

      file.write("rate ")
      for chan in channelList:
      if self.theInputs[chan]:
      file.write("{0:.4f} ".format(theRates[chan]))
      file.write("\n")
      file.write("------------\n")

   def numberOfSigChan(self, inputs):

      counter = 0

      if inputs['ggZZ']:
      counter += 1
      if inputs['ggZZ_signal']:
      counter += 1
      if inputs['ggZZ_interf']:
      counter += 1
      if inputs['VBF_offshell']:
      counter += 1

      return counter

   def numberOfBgChan(self, inputs):

      counter = 0

      if inputs['qqZZ']:
      counter += 1
      if inputs['ggZZbkg']:
      counter += 1
      if inputs['zjets']:
      counter += 1
      if inputs['ttbar']:
      counter += 1
      if inputs['zbb']:
      counter += 1

      return counter
