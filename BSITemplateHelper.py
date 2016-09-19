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

class BSITemplateHelper:

   def __init__(self, options, theMaker, theCategorizer, strBSIType, templateFileName, iCat, systName):
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.channel = theMaker.channel
      self.workspace = theMaker.workspace

      # RooRealVars from the datacard maker class
      self.muF = theMaker.muF # Could itself be a RooFormulaVar (e.g. Off-shell: muF = R*RF*x. On-shell: muF = R*RF)
      self.muV = theMaker.muV # Could itself be a RooFormulaVar (e.g. Off-shell: muV = R*RV*x. On-shell: muV = R*RV)
      self.kbkg_gg = theMaker.kbkg_gg
      self.kbkg_VBF = theMaker.kbkg_VBF
      self.fai1 = theMaker.fai1
      self.phiai1 = theMaker.phiai1
      self.phia1_gg = theMaker.phia1_gg # Could itself be a RooFormulaVar (e.g. phia1_gg = phia1+phi_SB_gg)
      self.phia1_VBF = theMaker.phia1_VBF # Could itself be a RooFormulaVar (e.g. phia1_VBF = phia1+phi_SB_VBF/2)

      self.varm4l = theMaker.CMS_zz4l_widthMass
      self.varKD = theMaker.CMS_zz4l_widthKD
      self.varKD2 = theMaker.CMS_zz4l_widthKD2

      self.low_M = options.mLow
      self.high_M = options.mHigh
      self.anomCoupl=options.anomCouplIndex
      self.isBkgSigOnly = options.isBkgSigOnly
      self.templateDir = options.templateDir
      self.dimensions = optionds.dimensions # Number of template dimensions>0
      self.ProjDim = options.ProjDim # The projected variable, -1 means do not project

      self.iCatScheme = theCategorizer.iCatScheme
      self.catNameList = theCategorizer.catNameList
      self.nCategories = theCategorizer.nCategories
      self.iCat = iCat
      if self.iCat>=self.nCategories:
         sys.exit("self.iCat={} >= self.nCategories={}!".format(self.iCat,self.nCategories))

      self.templateFileName = templateFileName
      self.systName = systName
      self.templateSuffix = "{0}_{1}_{2:.0f}_{3:.0f}TeV".format(self.systName,self.catNameList[self.iCat],self.channel,self.sqrts)


      # To be reset later
      self.nbinsx=(self.high_M - self.low_M) / 20
      self.nbinsy=30
      self.nbinsz=30
      self.templateXLow=self.low_M
      self.templateYLow=0
      self.templateZLow=0
      self.templateXHigh=self.high_M
      self.templateYHigh=1
      self.templateZhigh=1

      # Template file
      self.templateFile = None

      self.strBSIType = strBSIType # gg-like or VVH-like couplings structure
      self.isGGVVLikeCouplings = strBSIType.lower().startswith("gg")
   # Extended template lists
   # Bare SM
      # ggF
      self.gg_T_1 = None
      self.gg_T_2 = None
      self.gg_T_4_Re = None
      self.gg_T_4_Im = None
      # VBF
      self.VBF_T_1 = None
      self.VBF_T_2 = None
      self.VBF_T_4_Re = None
      self.VBF_T_4_Im = None
   # Signal ai**1 x a1**(2/4-1) real and imaginary parts
      # ggF
      self.gg_T_1_AC_1_Re = None
      self.gg_T_1_AC_1_Im = None
      # VBF
      self.VBF_T_1_AC_1_Re = None
      self.VBF_T_1_AC_1_Im = None
   # Signal ai**2 x a1**(2/4-2) real and imaginary parts
      # ggF
      self.gg_T_1_AC_2_Re = None
      # VBF
      self.VBF_T_1_AC_2_Re = None
      self.VBF_T_1_AC_2_Im = None
      self.VBF_T_1_AC_2_PosDef = None
   # Signal ai**3 x a1**1 real and imaginary parts
      # No ggF
      # VBF
      self.VBF_T_1_AC_3_Re = None
      self.VBF_T_1_AC_3_Im = None
   # Signal ai**4 x a1**0 real and imaginary parts
      # No ggF
      # VBF
      self.VBF_T_1_AC_4 = None

   # Interference ai**1 x a1**(1/2-1) real and imaginary parts
      # ggF
      self.gg_T_4_AC_1_Re = None
      self.gg_T_4_AC_1_Im = None
      # VBF
      self.VBF_T_4_AC_1_Re = None
      self.VBF_T_4_AC_1_Im = None
   # Interference ai**2 x a1**(2-2) real and imaginary parts
      # No ggF
      # VBF
      self.VBF_T_4_AC_2_Re = None
      self.VBF_T_4_AC_2_Im = None

   # BSI+AC FORMULAE
      # Bare formula strings
      self.ggSigFormula_list = []
      self.ggInterfFormula_list = []
      self.VBFSigFormula_list = []
      self.VBFInterfFormula_list = []
      # RooFormulaVars
      self.ggSigRFV_list = []
      self.ggInterfRFV_list = []
      self.VBFSigRFV_list = []
      self.VBFInterfRFV_list = []

   # PDF construction
      # Lists of template arguments
      self.ggSigFunctions_Args = []
      self.VBFSigFunctions_Args = []
      self.ggInterfFunctions_Args = []
      self.VBFInterfFunctions_Args = []
      self.ggHistFunc_Arg = None
      self.VBFHistFunc_Arg = None
      # p.d.f. lists
      self.ggPdf = None
      self.VBFPdf = None

   # Rate construction
      # Lists of rate arguments
      self.ggSigRates_RooFormulaVar = None
      self.ggInterfRates_RooFormulaVar = None
      self.ggBkgRates_RooFormulaVar = None
      self.VBFSigRates_RooFormulaVar = None
      self.VBFInterfRates_RooFormulaVar = None
      self.VBFBkgRates_RooFormulaVar = None
      # Rate lists
      self.ggTotalRate = None
      self.VBFTotalRate = None


      # Need to import these only once FIXME: to be moved to the maker
      self.workspace.importClassCode(AsymPow.Class(),True)
      self.workspace.importClassCode(AsymQuad.Class(),True)
      self.workspace.importClassCode(RooqqZZPdf_v2.Class(), True)
      self.workspace.importClassCode(RooFormulaVar.Class(), True)
      self.workspace.importClassCode(RooRealFlooredSumPdf.Class(),True)
      self.workspace.importClassCode(VerticalInterpPdf.Class(),True)


# Close the template files
   def close(self):
      self.templateFile.Close()


# Get shapes for each category
   def getTemplates(self,processName=None,templatePrefix="T_2D"):
      self.templateFile = ROOT.TFile.Open(self.templateFileName, "read")

      self.templatePrefix = templatePrefix
      if processName is not None:
         self.processName = processName
         self.templatePrefix = "{}_{}".format(self.templatePrefix,processName)
      elif not self.isGGVVLikeCouplings:
         self.processName = "VBF"
         self.templatePrefix = "{}_{}".format(self.templatePrefix,processName)
      else:
         self.processName = "gg"

      if self.isGGVVLikeCouplings:
         self.getTemplates_ggVVLike()
      else:
         self.getTemplates_vvVVLike()


   def getTemplates_ggVVLike(self):
#---------- SM SIGNAL AND BACKGROUND TEMPLATES -------------
# Bare SM
      self.gg_T_1 =
         ExtendedTemplate(
               self.templateFile.Get("{}_1".format(self.templatePrefix)).Clone("{}_1_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions, self.ProjDim,
               self.varm4l, self.varKD, self.varKD2
            )
      self.gg_T_2 =
         ExtendedTemplate(
               self.templateFile.Get("{}_2".format(self.templatePrefix)).Clone("{}_2_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions, self.ProjDim,
               self.varm4l, self.varKD, self.varKD2
            )
      if(not(self.isBkgSigOnly)):
         self.gg_T_4_Re =
            ExtendedTemplate(
                  self.templateFile.Get("{}_4_Re".format(self.templatePrefix)).Clone("{}_4_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
         if (self.anomCoupl==1):
            self.gg_T_4_Im =
               ExtendedTemplate(
                     self.templateFile.Get("{}_4_Im".format(self.templatePrefix)).Clone("{}_4_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
      # Special case: Get template properties from the gg bkg template
      if icat == 0:
         self.nbinsx=self.gg_T_2.theTemplate.GetNbinsX()
         self.nbinsy=self.gg_T_2.theTemplate.GetNbinsY()
         self.nbinsz=self.gg_T_2.theTemplate.GetNbinsZ()
         self.templateXLow=self.gg_T_2.theTemplate.GetXaxis().GetBinLowEdge(1)
         self.templateYLow=self.gg_T_2.theTemplate.GetYaxis().GetBinLowEdge(1)
         self.templateZLow=self.gg_T_2.theTemplate.GetZaxis().GetBinLowEdge(1)
         self.templateXHigh=self.gg_T_2.theTemplate.GetXaxis().GetBinUpEdge(self.nbinsx)
         self.templateYHigh=self.gg_T_2.theTemplate.GetYaxis().GetBinUpEdge(self.nbinsy)
         self.templateZhigh=self.gg_T_2.theTemplate.GetZaxis().GetBinUpEdge(self.nbinsz)
         self.blankTemplate = self.gg_T_2.theTemplate.Clone("blankTemplate")
         self.blankTemplate.Reset("M")

      if self.anomCoupl != 0:
#-----------------------------------------------------------------------#
#                        SIGNAL AC TERMS
#-----------------------------------------------------------------------#
# Signal ai**1 x a1**(2/4-1) real and imaginary parts
         self.gg_T_1_AC_1_Re =
            ExtendedTemplate(
                  self.templateFile.Get("{}_1_AC_1_Re".format(self.templatePrefix)).Clone("{}_1_AC_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
         if self.anomCoupl == 1:
            self.gg_T_1_AC_1_Im =
               ExtendedTemplate(
                     self.templateFile.Get("{}_1_AC_1_Im".format(self.templatePrefix)).Clone("{}_1_AC_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )

# Signal ai**2 x a1**(2/4-2) real and imaginary parts
         self.gg_T_1_AC_2_Re =
            ExtendedTemplate(
                  self.templateFile.Get("{}_1_AC_2_Re".format(self.templatePrefix)).Clone("{}_1_AC_2_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
         if self.anomCoupl == 1:
            continue # No ggF term
         elif self.anomCoupl == 2: # if self.anomCoupl == 2, PosDef = PosDef+Re
            continue # No ggF term

# Signal ai**3 x a1**1 real and imaginary parts
   # No ggF term

# Signal ai**4 x a1**0 real and imaginary parts
   # No ggF term

#-----------------------------------------------------------------------#
#                       INTERFERENCE AC TERMS
#-----------------------------------------------------------------------#
         if(not(self.isBkgSigOnly)):
# Interference ai**1 x a1**(1/2-1) real and imaginary parts
            self.gg_T_4_AC_1_Re =
               ExtendedTemplate(
                     self.templateFile.Get("{}_4_AC_1_Re".format(self.templatePrefix)).Clone("{}_4_AC_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
            if self.anomCoupl == 1:
               self.gg_T_4_AC_1_Im =
                  ExtendedTemplate(
                        self.templateFile.Get("{}_4_AC_1_Im".format(self.templatePrefix)).Clone("{}_4_AC_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )

# Interference ai**2 x a1**(2-2) real and imaginary parts
   # No ggF term

# FORMULAE
      # In signals and interferences, @0==muF/V; additionally in interferences, @1==kbkg_gg/VBF
      if self.anomCoupl == 1: # Full parameterization with fai1=[-1, 1], and phases phiai1 and phia1_gg, phia1_VBF (phia1 are different since the bkg phase could be different)
         # 0-3 are templates, @1==fai1, @2==phiai1
         self.ggSigFormula_list.append("@0*(1-abs(@1))")
         self.ggSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))*cos(@2)")
         self.ggSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))*sin(@2)")
         self.ggSigFormula_list.append("@0*abs(@1)")
         if(not(self.isBkgSigOnly)):
            # 0-3 are templates, @2==fai1, @3==phiai1, @4==phia1 (gg)
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sqrt(1-abs(@2))*cos(@4)")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sqrt(1-abs(@2))*sin(@4)")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))*cos(@3+@4)")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))*sin(@3+@4)")
      elif self.anomCoupl == 2: # No phases, just fai1=[-1, 1]
         # 0-2 are templates, @1==fai1
         self.ggSigFormula_list.append("@0*(1-abs(@1))")
         self.ggSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))")
         self.ggSigFormula_list.append("@0*abs(@1)")
         if(not(self.isBkgSigOnly)):
            # 0-1 are templates, @2==fai1
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sqrt(1-abs(@2))")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))")
      else: # No ai1 dependence
         self.ggSigFormula_list.append("@0")
         if(not(self.isBkgSigOnly)):
            self.ggInterfFormula_list.append("sqrt(@0*@1)")

      for irfv in range(0,len(self.ggSigFormula_list)):
         rfvname = "{0}Sig_AC_{1:.0f}_Coef".format(self.processName,irfv)
         rfvargs = ROOT.RooArgList()
         rfvargs.add(self.muF)
         if self.anomCoupl == 1:
            rfvargs.add(self.fai1)
            rfvargs.add(self.phiai1)
         elif self.anomCoupl == 2:
            rfvargs.add(self.fai1)
         if rfvargs.getSize()>0:
            seg_rfv = ROOT.RooFormulaVar( rfvname , self.ggSigFormula_list[irfv] , rfvargs )
            self.ggSigRFV_list.append(seg_rfv)
         else:
            seg_rfv = ROOT.RooRealVar( rfvname , rfvname , 1.0 )
            self.ggSigRFV_list.append(seg_rfv)

      for irfv in range(0,len(self.ggInterfFormula_list)):
         rfvname = "{0}Interf_AC_{1:.0f}_Coef".format(self.processName,irfv)
         rfvargs = ROOT.RooArgList()
         rfvargs.add(self.muF)
         rfvargs.add(self.kbkg_gg)
         if self.anomCoupl == 1:
            rfvargs.add(self.fai1)
            rfvargs.add(self.phiai1)
            rfvargs.add(self.phia1_gg)
         elif self.anomCoupl == 2:
            rfvargs.add(self.fai1)
         if rfvargs.getSize()>0:
            seg_rfv = ROOT.RooFormulaVar( rfvname , self.ggInterfFormula_list[irfv] , rfvargs )
            self.ggInterfRFV_list.append(seg_rfv)
         else:
            seg_rfv = ROOT.RooRealVar( rfvname , rfvname , 1.0 )
            self.ggInterfRFV_list.append(seg_rfv)

# Lists of template arguments
      self.ggSigFunctions_Args.append(self.gg_T_1)
      if self.anomCoupl == 1:
         self.ggSigFunctions_Args.append(self.gg_T_1_AC_1_Re)
         self.ggSigFunctions_Args.append(self.gg_T_1_AC_1_Im)
         self.ggSigFunctions_Args.append(self.gg_T_1_AC_2_Re)
      elif self.anomCoupl == 2:
         self.ggSigFunctions_Args.append(self.gg_T_1_AC_1_Re)
         self.ggSigFunctions_Args.append(self.gg_T_1_AC_2_Re)
      if len(self.ggSigFunctions_Args)!=len(self.ggSigRFV_list):
         sys.exit("Number of {0}Sig templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.ggSigFunctions_Args),len(self.ggSigRFV_list)))

      if(not(self.isBkgSigOnly)):
         self.ggInterfFunctions_Args.append(self.gg_T_4_Re)
         if self.anomCoupl == 1:
            self.ggInterfFunctions_Args.append(self.gg_T_4_Im)
            self.ggInterfFunctions_Args.append(self.gg_T_4_AC_1_Re)
            self.ggInterfFunctions_Args.append(self.gg_T_4_AC_1_Im)
         elif self.anomCoupl == 2:
            self.ggInterfFunctions_Args.append(self.gg_T_4_AC_1_Re)
         if len(self.ggInterfFunctions_Args)!=len(self.ggInterfRFV_list):
            sys.exit("Number of {0}Interf templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.ggInterfFunctions_Args),len(self.ggInterfRFV_list)))

      self.ggHistFunc_Arg = ROOT.RooArgList()
      for var in self.ggSigFunctions_Args:
         self.ggHistFunc_Arg.add(var.theHistFunc)
      for var in self.ggInterfFunctions_Args:
         self.ggHistFunc_Arg.add(var.theHistFunc)
      self.ggHistFunc_Arg.add(self.gg_T_2.theHistFunc)

# Construct the p.d.f.'s
      rfvargs = ROOT.RooArgList()
      for var in self.ggSigRFV_list:
         rfvargs.add(var)
      for var in self.ggInterfRFV_list:
         rfvargs.add(var)
      rfvargs.add(self.kbkg_gg)
      PdfName = "{}Pdf_{}".format(self.processName,self.templateSuffix))
      self.ggPdf = ROOT.RooRealSumPdf(
         PdfName, PdfName,
         self.ggHistFunc_Arg,rfvargs
      )

# Lists of rate FormulaVars
# Each signal, bkg and interf is constructed separately to be able to count their contributions in the end
      rfvargs = ROOT.RooArgList()
      strformula = ""
      for ivar in range(0,len(self.ggSigFunctions_Args)):
         rfvargs.add(self.ggSigRFV_list[ivar])
         rfvargs.add(self.ggSigFunctions_Args[ivar].theRate)
         if ivar==0:
            strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
         else:
            strformula = "{2} + {0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1,strformula)
      rfvname = "{}SigRate_{}".format(self.processName,self.templateSuffix)
      self.ggSigRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      if len(self.ggInterfFunctions_Args)>0:
         rfvargs = ROOT.RooArgList()
         strformula = ""
         for ivar in range(0,len(self.ggInterfFunctions_Args)):
            rfvargs.add(self.ggInterfRFV_list[ivar])
            rfvargs.add(self.ggInterfFunctions_Args[ivar].theRate)
            if ivar==0:
               strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
            else:
               strformula = "{2} + {0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1,strformula)
         rfvname = "{}InterfRate_{}".format(self.processName,self.templateSuffix)
         self.ggInterfRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      rfvargs = ROOT.RooArgList()
      strformula = "@0*@1"
      rfvargs.add(self.kbkg_gg)
      rfvargs.add(self.gg_T_2.theRate)
      rfvname = "{}BkgRate_{}".format(self.processName,self.templateSuffix)
      self.ggBkgRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

# Construct total rates
      rfvargs = ROOT.RooArgList()
      strformula = "@0+@1"
      rfvargs.add(self.ggSigRates_RooFormulaVar)
      rfvargs.add(self.ggBkgRates_RooFormulaVar)
      if(not(self.isBkgSigOnly)):
         strformula = "@0+@1+@2"
         rfvargs.add(self.ggInterfRates_RooFormulaVar)
      rfvname = "{}TotalRate_{}".format(self.processName,self.templateSuffix)
      self.ggTotalRate = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )


   def getTemplates_vvVVLike(self):
#---------- SM SIGNAL AND BACKGROUND TEMPLATES -------------
# Bare SM
      self.VBF_T_1 =
         ExtendedTemplate(
               self.templateFile.Get("{}_1".format(self.templatePrefix)).Clone("{}_1_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions, self.ProjDim,
               self.varm4l, self.varKD, self.varKD2
            )
      self.VBF_T_2 =
         ExtendedTemplate(
               self.templateFile.Get("{}_2".format(self.templatePrefix)).Clone("{}_2_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions, self.ProjDim,
               self.varm4l, self.varKD, self.varKD2
            )
      if(not(self.isBkgSigOnly)):
         self.VBF_T_4_Re =
            ExtendedTemplate(
                  self.templateFile.Get("{}_4_Re".format(self.templatePrefix)).Clone("{}_4_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
         if (self.anomCoupl==1):
            self.VBF_T_4_Im =
               ExtendedTemplate(
                     self.templateFile.Get("{}_4_Im".format(self.templatePrefix)).Clone("{}_4_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )

      if self.anomCoupl != 0:
#-----------------------------------------------------------------------#
#                        SIGNAL AC TERMS
#-----------------------------------------------------------------------#
# Signal ai**1 x a1**(2/4-1) real and imaginary parts
         self.VBF_T_1_AC_1_Re =
            ExtendedTemplate(
                  self.templateFile.Get("{}_1_AC_1_Re".format(self.templatePrefix)).Clone("{}_1_AC_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
         if self.anomCoupl == 1:
            self.VBF_T_1_AC_1_Im =
               ExtendedTemplate(
                     self.templateFile.Get("{}_1_AC_1_Im".format(self.templatePrefix)).Clone("{}_1_AC_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )

# Signal ai**2 x a1**(2/4-2) real and imaginary parts
         if self.anomCoupl == 1:
            self.VBF_T_1_AC_2_PosDef =
               ExtendedTemplate(
                     self.templateFile.Get("{}_1_AC_2_PosDef".format(self.templatePrefix)).Clone("{}_1_AC_2_PosDef_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
            self.VBF_T_1_AC_2_Re =
               ExtendedTemplate(
                     self.templateFile.Get("{}_1_AC_2_Re".format(self.templatePrefix)).Clone("{}_1_AC_2_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
            self.VBF_T_1_AC_2_Im =
               ExtendedTemplate(
                     self.templateFile.Get("{}_1_AC_2_Im".format(self.templatePrefix)).Clone("{}_1_AC_2_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
         elif self.anomCoupl == 2: # if self.anomCoupl == 2, PosDef = PosDef+Re
            tmpTpl = self.templateFile.Get("{}_1_AC_2_PosDef".format(self.templatePrefix)).Clone("{}_1_AC_2_PosDef_{}".format(self.templatePrefix,self.templateSuffix))
            tmpTpl.Add(self.templateFile.Get("{}_1_AC_2_Re".format(self.templatePrefix)).Clone("{}_1_AC_2_Re_{}".format(self.templatePrefix,self.templateSuffix)))
            self.VBF_T_1_AC_2_PosDef =
               ExtendedTemplate(
                     tmpTpl,
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )

# Signal ai**3 x a1**1 real and imaginary parts
         self.VBF_T_1_AC_3_Re =
            ExtendedTemplate(
                  self.templateFile.Get("{}_1_AC_3_Re".format(self.templatePrefix)).Clone("{}_1_AC_3_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )

         if self.anomCoupl == 1:
            self.VBF_T_1_AC_3_Im =
               ExtendedTemplate(
                     self.templateFile.Get("{}_1_AC_3_Im".format(self.templatePrefix)).Clone("{}_1_AC_3_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )

# Signal ai**4 x a1**0 real and imaginary parts
         self.VBF_T_1_AC_4 =
            ExtendedTemplate(
                  self.templateFile.Get("{}_1_AC_4".format(self.templatePrefix)).Clone("{}_1_AC_4_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )

#-----------------------------------------------------------------------#
#                       INTERFERENCE AC TERMS
#-----------------------------------------------------------------------#
         if(not(self.isBkgSigOnly)):
# Interference ai**1 x a1**(1/2-1) real and imaginary parts
            self.VBF_T_4_AC_1_Re =
               ExtendedTemplate(
                     self.templateFile.Get("{}_4_AC_1_Re".format(self.templatePrefix)).Clone("{}_4_AC_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
            if self.anomCoupl == 1:
               self.VBF_T_4_AC_1_Im =
                  ExtendedTemplate(
                        self.templateFile.Get("{}_4_AC_1_Im".format(self.templatePrefix)).Clone("{}_4_AC_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )

# Interference ai**2 x a1**(2-2) real and imaginary parts
            self.VBF_T_4_AC_2_Re =
               ExtendedTemplate(
                     self.templateFile.Get("{}_4_AC_2_Re".format(self.templatePrefix)).Clone("{}_4_AC_2_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
            if self.anomCoupl == 1:
               self.VBF_T_4_AC_2_Im =
                  ExtendedTemplate(
                        self.templateFile.Get("{}_4_AC_2_Im".format(self.templatePrefix)).Clone("{}_4_AC_2_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )

# FORMULAE
      # In signals and interferences, @0==muV; additionally in interferences, @1==kbkg_VBF
      if self.anomCoupl == 1: # Full parameterization with fai1=[-1, 1], and phases phiai1 and phia1_gg, phia1_VBF (phia1 are different since the bkg phase could be different)
         # 0-8 are templates, @1==fai1, @2==phiai1
         self.VBFSigFormula_list.append("@0*pow((1-abs(@1)),2)")
         self.VBFSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1))*pow(sqrt(1-abs(@1)),3)*cos(@2)")
         self.VBFSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1))*pow(sqrt(1-abs(@1)),3)*sin(@2)")
         self.VBFSigFormula_list.append("@0*abs(@1)*(1-abs(@1))")
         self.VBFSigFormula_list.append("@0*abs(@1)*(1-abs(@1))*cos(2*@2)")
         self.VBFSigFormula_list.append("@0*abs(@1)*(1-abs(@1))*sin(2*@2)")
         self.VBFSigFormula_list.append("@0*sign(@1)*pow(sqrt(abs(@1)),3)*sqrt(1-abs(@1))*cos(@2)")
         self.VBFSigFormula_list.append("@0*sign(@1)*pow(sqrt(abs(@1)),3)*sqrt(1-abs(@1))*sin(@2)")
         self.VBFSigFormula_list.append("@0*pow(@1,2)")
         if(not(self.isBkgSigOnly)):
            # 0-5 are templates, @2==fai1, @3==phiai1, @4==phia1 (VBF)
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*(1-abs(@2))*cos(2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*(1-abs(@2))*sin(2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2)*(1-abs(@2)))*cos(@3+2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2)*(1-abs(@2)))*sin(@3+2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*abs(@2)*cos(2*(@3+@4))")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*abs(@2)*sin(2*(@3+@4))")
      elif self.anomCoupl == 2: # No phases, just fai1=[-1, 1]
         # 0-4 are templates, @1==fai1
         self.VBFSigFormula_list.append("@0*pow((1-abs(@1)),2)")
         self.VBFSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1))*pow(sqrt(1-abs(@1)),3)")
         self.VBFSigFormula_list.append("@0*abs(@1)*(1-abs(@1))")
         self.VBFSigFormula_list.append("@0*sign(@1)*pow(sqrt(abs(@1)),3)*sqrt(1-abs(@1))")
         self.VBFSigFormula_list.append("@0*pow(@1,2)")
         if(not(self.isBkgSigOnly)):
            # 0-2 are templates, @2==fai1
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*(1-abs(@2))")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2)*(1-abs(@2)))")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*abs(@2)")
      else: # No ai1 dependence
         self.VBFSigFormula_list.append("@0")
         if(not(self.isBkgSigOnly)):
            self.VBFInterfFormula_list.append("sqrt(@0*@1)")

      for irfv in range(0,len(self.VBFSigFormula_list)):
         rfvname = "{0}Sig_AC_{1:.0f}_Coef".format(self.processName,irfv)
         rfvargs = ROOT.RooArgList()
         rfvargs.add(self.muV)
         if self.anomCoupl == 1:
            rfvargs.add(self.fai1)
            rfvargs.add(self.phiai1)
         elif self.anomCoupl == 2:
            rfvargs.add(self.fai1)
         if rfvargs.getSize()>0:
            seg_rfv = ROOT.RooFormulaVar( rfvname , self.VBFSigFormula_list[irfv] , rfvargs )
            self.VBFSigRFV_list.append(seg_rfv)
         else:
            seg_rfv = ROOT.RooRealVar( rfvname , rfvname , 1.0 )
            self.VBFSigRFV_list.append(seg_rfv)

      for irfv in range(0,len(self.VBFInterfFormula_list)):
         rfvname = "{0}Interf_AC_{1:.0f}_Coef".format(self.processName,irfv)
         rfvargs = ROOT.RooArgList()
         rfvargs.add(self.muV)
         rfvargs.add(self.kbkg_VBF)
         if self.anomCoupl == 1:
            rfvargs.add(self.fai1)
            rfvargs.add(self.phiai1)
            rfvargs.add(self.phia1_VBF)
         elif self.anomCoupl == 2:
            rfvargs.add(self.fai1)
         if rfvargs.getSize()>0:
            seg_rfv = ROOT.RooFormulaVar( rfvname , self.VBFInterfFormula_list[irfv] , rfvargs )
            self.VBFInterfRFV_list.append(seg_rfv)
         else:
            seg_rfv = ROOT.RooRealVar( rfvname , rfvname , 1.0 )
            self.VBFInterfRFV_list.append(seg_rfv)

# Lists of template arguments
      self.VBFSigFunctions_Args.append(self.VBF_T_1)
      if self.anomCoupl == 1:
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_1_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_1_Im)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_2_PosDef)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_2_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_2_Im)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_3_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_3_Im)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_4)
      elif self.anomCoupl == 2:
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_1_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_2_PosDef)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_3_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_1_AC_4)
      if len(self.VBFSigFunctions_Args)!=len(self.VBFSigRFV_list):
         sys.exit("Number of {0}Sig templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.VBFSigFunctions_Args),len(self.VBFSigRFV_list)))

      if(not(self.isBkgSigOnly)):
         self.VBFInterfFunctions_Args.append(self.VBF_T_4_Re)
         if self.anomCoupl == 1:
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_Im)
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_AC_1_Re)
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_AC_1_Im)
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_AC_2_Re)
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_AC_2_Im)
         elif self.anomCoupl == 2:
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_AC_1_Re)
            self.VBFInterfFunctions_Args.append(self.VBF_T_4_AC_2_Re)
         if len(self.VBFInterfFunctions_Args)!=len(self.VBFInterfRFV_list):
            sys.exit("Number of {0}Interf templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.VBFInterfFunctions_Args),len(self.VBFInterfRFV_list)))

      self.VBFHistFunc_Arg = ROOT.RooArgList()
      for var in self.VBFSigFunctions_Args:
         self.VBFHistFunc_Arg.add(var.theHistFunc)
      for var in self.VBFInterfFunctions_Args:
         self.VBFHistFunc_Arg.add(var.theHistFunc)
      self.VBFHistFunc_Arg.add(self.VBF_T_2.theHistFunc)

# Construct the p.d.f.'s
      rfvargs = ROOT.RooArgList()
      for var in self.VBFSigRFV_list:
         rfvargs.add(var)
      for var in self.VBFInterfRFV_list:
         rfvargs.add(var)
      rfvargs.add(self.kbkg_VBF)
      PdfName = "{}Pdf_{}".format(self.processName,self.templateSuffix))
      self.VBFPdf = ROOT.RooRealSumPdf(
         PdfName, PdfName,
         self.VBFHistFunc_Arg,rfvargs
      )

# Lists of rate FormulaVars
# Each signal, bkg and interf is constructed separately to be able to count their contributions in the end
      rfvargs = ROOT.RooArgList()
      strformula = ""
      for ivar in range(0,len(self.VBFSigFunctions_Args)):
         rfvargs.add(self.VBFSigRFV_list[ivar])
         rfvargs.add(self.VBFSigFunctions_Args[ivar].theRate)
         if ivar==0:
            strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
         else:
            strformula = "{2} + {0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1,strformula)
      rfvname = "{}SigRate_{}".format(self.processName,self.templateSuffix)
      self.VBFSigRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      if len(self.VBFInterfFunctions_Args)>0:
         rfvargs = ROOT.RooArgList()
         strformula = ""
         for ivar in range(0,len(self.VBFInterfFunctions_Args)):
            rfvargs.add(self.VBFInterfRFV_list[ivar])
            rfvargs.add(self.VBFInterfFunctions_Args[ivar].theRate)
            if ivar==0:
               strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
            else:
               strformula = "{2} + {0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1,strformula)
         rfvname = "{}InterfRate_{}".format(self.processName,self.templateSuffix)
         self.VBFInterfRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      rfvargs = ROOT.RooArgList()
      strformula = "@0*@1"
      rfvargs.add(self.kbkg_VBF)
      rfvargs.add(self.VBF_T_2.theRate)
      rfvname = "{}BkgRate_{}".format(self.processName,self.templateSuffix)
      self.VBFBkgRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

# Construct total rates
      rfvargs = ROOT.RooArgList()
      strformula = "@0+@1"
      rfvargs.add(self.VBFSigRates_RooFormulaVar)
      rfvargs.add(self.VBFBkgRates_RooFormulaVar)
      if(not(self.isBkgSigOnly)):
         strformula = "@0+@1+@2"
         rfvargs.add(self.VBFInterfRates_RooFormulaVar)
      rfvname = "{}TotalRate_{}".format(self.processName,self.templateSuffix)
      self.VBFTotalRate = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )



