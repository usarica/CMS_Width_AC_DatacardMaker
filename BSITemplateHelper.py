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

   def __init__(self, options, theMaker, theCategorizer):
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.channel = theMaker.channel
      self.workspace = theMAker.workspace

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

      self.templateFileList = []

      # Need to import these only once FIXME: to be moved to the maker
      self.workspace.importClassCode(AsymPow.Class(),True)
      self.workspace.importClassCode(AsymQuad.Class(),True)
      self.workspace.importClassCode(RooqqZZPdf_v2.Class(), True)
      self.workspace.importClassCode(RooFormulaVar.Class(), True)
      self.workspace.importClassCode(RooRealFlooredSumPdf.Class(),True)
      self.workspace.importClassCode(VerticalInterpPdf.Class(),True)


# Get shapes for each category
   def getTemplates(self, templateFileAppendName, systName):
      templateNameMain = "{0}/{1}".format(self.templateDir, templateFileAppendName)
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


#---------- SM SIGNAL AND BACKGROUND TEMPLATES -------------

# Bare SM
   # Templates
      # ggF
      self.gg_T_1_list = []
      self.gg_T_2_list = []
      self.gg_T_4_Re_list = []
      self.gg_T_4_Im_list = []
      # VBF
      self.VBF_T_1_list = []
      self.VBF_T_2_list = []
      self.VBF_T_4_Re_list = []
      self.VBF_T_4_Im_list = []

      for icat in range(0,self.nCategories):
         self.gg_T_1_list.append(
            ExtendedTemplate(
                  self.templateFileList[icat].Get("T_2D_1").Clone("T_2D_gg_1_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
            )
         self.gg_T_2_list.append(
            ExtendedTemplate(
                  self.templateFileList[icat].Get("T_2D_2").Clone("T_2D_gg_2_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
            )
         if(not(self.isBkgSigOnly)):
            self.gg_T_4_Re_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_4_Re").Clone("T_2D_gg_4_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )

         self.VBF_T_1_list.append(
            ExtendedTemplate(
                  self.templateFileList[icat].Get("T_2D_VBF_1").Clone("T_2D_VBF_1_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
            )
         self.VBF_T_2_list.append(
            ExtendedTemplate(
                  self.templateFileList[icat].Get("T_2D_VBF_2").Clone("T_2D_VBF_2_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                  self.dimensions, self.ProjDim,
                  self.varm4l, self.varKD, self.varKD2
               )
            )
         if(not(self.isBkgSigOnly)):
            self.VBF_T_4_Re_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_VBF_4_Re").Clone("T_2D_VBF_4_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )

         if (self.anomCoupl==1 and not(self.isBkgSigOnly)):
            self.gg_T_4_Im_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_4_Im").Clone("T_2D_gg_4_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )
            self.VBF_T_4_Im_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_VBF_4_Im").Clone("T_2D_VBF_4_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )

         # Special case: Get template properties from the gg bkg template
         if icat == 0:
            self.nbinsx=self.gg_T_2_list[icat].theTemplate.GetNbinsX()
            self.nbinsy=self.gg_T_2_list[icat].theTemplate.GetNbinsY()
            self.nbinsz=self.gg_T_2_list[icat].theTemplate.GetNbinsZ()
            self.templateXLow=self.gg_T_2_list[icat].theTemplate.GetXaxis().GetBinLowEdge(1)
            self.templateYLow=self.gg_T_2_list[icat].theTemplate.GetYaxis().GetBinLowEdge(1)
            self.templateZLow=self.gg_T_2_list[icat].theTemplate.GetZaxis().GetBinLowEdge(1)
            self.templateXHigh=self.gg_T_2_list[icat].theTemplate.GetXaxis().GetBinUpEdge(self.nbinsx)
            self.templateYHigh=self.gg_T_2_list[icat].theTemplate.GetYaxis().GetBinUpEdge(self.nbinsy)
            self.templateZhigh=self.gg_T_2_list[icat].theTemplate.GetZaxis().GetBinUpEdge(self.nbinsz)
            self.blankTemplate = self.gg_T_2_list[icat].theTemplate.Clone("blankTemplate")
            self.blankTemplate.Reset("M")


#-----------------------------------------------------------------------#
#                        SIGNAL AC TERMS
#-----------------------------------------------------------------------#
# Signal ai**1 x a1**(2/4-1) real and imaginary parts
   # Templates
      # ggF
      self.gg_T_1_AC_1_Re_list = []
      self.gg_T_1_AC_1_Im_list = []
      # VBF
      self.VBF_T_1_AC_1_Re_list = []
      self.VBF_T_1_AC_1_Im_list = []

      if self.anomCoupl != 0:
         for icat in range(0,self.nCategories):
            self.gg_T_1_AC_1_Re_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_1_AC_1_Re").Clone("T_2D_gg_1_AC_1_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )
            self.VBF_T_1_AC_1_Re_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_VBF_1_AC_1_Re").Clone("T_2D_VBF_1_AC_1_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )

            if self.anomCoupl == 1:
               self.gg_T_1_AC_1_Im_list.append(
                  ExtendedTemplate(
                        self.templateFileList[icat].Get("T_2D_1_AC_1_Im").Clone("T_2D_gg_1_AC_1_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )
               self.VBF_T_1_AC_1_Im_list.append(
                  ExtendedTemplate(
                        self.templateFileList[icat].Get("T_2D_VBF_1_AC_1_Im").Clone("T_2D_VBF_1_AC_1_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )


# Signal ai**2 x a1**(2/4-2) real and imaginary parts
   # Templates
      # ggF
      self.gg_T_1_AC_2_Re_list = []
      # VBF
      self.VBF_T_1_AC_2_Re_list = []
      self.VBF_T_1_AC_2_Im_list = []
      self.VBF_T_1_AC_2_PosDef_list = []

      if self.anomCoupl != 0:
         for icat in range(0,self.nCategories):
            self.gg_T_1_AC_2_Re_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_1_AC_2_Re").Clone("T_2D_gg_1_AC_2_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )

            if self.anomCoupl == 1:
               self.VBF_T_1_AC_2_PosDef_list.append(
                  ExtendedTemplate(
                        self.templateFileList[icat].Get("T_2D_VBF_1_AC_2_PosDef").Clone("T_2D_VBF_1_AC_2_PosDef_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )

               self.VBF_T_1_AC_2_Re_list.append(
                  ExtendedTemplate(
                        self.templateFileList[icat].Get("T_2D_VBF_1_AC_2_Re").Clone("T_2D_VBF_1_AC_2_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )

               self.VBF_T_1_AC_2_Im_list.append(
                  ExtendedTemplate(
                        self.templateFileList[icat].Get("T_2D_VBF_1_AC_2_Im").Clone("T_2D_VBF_1_AC_2_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )
            elif self.anomCoupl == 2: # if self.anomCoupl == 2, PosDef = PosDef+Re
               tmpTpl = self.templateFileList[icat].Get("T_2D_VBF_1_AC_2_PosDef").Clone("T_2D_VBF_1_AC_2_PosDef_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts))
               tmpTpl.Add(self.templateFileList[icat].Get("T_2D_VBF_1_AC_2_Re").Clone("T_2D_VBF_1_AC_2_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)))
               self.VBF_T_1_AC_2_PosDef_list.append(
                  ExtendedTemplate(
                        tmpTpl,
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )


# Signal ai**3 x a1**1 real and imaginary parts
   # Templates
      # No ggF
      # VBF
      self.VBF_T_1_AC_3_Re_list = []
      self.VBF_T_1_AC_3_Im_list = []

      if self.anomCoupl != 0:
         for icat in range(0,self.nCategories):
            self.VBF_T_1_AC_3_Re_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_VBF_1_AC_3_Re").Clone("T_2D_VBF_1_AC_3_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )

            if self.anomCoupl == 1:
               self.VBF_T_1_AC_3_Im_list.append(
                  ExtendedTemplate(
                        self.templateFileList[icat].Get("T_2D_VBF_1_AC_3_Im").Clone("T_2D_VBF_1_AC_3_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                        self.dimensions, self.ProjDim,
                        self.varm4l, self.varKD, self.varKD2
                     )
                  )


# Signal ai**4 x a1**0 real and imaginary parts
   # Templates
      # No ggF
      # VBF
      self.VBF_T_1_AC_4_list = []

      if self.anomCoupl != 0:
         for icat in range(0,self.nCategories):
            self.VBF_T_1_AC_4_list.append(
               ExtendedTemplate(
                     self.templateFileList[icat].Get("T_2D_VBF_1_AC_4").Clone("T_2D_VBF_1_AC_4_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                     self.dimensions, self.ProjDim,
                     self.varm4l, self.varKD, self.varKD2
                  )
               )


#-----------------------------------------------------------------------#
#                       INTERFERENCE AC TERMS
#-----------------------------------------------------------------------#
      if(not(self.isBkgSigOnly)):
# Interference ai**1 x a1**(1/2-1) real and imaginary parts
         # Templates
            # ggF
            self.gg_T_4_AC_1_Re_list = []
            self.gg_T_4_AC_1_Im_list = []
            # VBF
            self.VBF_T_4_AC_1_Re_list = []
            self.VBF_T_4_AC_1_Im_list = []

            if self.anomCoupl != 0:
               for icat in range(0,self.nCategories):
                  self.gg_T_4_AC_1_Re_list.append(
                     ExtendedTemplate(
                           self.templateFileList[icat].Get("T_2D_4_AC_1_Re").Clone("T_2D_gg_4_AC_1_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                           self.dimensions, self.ProjDim,
                           self.varm4l, self.varKD, self.varKD2
                        )
                     )
                  self.VBF_T_4_AC_1_Re_list.append(
                     ExtendedTemplate(
                           self.templateFileList[icat].Get("T_2D_VBF_4_AC_1_Re").Clone("T_2D_VBF_4_AC_1_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                           self.dimensions, self.ProjDim,
                           self.varm4l, self.varKD, self.varKD2
                        )
                     )

                  if self.anomCoupl == 1:
                     self.gg_T_4_AC_1_Im_list.append(
                        ExtendedTemplate(
                              self.templateFileList[icat].Get("T_2D_4_AC_1_Im").Clone("T_2D_gg_4_AC_1_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                              self.dimensions, self.ProjDim,
                              self.varm4l, self.varKD, self.varKD2
                           )
                        )
                     self.VBF_T_4_AC_1_Im_list.append(
                        ExtendedTemplate(
                              self.templateFileList[icat].Get("T_2D_VBF_4_AC_1_Im").Clone("T_2D_VBF_4_AC_1_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                              self.dimensions, self.ProjDim,
                              self.varm4l, self.varKD, self.varKD2
                           )
                        )


# Interference ai**2 x a1**(2-2) real and imaginary parts
         # Templates
            # No ggF
            # VBF
            self.VBF_T_4_AC_2_Re_list = []
            self.VBF_T_4_AC_2_Im_list = []

            if self.anomCoupl != 0:
               for icat in range(0,self.nCategories):
                  self.VBF_T_4_AC_2_Re_list.append(
                     ExtendedTemplate(
                           self.templateFileList[icat].Get("T_2D_VBF_4_AC_2_Re").Clone("T_2D_VBF_4_AC_2_Re_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                           self.dimensions, self.ProjDim,
                           self.varm4l, self.varKD, self.varKD2
                        )
                     )

                  if self.anomCoupl == 1:
                     self.VBF_T_4_AC_2_Im_list.append(
                        ExtendedTemplate(
                              self.templateFileList[icat].Get("T_2D_VBF_4_AC_2_Im").Clone("T_2D_VBF_4_AC_2_Im_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)),
                              self.dimensions, self.ProjDim,
                              self.varm4l, self.varKD, self.varKD2
                           )
                        )


# FORMULAE
      self.ggSigFormula_list = []
      self.ggInterfFormula_list = []
      self.VBFSigFormula_list = []
      self.VBFInterfFormula_list = []

      self.ggSigRFV_list = []
      self.ggInterfRFV_list = []
      self.VBFSigRFV_list = []
      self.VBFInterfRFV_list = []

      # In signals and interferences, @0==muF/V; additionally in interferences, @1==kbkg_gg/VBF
      if self.anomCoupl == 1: # Full parameterization with fai1=[-1, 1], and phases phiai1 and phia1_gg, phia1_VBF (phia1 are different since the bkg phase could be different)
         # 0-3 are templates, @1==fai1, @2==phiai1
         self.ggSigFormula_list.append("@0*(1-abs(@1))")
         self.ggSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))*cos(@2)")
         self.ggSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))*sin(@2)")
         self.ggSigFormula_list.append("@0*abs(@1)")

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
            # 0-3 are templates, @2==fai1, @3==phiai1, @4==phia1 (gg)
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sqrt(1-abs(@2))*cos(@4)")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sqrt(1-abs(@2))*sin(@4)")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))*cos(@3+@4)")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))*sin(@3+@4)")

            # 0-5 are templates, @2==fai1, @3==phiai1, @4==phia1 (VBF)
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*(1-abs(@2))*cos(2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*(1-abs(@2))*sin(2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2)*(1-abs(@2)))*cos(@3+2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2)*(1-abs(@2)))*sin(@3+2*@4)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*abs(@2)*cos(2*(@3+@4))")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*abs(@2)*sin(2*(@3+@4))")
      elif self.anomCoupl == 2: # No phases, just fai1=[-1, 1]
         # 0-2 are templates, @1==fai1
         self.ggSigFormula_list.append("@0*(1-abs(@1))")
         self.ggSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))")
         self.ggSigFormula_list.append("@0*abs(@1)")

         # 0-4 are templates, @1==fai1
         self.VBFSigFormula_list.append("@0*pow((1-abs(@1)),2)")
         self.VBFSigFormula_list.append("@0*sign(@1)*sqrt(abs(@1))*pow(sqrt(1-abs(@1)),3)")
         self.VBFSigFormula_list.append("@0*abs(@1)*(1-abs(@1))")
         self.VBFSigFormula_list.append("@0*sign(@1)*pow(sqrt(abs(@1)),3)*sqrt(1-abs(@1))")
         self.VBFSigFormula_list.append("@0*pow(@1,2)")

         if(not(self.isBkgSigOnly)):
            # 0-1 are templates, @2==fai1
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sqrt(1-abs(@2))")
            self.ggInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))")

            # 0-2 are templates, @2==fai1
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*(1-abs(@2))")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*sign(@2)*sqrt(abs(@2)*(1-abs(@2)))")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)*abs(@2)")
      else: # No ai1 dependence
         self.ggSigFormula_list.append("@0")
         self.VBFSigFormula_list.append("@0")
         if(not(self.isBkgSigOnly)):
         self.ggInterfFormula_list.append("sqrt(@0*@1)")
            self.VBFInterfFormula_list.append("sqrt(@0*@1)")

      for irfv in range(0,len(self.ggSigFormula_list)):
         rfvname = "ggSig_AC_{0:.0f}_Coef".format(irfv)
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

      for irfv in range(0,len(self.VBFSigFormula_list)):
         rfvname = "VBFSig_AC_{0:.0f}_Coef".format(irfv)
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

      for irfv in range(0,len(self.ggInterfFormula_list)):
         rfvname = "ggInterf_AC_{0:.0f}_Coef".format(irfv)
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

      for irfv in range(0,len(self.VBFInterfFormula_list)):
         rfvname = "VBFInterf_AC_{0:.0f}_Coef".format(irfv)
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
      self.ggSigFunctions_Args_list = []
      self.VBFSigFunctions_Args_list = []
      self.ggInterfFunctions_Args_list = []
      self.VBFInterfFunctions_Args_list = []

      self.ggSigFunctions_Args_list.append(self.gg_T_1_list)
      if self.anomCoupl == 1:
         self.ggSigFunctions_Args_list.append(self.gg_T_1_AC_1_Re_list)
         self.ggSigFunctions_Args_list.append(self.gg_T_1_AC_1_Im_list)
         self.ggSigFunctions_Args_list.append(self.gg_T_1_AC_2_Re_list)
      elif self.anomCoupl == 2:
         self.ggSigFunctions_Args_list.append(self.gg_T_1_AC_1_Re_list)
         self.ggSigFunctions_Args_list.append(self.gg_T_1_AC_2_Re_list)
      if len(self.ggSigFunctions_Args_list)!=len(self.ggSigRFV_list):
         sys.exit("Number of ggSig templates {0:.0f} is not equal to number of funcficients {1:.0f}!".format(len(self.ggSigFunctions_Args_list),len(self.ggSigRFV_list)))

      self.VBFSigFunctions_Args_list.append(self.VBF_T_1_list)
      if self.anomCoupl == 1:
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_1_Re_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_1_Im_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_2_PosDef_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_2_Re_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_2_Im_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_3_Re_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_3_Im_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_4_list)
      elif self.anomCoupl == 2:
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_1_Re_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_2_PosDef_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_3_Re_list)
         self.VBFSigFunctions_Args_list.append(self.VBF_T_1_AC_4_list)
      if len(self.VBFSigFunctions_Args_list)!=len(self.VBFSigRFV_list):
         sys.exit("Number of VBFSig templates {0:.0f} is not equal to number of funcficients {1:.0f}!".format(len(self.VBFSigFunctions_Args_list),len(self.VBFSigRFV_list)))

      if(not(self.isBkgSigOnly)):
         self.ggInterfFunctions_Args_list.append(self.gg_T_4_Re_list)
         if self.anomCoupl == 1:
            self.ggInterfFunctions_Args_list.append(self.gg_T_4_Im_list)
            self.ggInterfFunctions_Args_list.append(self.gg_T_4_AC_1_Re_list)
            self.ggInterfFunctions_Args_list.append(self.gg_T_4_AC_1_Im_list)
         elif self.anomCoupl == 2:
            self.ggInterfFunctions_Args_list.append(self.gg_T_4_AC_1_Re_list)
         if len(self.ggInterfFunctions_Args_list)!=len(self.ggInterfRFV_list):
            sys.exit("Number of ggInterf templates {0:.0f} is not equal to number of funcficients {1:.0f}!".format(len(self.ggInterfFunctions_Args_list),len(self.ggInterfRFV_list)))

         self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_Re_list)
         if self.anomCoupl == 1:
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_Im_list)
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_AC_1_Re_list)
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_AC_1_Im_list)
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_AC_2_Re_list)
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_AC_2_Im_list)
         elif self.anomCoupl == 2:
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_AC_1_Re_list)
            self.VBFInterfFunctions_Args_list.append(self.VBF_T_4_AC_2_Re_list)
         if len(self.VBFInterfFunctions_Args_list)!=len(self.VBFInterfRFV_list):
            sys.exit("Number of VBFInterf templates {0:.0f} is not equal to number of funcficients {1:.0f}!".format(len(self.VBFInterfFunctions_Args_list),len(self.VBFInterfRFV_list)))


      self.ggHistFunc_Arg_list = []
      self.VBFHistFunc_Arg_list = []
      for icat in range(0,self.nCategories):
         ral = ROOT.RooArgList()
         for ivar in range(0,len(self.ggSigFunctions_Args_list)):
            ral.add(self.ggSigFunctions_Args_list[ivar][icat].theHistFunc)
         for ivar in range(0,len(self.ggInterfFunctions_Args_list)):
            ral.add(self.ggInterfFunctions_Args_list[ivar][icat].theHistFunc)
         ral.add(self.gg_T_2_list[icat].theHistFunc)
         self.ggHistFunc_Arg_list.append(ral)
      for icat in range(0,self.nCategories):
         ral = ROOT.RooArgList()
         for ivar in range(0,len(self.VBFSigFunctions_Args_list)):
            ral.add(self.VBFSigFunctions_Args_list[ivar][icat].theHistFunc)
         for ivar in range(0,len(self.VBFInterfFunctions_Args_list)):
            ral.add(self.VBFInterfFunctions_Args_list[ivar][icat].theHistFunc)
         ral.add(self.VBF_T_2_list[icat].theHistFunc)
         self.VBFHistFunc_Arg_list.append(ral)


# Construct the p.d.f.'s
      self.ggPdf_list = []
      self.VBFPdf_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         for ivar in range(0,len(self.ggSigRFV_list)):
            rfvargs.add(self.ggSigRFV_list[ivar])
         for ivar in range(0,len(self.ggInterfRFV_list)):
            rfvargs.add(self.ggInterfRFV_list[ivar])
         rfvargs.add(self.kbkg_gg)
         PdfName = "ggPdf_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts))
         ggPdf = ROOT.RooRealSumPdf(
            PdfName, PdfName,
            self.ggHistFunc_Arg_list[icat],ggZZ_funcficients
         )
         self.ggPdf_list.append(ggPdf)
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         for ivar in range(0,len(self.VBFSigRFV_list)):
            rfvargs.add(self.VBFSigRFV_list[ivar])
         for ivar in range(0,len(self.VBFInterfRFV_list)):
            rfvargs.add(self.VBFInterfRFV_list[ivar])
         rfvargs.add(self.kbkg_VBF)
         PdfName = "VBFPdf_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts))
         VBFPdf = ROOT.RooRealSumPdf(
            PdfName, PdfName,
            self.VBFHistFunc_Arg_list[icat],VBFZZ_funcficients
         )
         self.VBFPdf_list.append(VBFPdf)


# Lists of rate FormulaVars
# Each signal, bkg and interf is constructed separately to be able to count their contributions in the end
      self.ggSigRates_RooFormulaVar_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         strformula = ""
         for ivar in range(0,len(self.ggSigFunctions_Args_list)):
            rfvargs.add(self.ggSigRFV_list[ivar])
            rfvargs.add(self.ggSigFunctions_Args_list[ivar][icat].theRate)
            if ivar==0:
               strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
            else:
               strformula = "{0:.0f}*{1:,0f} + {2}".format(2*ivar,2*ivar+1,strformula)
         rfvname = "ggSigRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
         rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
         self.ggSigRates_RooFormulaVar_list.append(rfv)

      self.ggInterfRates_RooFormulaVar_list = []
      if len(self.ggInterfFunctions_Args_list)>0:
         for icat in range(0,self.nCategories):
            rfvargs = ROOT.RooArgList()
            strformula = ""
            for ivar in range(0,len(self.ggInterfFunctions_Args_list)):
               rfvargs.add(self.ggInterfRFV_list[ivar])
               rfvargs.add(self.ggInterfFunctions_Args_list[ivar][icat].theRate)
               if ivar==0:
                  strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
               else:
                  strformula = "{0:.0f}*{1:,0f} + {2}".format(2*ivar,2*ivar+1,strformula)
            rfvname = "ggInterfRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
            rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
            self.ggInterfRates_RooFormulaVar_list.append(rfv)

      self.ggBkgRates_RooFormulaVar_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         strformula = "@0*@1"
         rfvargs.add(self.kbkg_gg)
         rfvargs.add(self.gg_T_2_list[icat].theRate)
         rfvname = "ggBkgRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
         rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
         self.ggBkgRates_RooFormulaVar_list.append(rfv)

      self.VBFSigRates_RooFormulaVar_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         strformula = ""
         for ivar in range(0,len(self.VBFSigFunctions_Args_list)):
            rfvargs.add(self.VBFSigRFV_list[ivar])
            rfvargs.add(self.VBFSigFunctions_Args_list[ivar][icat].theRate)
            if ivar==0:
               strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
            else:
               strformula = "{0:.0f}*{1:,0f} + {2}".format(2*ivar,2*ivar+1,strformula)
         rfvname = "VBFSigRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
         rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
         self.VBFSigRates_RooFormulaVar_list.append(rfv)

      self.VBFInterfRates_RooFormulaVar_list = []
      if len(self.VBFInterfFunctions_Args_list)>0:
         for icat in range(0,self.nCategories):
            rfvargs = ROOT.RooArgList()
            strformula = ""
            for ivar in range(0,len(self.VBFInterfFunctions_Args_list)):
               rfvargs.add(self.VBFInterfRFV_list[ivar])
               rfvargs.add(self.VBFInterfFunctions_Args_list[ivar][icat].theRate)
               if ivar==0:
                  strformula = "{0:.0f}*{1:,0f}".format(2*ivar,2*ivar+1)
               else:
                  strformula = "{0:.0f}*{1:,0f} + {2}".format(2*ivar,2*ivar+1,strformula)
            rfvname = "VBFInterfRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
            self.VBFInterfRates_RooFormulaVar_list.append(rfv)

      self.VBFBkgRates_RooFormulaVar_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         strformula = "@0*@1"
         rfvargs.add(self.kbkg_VBF)
         rfvargs.add(self.VBF_T_2_list[icat].theRate)
         rfvname = "VBFBkgRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
         rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
         self.VBFBkgRates_RooFormulaVar_list.append(rfv)


# Construct total rates
      self.ggTotalRate_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         strformula = "@0+@1"
         if(not(self.isBkgSigOnly)):
            strformula = "@0+@1+@2"
         rfvargs.add(self.ggSigRates_RooFormulaVar_list[icat])
         rfvargs.add(self.ggInterfRates_RooFormulaVar_list[icat])
         rfvargs.add(self.ggBkgRates_RooFormulaVar_list[icat])
         rfvname = "ggTotalRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
         rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
         self.ggTotalRate_list.append(rfv)

      self.VBFTotalRate_list = []
      for icat in range(0,self.nCategories):
         rfvargs = ROOT.RooArgList()
         strformula = "@0+@1"
         if(not(self.isBkgSigOnly)):
            strformula = "@0+@1+@2"
         rfvargs.add(self.VBFSigRates_RooFormulaVar_list[icat])
         rfvargs.add(self.VBFInterfRates_RooFormulaVar_list[icat])
         rfvargs.add(self.VBFBkgRates_RooFormulaVar_list[icat])
         rfvname = "VBFTotalRate_{0}_{1}_{2}_{3}TeV".format(systName,self.catNameList[icat],self.channel,self.sqrts)
         rfv = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
         self.VBFTotalRate_list.append(rfv)



