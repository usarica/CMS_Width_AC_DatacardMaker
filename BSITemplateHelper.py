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

class BSITemplateHelper:
   def __init__(self, options, theMaker, theEqnsMaker, theCategorizer, procname, proctype, templateFileName, iCat, systName):
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.channel = theMaker.channel
      self.theChannelName = theMaker.theChannelName
      self.workspace = theMaker.workspace

      # RooRealVars from the equations maker class
      self.muF = theEqnsMaker.rrvars["muF"] # Could itself be a RooFormulaVar (e.g. Off-shell: muF = R*RF*x. On-shell: muF = R*RF)
      self.muV = theEqnsMaker.rrvars["muV"] # Could itself be a RooFormulaVar (e.g. Off-shell: muV = R*RV*x. On-shell: muV = R*RV)
      self.kbkg_gg = theEqnsMaker.rrvars["kbkg_gg"]
      self.kbkg_VBF = theEqnsMaker.rrvars["kbkg_VBF"]
      self.fai1 = theEqnsMaker.rrvars["fai1"]
      self.phiai1 = theEqnsMaker.rrvars["phiai1"]
      self.phia1_gg = theEqnsMaker.rrvars["phia1_gg"] # Could itself be a RooFormulaVar (e.g. phia1_gg = phia1+phi_SB_gg)
      self.phia1_VBF = theEqnsMaker.rrvars["phia1_VBF"] # Could itself be a RooFormulaVar (e.g. phia1_VBF = phia1+phi_SB_VBF/2)

      self.mass = theMaker.mass
      self.KD1 = theMaker.KD1
      self.KD2 = theMaker.KD2
      self.KD3 = theMaker.KD3

      self.mLow = options.mLow
      self.mHigh = options.mHigh
      self.anomCoupl = options.anomCouplIndex
      self.templateDir = options.templateDir
      self.dimensions = options.dimensions # Number of template dimensions>0


      self.iCatScheme = theCategorizer.iCatScheme
      self.catNameList = theCategorizer.catNameList
      self.nCategories = theCategorizer.nCategories
      self.iCat = iCat
      if self.iCat>=self.nCategories:
         sys.exit("self.iCat={} >= self.nCategories={}!".format(self.iCat,self.nCategories))

      self.templateFileName = templateFileName
      self.systName = systName
      self.templateSuffix = "{0}_{1}_{2}_{3:.0f}TeV".format(self.systName,self.catNameList[self.iCat],self.theChannelName,self.sqrts)

      # To be reset later
      self.nbinsx=(self.mHigh - self.mLow) / 20
      self.nbinsy=30
      self.nbinsz=30
      self.templateXLow=self.mLow
      self.templateYLow=0
      self.templateZLow=0
      self.templateXHigh=self.mHigh
      self.templateYHigh=1
      self.templateZhigh=1

      # Template file
      self.templateFile = None

      self.procname = procname
      self.proctype = proctype
      self.isGGVVLikeCouplings = self.procname.lower().startswith("gg") or self.procname.lower().startswith("tt")
      self.isSigOnly = self.proctype==1

   # Extended template lists
   # Bare SM
      # ggF
      self.gg_T_Sig = None
      self.gg_T_Bkg = None
      self.gg_T_Int_Re = None
      self.gg_T_Int_Im = None
      # VBF
      self.VBF_T_Sig = None
      self.VBF_T_Bkg = None
      self.VBF_T_Int_Re = None
      self.VBF_T_Int_Im = None
   # Signal ai**1 x a1**(2/4-1) real and imaginary parts
      # ggF
      self.gg_T_Sig_ai1_1_Re = None
      self.gg_T_Sig_ai1_1_Im = None
      # VBF
      self.VBF_T_Sig_ai1_1_Re = None
      self.VBF_T_Sig_ai1_1_Im = None
   # Signal ai**2 x a1**(2/4-2) real and imaginary parts
      # ggF
      self.gg_T_Sig_ai1_2_Re = None
      # VBF
      self.VBF_T_Sig_ai1_2_Re = None
      self.VBF_T_Sig_ai1_2_Im = None
      self.VBF_T_Sig_ai1_2_PosDef = None
   # Signal ai**3 x a1**1 real and imaginary parts
      # No ggF
      # VBF
      self.VBF_T_Sig_ai1_3_Re = None
      self.VBF_T_Sig_ai1_3_Im = None
   # Signal ai**4 x a1**0 real and imaginary parts
      # No ggF
      # VBF
      self.VBF_T_Sig_ai1_4 = None

   # Interference ai**1 x a1**(1/2-1) real and imaginary parts
      # ggF
      self.gg_T_Int_ai1_1_Re = None
      self.gg_T_Int_ai1_1_Im = None
      # VBF
      self.VBF_T_Int_ai1_1_Re = None
      self.VBF_T_Int_ai1_1_Im = None
   # Interference ai**2 x a1**(2-2) real and imaginary parts
      # No ggF
      # VBF
      self.VBF_T_Int_ai1_2_Re = None
      self.VBF_T_Int_ai1_2_Im = None

   # BSI+AC FORMULAE
      # Bare formula strings
      self.ggSigFormula_list = theEqnsMaker.ggSigFormula_list
      self.ggInterfFormula_list = theEqnsMaker.ggInterfFormula_list
      self.VBFSigFormula_list = theEqnsMaker.VBFSigFormula_list
      self.VBFInterfFormula_list = theEqnsMaker.VBFInterfFormula_list
      # RooFormulaVars
      self.ggSigRFV_list = theEqnsMaker.ggSigRFV_list
      self.ggInterfRFV_list = theEqnsMaker.ggInterfRFV_list
      self.VBFSigRFV_list = theEqnsMaker.VBFSigRFV_list
      self.VBFInterfRFV_list = theEqnsMaker.VBFInterfRFV_list

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


# Import the pdf
   def importToWorkspace(self):
      if self.ggPdf is not None:
         getattr(self.workspace, 'import')(self.ggPdf, ROOT.RooFit.RecycleConflictNodes())
      if self.ggTotalRate is not None:
         getattr(self.workspace, 'import')(self.ggTotalRate, ROOT.RooFit.RecycleConflictNodes())
      if self.VBFPdf is not None:
         getattr(self.workspace, 'import')(self.VBFPdf, ROOT.RooFit.RecycleConflictNodes())
      if self.VBFTotalRate is not None:
         getattr(self.workspace, 'import')(self.VBFTotalRate, ROOT.RooFit.RecycleConflictNodes())


# Open the template files
   def openFile(self):
      print "Opening file ",self.templateFileName
      self.templateFile = ROOT.TFile.Open(self.templateFileName, "read")
      if self.templateFile is None:
         raise RuntimeError("BSITemplateHelper file {} is None!".format(self.templateFileName))
      elif self.templateFile.IsZombie():
         raise RuntimeError("BSITemplateHelper could not open file {}!".format(self.templateFileName))
# Close the template files
   def close(self):
      if self.templateFile is not None:
         if self.templateFile.IsOpen():
            self.templateFile.Close()


   def getThePdf(self):
      if self.ggPdf is not None:
         return self.ggPdf
      else:
         return self.VBFPdf
   def getTheRate(self):
      if self.ggTotalRate is not None:
         return self.ggTotalRate
      else:
         return self.VBFTotalRate


# Get shapes for each category
   def getTemplates(self,processName=None,templatePrefix="T"):
      self.openFile()

      self.templatePrefix = templatePrefix
      if processName is not None:
         self.processName = processName
         self.templatePrefix = "{}_{}".format(self.templatePrefix,processName)
      elif not self.isGGVVLikeCouplings:
         self.processName = "VBF"
         self.templatePrefix = "{}_{}".format(self.templatePrefix,processName)
      else:
         self.processName = "gg"
         self.templatePrefix = "{}_{}".format(self.templatePrefix,processName)

      if self.isGGVVLikeCouplings:
         self.getTemplates_ggVVLike()
      else:
         self.getTemplates_vvVVLike()


   def getTemplates_ggVVLike(self):
#---------- SM SIGNAL AND BACKGROUND TEMPLATES -------------
# Bare SM
      self.gg_T_Sig = ExtendedTemplate(
               self.templateFile.Get("{}_Sig".format(self.templatePrefix)).Clone("{}_Sig_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions,
               self.KD1, self.KD2, self.KD3
            )
      if(not(self.isSigOnly)):
         self.gg_T_Bkg = ExtendedTemplate(
                  self.templateFile.Get("{}_Bkg".format(self.templatePrefix)).Clone("{}_Bkg_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         self.gg_T_Int_Re = ExtendedTemplate(
                  self.templateFile.Get("{}_Int_Re".format(self.templatePrefix)).Clone("{}_Int_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         if (self.anomCoupl==1):
            self.gg_T_Int_Im = ExtendedTemplate(
                     self.templateFile.Get("{}_Int_Im".format(self.templatePrefix)).Clone("{}_Int_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )

      if self.anomCoupl != 0:
#-----------------------------------------------------------------------#
#                        SIGNAL AC TERMS
#-----------------------------------------------------------------------#
# Signal ai**1 x a1**(2/4-1) real and imaginary parts
         self.gg_T_Sig_ai1_1_Re = ExtendedTemplate(
                  self.templateFile.Get("{}_Sig_ai1_1_Re".format(self.templatePrefix)).Clone("{}_Sig_ai1_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         if self.anomCoupl == 1:
            self.gg_T_Sig_ai1_1_Im = ExtendedTemplate(
                     self.templateFile.Get("{}_Sig_ai1_1_Im".format(self.templatePrefix)).Clone("{}_Sig_ai1_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )

# Signal ai**2 x a1**(2/4-2) real and imaginary parts
         self.gg_T_Sig_ai1_2_Re = ExtendedTemplate(
                  self.templateFile.Get("{}_Sig_ai1_2".format(self.templatePrefix)).Clone("{}_Sig_ai1_2_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         #if self.anomCoupl == 1:
         #   continue # No ggF term
         #elif self.anomCoupl == 2: # if self.anomCoupl == 2, PosDef = PosDef+Re
         #   continue # No ggF term

# Signal ai**3 x a1**1 real and imaginary parts
   # No ggF term

# Signal ai**4 x a1**0 real and imaginary parts
   # No ggF term

#-----------------------------------------------------------------------#
#                       INTERFERENCE AC TERMS
#-----------------------------------------------------------------------#
         if(not(self.isSigOnly)):
# Interference ai**1 x a1**(1/2-1) real and imaginary parts
            self.gg_T_Int_ai1_1_Re = ExtendedTemplate(
                     self.templateFile.Get("{}_Int_ai1_1_Re".format(self.templatePrefix)).Clone("{}_Int_ai1_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )
            if self.anomCoupl == 1:
               self.gg_T_Int_ai1_1_Im = ExtendedTemplate(
                        self.templateFile.Get("{}_Int_ai1_1_Im".format(self.templatePrefix)).Clone("{}_Int_ai1_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                        self.dimensions,
                        self.KD1, self.KD2, self.KD3
                     )

# Interference ai**2 x a1**(2-2) real and imaginary parts
   # No ggF term

# Lists of template arguments
      self.ggSigFunctions_Args.append(self.gg_T_Sig)
      if self.anomCoupl == 1:
         self.ggSigFunctions_Args.append(self.gg_T_Sig_ai1_1_Re)
         self.ggSigFunctions_Args.append(self.gg_T_Sig_ai1_1_Im)
         self.ggSigFunctions_Args.append(self.gg_T_Sig_ai1_2_Re)
      elif self.anomCoupl == 2:
         self.ggSigFunctions_Args.append(self.gg_T_Sig_ai1_1_Re)
         self.ggSigFunctions_Args.append(self.gg_T_Sig_ai1_2_Re)
      if len(self.ggSigFunctions_Args)!=len(self.ggSigRFV_list):
         sys.exit("Number of {0}Sig templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.ggSigFunctions_Args),len(self.ggSigRFV_list)))

      if(not(self.isSigOnly)):
         self.ggInterfFunctions_Args.append(self.gg_T_Int_Re)
         if self.anomCoupl == 1:
            self.ggInterfFunctions_Args.append(self.gg_T_Int_Im)
            self.ggInterfFunctions_Args.append(self.gg_T_Int_ai1_1_Re)
            self.ggInterfFunctions_Args.append(self.gg_T_Int_ai1_1_Im)
         elif self.anomCoupl == 2:
            self.ggInterfFunctions_Args.append(self.gg_T_Int_ai1_1_Re)
         if len(self.ggInterfFunctions_Args)!=len(self.ggInterfRFV_list):
            sys.exit("Number of {0}Interf templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.ggInterfFunctions_Args),len(self.ggInterfRFV_list)))

      self.ggHistFunc_Arg = ROOT.RooArgList()
      for var in self.ggSigFunctions_Args:
         self.ggHistFunc_Arg.add(var.theHistFunc)
      for var in self.ggInterfFunctions_Args:
         self.ggHistFunc_Arg.add(var.theHistFunc)
      if(not(self.isSigOnly)):
         self.ggHistFunc_Arg.add(self.gg_T_Bkg.theHistFunc)

# Construct the p.d.f.'s
      rfvargs = ROOT.RooArgList()
      for var in self.ggSigRFV_list:
         rfvargs.add(var)
      if(not(self.isSigOnly)):
         for var in self.ggInterfRFV_list:
            rfvargs.add(var)
         rfvargs.add(self.kbkg_gg)
      PdfName = "{}Pdf_{}".format(self.processName,self.templateSuffix)
      self.ggPdf = ROOT.RooRealFlooredSumPdf(
         PdfName, PdfName,
         self.ggHistFunc_Arg, rfvargs
      )
      self.ggHistFunc_Arg.Print("v")
      rfvargs.Print("v")
      print PdfName,"value =",self.ggPdf.getVal()

# Lists of rate FormulaVars
# Each signal, bkg and interf is constructed separately to be able to count their contributions in the end
      rfvargs = ROOT.RooArgList()
      strformula = ""
      for ivar in range(0,len(self.ggSigFunctions_Args)):
         rfvargs.add(self.ggSigRFV_list[ivar])
         rfvargs.add(self.ggSigFunctions_Args[ivar].theRate)
         if ivar==0:
            strformula = "@{0:.0f}*@{1:.0f}".format(2*ivar,2*ivar+1)
         else:
            strformula = "{2} + @{0:.0f}*@{1:.0f}".format(2*ivar,2*ivar+1,strformula)
      rfvname = "{}SigRate_{}".format(self.processName,self.templateSuffix)
      self.ggSigRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      if len(self.ggInterfFunctions_Args)>0:
         rfvargs = ROOT.RooArgList()
         strformula = ""
         for ivar in range(0,len(self.ggInterfFunctions_Args)):
            rfvargs.add(self.ggInterfRFV_list[ivar])
            rfvargs.add(self.ggInterfFunctions_Args[ivar].theRate)
            if ivar==0:
               strformula = "@{0:.0f}*@{1:,0f}".format(2*ivar,2*ivar+1)
            else:
               strformula = "{2} + @{0:.0f}*@{1:,0f}".format(2*ivar,2*ivar+1,strformula)
         rfvname = "{}InterfRate_{}".format(self.processName,self.templateSuffix)
         self.ggInterfRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      self.ggBkgRates_RooFormulaVar = None
      if(not(self.isSigOnly)):
         rfvargs = ROOT.RooArgList()
         strformula = "@0*@1"
         rfvargs.add(self.kbkg_gg)
         rfvargs.add(self.gg_T_Bkg.theRate)
         rfvname = "{}BkgRate_{}".format(self.processName,self.templateSuffix)
         self.ggBkgRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

# Construct total rates
      rfvargs = ROOT.RooArgList()
      strformula = "@0"
      rfvargs.add(self.ggSigRates_RooFormulaVar)
      if(not(self.isSigOnly)):
         strformula = "@0+@1+@2"
         rfvargs.add(self.ggInterfRates_RooFormulaVar)
         rfvargs.add(self.ggBkgRates_RooFormulaVar)
      rfvname = "{}TotalRate_{}".format(self.processName,self.templateSuffix)
      self.ggTotalRate = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
      print rfvname,"total rate =",self.ggTotalRate.getVal()


   def getTemplates_vvVVLike(self):
#---------- SM SIGNAL AND BACKGROUND TEMPLATES -------------
# Bare SM
      self.VBF_T_Sig = ExtendedTemplate(
               self.templateFile.Get("{}_Sig".format(self.templatePrefix)).Clone("{}_Sig_{}".format(self.templatePrefix,self.templateSuffix)),
               self.dimensions,
               self.KD1, self.KD2, self.KD3
            )
      if(not(self.isSigOnly)):
         self.VBF_T_Bkg = ExtendedTemplate(
                  self.templateFile.Get("{}_Bkg".format(self.templatePrefix)).Clone("{}_Bkg_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         self.VBF_T_Int_Re = ExtendedTemplate(
                  self.templateFile.Get("{}_Int_Re".format(self.templatePrefix)).Clone("{}_Int_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         if (self.anomCoupl==1):
            self.VBF_T_Int_Im = ExtendedTemplate(
                     self.templateFile.Get("{}_Int_Im".format(self.templatePrefix)).Clone("{}_Int_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )

      if self.anomCoupl != 0:
#-----------------------------------------------------------------------#
#                        SIGNAL AC TERMS
#-----------------------------------------------------------------------#
# Signal ai**1 x a1**(2/4-1) real and imaginary parts
         self.VBF_T_Sig_ai1_1_Re = ExtendedTemplate(
                  self.templateFile.Get("{}_Sig_ai1_1_Re".format(self.templatePrefix)).Clone("{}_Sig_ai1_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )
         if self.anomCoupl == 1:
            self.VBF_T_Sig_ai1_1_Im = ExtendedTemplate(
                     self.templateFile.Get("{}_Sig_ai1_1_Im".format(self.templatePrefix)).Clone("{}_Sig_ai1_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )

# Signal ai**2 x a1**(2/4-2) real and imaginary parts
         if self.anomCoupl == 1:
            self.VBF_T_Sig_ai1_2_PosDef = ExtendedTemplate(
                     self.templateFile.Get("{}_Sig_ai1_2_PosDef".format(self.templatePrefix)).Clone("{}_Sig_ai1_2_PosDef_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )
            self.VBF_T_Sig_ai1_2_Re = ExtendedTemplate(
                     self.templateFile.Get("{}_Sig_ai1_2_Re".format(self.templatePrefix)).Clone("{}_Sig_ai1_2_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )
            self.VBF_T_Sig_ai1_2_Im = ExtendedTemplate(
                     self.templateFile.Get("{}_Sig_ai1_2_Im".format(self.templatePrefix)).Clone("{}_Sig_ai1_2_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )
         elif self.anomCoupl == 2: # if self.anomCoupl == 2, PosDef = PosDef+Re
            tmpTpl = self.templateFile.Get("{}_Sig_ai1_2_PosDef".format(self.templatePrefix))
            tmpTpl2 = self.templateFile.Get("{}_Sig_ai1_2_Re".format(self.templatePrefix))
            tmpTpl3 = None
            if tmpTpl:
               tmpTpl3 = tmpTpl.Clone("{}_Sig_ai1_2_PosDef_{}".format(self.templatePrefix,self.templateSuffix))
               if tmpTpl2:
                  tmpTpl3.Add(tmpTpl2)
            elif tmpTpl2:
               tmpTpl3 = tmpTpl2.Clone("{}_Sig_ai1_2_PosDef_{}".format(self.templatePrefix,self.templateSuffix))
            self.VBF_T_Sig_ai1_2_PosDef = ExtendedTemplate(
                     tmpTpl3,
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )

# Signal ai**3 x a1**1 real and imaginary parts
         self.VBF_T_Sig_ai1_3_Re = ExtendedTemplate(
                  self.templateFile.Get("{}_Sig_ai1_3_Re".format(self.templatePrefix)).Clone("{}_Sig_ai1_3_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )

         if self.anomCoupl == 1:
            self.VBF_T_Sig_ai1_3_Im = ExtendedTemplate(
                     self.templateFile.Get("{}_Sig_ai1_3_Im".format(self.templatePrefix)).Clone("{}_Sig_ai1_3_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )

# Signal ai**4 x a1**0 real and imaginary parts
         self.VBF_T_Sig_ai1_4 = ExtendedTemplate(
                  self.templateFile.Get("{}_Sig_ai1_4".format(self.templatePrefix)).Clone("{}_Sig_ai1_4_{}".format(self.templatePrefix,self.templateSuffix)),
                  self.dimensions,
                  self.KD1, self.KD2, self.KD3
               )

#-----------------------------------------------------------------------#
#                       INTERFERENCE AC TERMS
#-----------------------------------------------------------------------#
         if(not(self.isSigOnly)):
# Interference ai**1 x a1**(1/2-1) real and imaginary parts
            self.VBF_T_Int_ai1_1_Re = ExtendedTemplate(
                     self.templateFile.Get("{}_Int_ai1_1_Re".format(self.templatePrefix)).Clone("{}_Int_ai1_1_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )
            if self.anomCoupl == 1:
               self.VBF_T_Int_ai1_1_Im = ExtendedTemplate(
                        self.templateFile.Get("{}_Int_ai1_1_Im".format(self.templatePrefix)).Clone("{}_Int_ai1_1_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                        self.dimensions,
                        self.KD1, self.KD2, self.KD3
                     )

# Interference ai**2 x a1**(2-2) real and imaginary parts
            self.VBF_T_Int_ai1_2_Re = ExtendedTemplate(
                     self.templateFile.Get("{}_Int_ai1_2_Re".format(self.templatePrefix)).Clone("{}_Int_ai1_2_Re_{}".format(self.templatePrefix,self.templateSuffix)),
                     self.dimensions,
                     self.KD1, self.KD2, self.KD3
                  )
            if self.anomCoupl == 1:
               self.VBF_T_Int_ai1_2_Im = ExtendedTemplate(
                        self.templateFile.Get("{}_Int_ai1_2_Im".format(self.templatePrefix)).Clone("{}_Int_ai1_2_Im_{}".format(self.templatePrefix,self.templateSuffix)),
                        self.dimensions,
                        self.KD1, self.KD2, self.KD3
                     )

# Lists of template arguments
      self.VBFSigFunctions_Args.append(self.VBF_T_Sig)
      if self.anomCoupl == 1:
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_1_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_1_Im)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_2_PosDef)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_2_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_2_Im)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_3_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_3_Im)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_4)
      elif self.anomCoupl == 2:
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_1_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_2_PosDef)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_3_Re)
         self.VBFSigFunctions_Args.append(self.VBF_T_Sig_ai1_4)
      if len(self.VBFSigFunctions_Args)!=len(self.VBFSigRFV_list):
         sys.exit("Number of {0}Sig templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.VBFSigFunctions_Args),len(self.VBFSigRFV_list)))

      if(not(self.isSigOnly)):
         self.VBFInterfFunctions_Args.append(self.VBF_T_Int_Re)
         if self.anomCoupl == 1:
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_Im)
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_ai1_1_Re)
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_ai1_1_Im)
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_ai1_2_Re)
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_ai1_2_Im)
         elif self.anomCoupl == 2:
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_ai1_1_Re)
            self.VBFInterfFunctions_Args.append(self.VBF_T_Int_ai1_2_Re)
         if len(self.VBFInterfFunctions_Args)!=len(self.VBFInterfRFV_list):
            sys.exit("Number of {0}Interf templates {1:.0f} is not equal to number of funcficients {2:.0f}!".format(self.processName,len(self.VBFInterfFunctions_Args),len(self.VBFInterfRFV_list)))

      self.VBFHistFunc_Arg = ROOT.RooArgList()
      for var in self.VBFSigFunctions_Args:
         self.VBFHistFunc_Arg.add(var.theHistFunc)
      for var in self.VBFInterfFunctions_Args:
         self.VBFHistFunc_Arg.add(var.theHistFunc)
      if(not(self.isSigOnly)):
         self.VBFHistFunc_Arg.add(self.VBF_T_Bkg.theHistFunc)

# Construct the p.d.f.'s
      rfvargs = ROOT.RooArgList()
      for var in self.VBFSigRFV_list:
         rfvargs.add(var)
      if(not(self.isSigOnly)):
         for var in self.VBFInterfRFV_list:
            rfvargs.add(var)
         rfvargs.add(self.kbkg_VBF)
      PdfName = "{}Pdf_{}".format(self.processName,self.templateSuffix)
      self.VBFPdf = ROOT.RooRealFlooredSumPdf(
         PdfName, PdfName,
         self.VBFHistFunc_Arg, rfvargs
      )
      self.VBFHistFunc_Arg.Print("v")
      rfvargs.Print("v")
      print PdfName,"value =",self.VBFPdf.getVal()

# Lists of rate FormulaVars
# Each signal, bkg and interf is constructed separately to be able to count their contributions in the end
      rfvargs = ROOT.RooArgList()
      strformula = ""
      for ivar in range(0,len(self.VBFSigFunctions_Args)):
         rfvargs.add(self.VBFSigRFV_list[ivar])
         rfvargs.add(self.VBFSigFunctions_Args[ivar].theRate)
         if ivar==0:
            strformula = "@{0:.0f}*@{1:.0f}".format(2*ivar,2*ivar+1)
         else:
            strformula = "{2} + @{0:.0f}*@{1:.0f}".format(2*ivar,2*ivar+1,strformula)
      rfvname = "{}SigRate_{}".format(self.processName,self.templateSuffix)
      self.VBFSigRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      if len(self.VBFInterfFunctions_Args)>0:
         rfvargs = ROOT.RooArgList()
         strformula = ""
         for ivar in range(0,len(self.VBFInterfFunctions_Args)):
            rfvargs.add(self.VBFInterfRFV_list[ivar])
            rfvargs.add(self.VBFInterfFunctions_Args[ivar].theRate)
            if ivar==0:
               strformula = "@{0:.0f}*@{1:.0f}".format(2*ivar,2*ivar+1)
            else:
               strformula = "{2} + @{0:.0f}*@{1:.0f}".format(2*ivar,2*ivar+1,strformula)
         rfvname = "{}InterfRate_{}".format(self.processName,self.templateSuffix)
         self.VBFInterfRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

      self.VBFBkgRates_RooFormulaVar = None
      if(not(self.isSigOnly)):
         rfvargs = ROOT.RooArgList()
         strformula = "@0*@1"
         rfvargs.add(self.kbkg_VBF)
         rfvargs.add(self.VBF_T_Bkg.theRate)
         rfvname = "{}BkgRate_{}".format(self.processName,self.templateSuffix)
         self.VBFBkgRates_RooFormulaVar = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )

# Construct total rates
      rfvargs = ROOT.RooArgList()
      strformula = "@0"
      rfvargs.add(self.VBFSigRates_RooFormulaVar)
      if(not(self.isSigOnly)):
         strformula = "@0+@1+@2"
         rfvargs.add(self.VBFInterfRates_RooFormulaVar)
         rfvargs.add(self.VBFBkgRates_RooFormulaVar)
      rfvname = "{}TotalRate_{}".format(self.processName,self.templateSuffix)
      self.VBFTotalRate = ROOT.RooFormulaVar( rfvname , strformula , rfvargs )
      print rfvname,"total rate =",self.VBFTotalRate.getVal()



