#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array


# Container class for the equations pertaining to different datacard models
# Note: This class is independent of the input card, so luminosity and other variables common per-sqrts have to be created somewhere else.
class EquationsMaker:
   def __init__(self, options):
      self.low_M = options.mLow
      self.high_M = options.mHigh
      self.anomCoupl = options.anomCouplIndex
      self.isBkgSigOnly = options.isBkgSigOnly
      self.GHmodel = options.GHmodel
      self.GHrefval = options.GHrefval

   # Variables for the template dimensions
      varname = "CMS_zz4l_widthMass"
      self.varm4l = ROOT.RooRealVar(varname, varname, self.low_M, self.high_M)
      self.varm4l.setBins(69) # To be reset later
      varname = "CMS_zz4l_widthKD"
      self.varKD = ROOT.RooRealVar(varname, varname, 0., 1.)
      self.varKD.setBins(30) # To be reset later
      varname = "CMS_zz4l_widthKD2"
      self.varKD2 = ROOT.RooRealVar(varname, varname, 0., 1.)
      self.varKD2.setBins(30) # To be reset later
      varname = "CMS_zz4l_widthKDint"
      self.varKDint = ROOT.RooRealVar(varname, varname, -1., 1.)
      self.varKDint.setBins(30) # To be reset later

   # Variables for signal and bkg strength
      self.muF = None
      self.muV = None

      varname = "R"
      self.R = ROOT.RooRealVar(varname, varname, 1, 0, 100)
      self.R.setVal(1)
      self.R.setBins(100)
      varname = "RF"
      self.RF = ROOT.RooRealVar(varname, varname, 1, 0, 100)
      self.RF.setVal(1)
      self.RF.setBins(100)
      varname = "RV"
      self.RV = ROOT.RooRealVar(varname, varname, 1, 0, 100)
      self.RV.setVal(1)
      self.RV.setBins(100)
      varname = "GHratio"
      self.GHratio = ROOT.RooRealVar(varname, varname, 1., 0., 50.)
      self.GHratio.setBins(500)
      varname = "GHrefval"
      self.GH = ROOT.RooRealVar(varname, varname, self.GHrefval)
      self.GH.setConstant(True)
      varname = "kbkg_gg"
      self.kbkg_gg = ROOT.RooRealVar(varname, varname, 1., 0., 2.)
      self.kbkg_gg.setBins(200)
      varname = "kbkg_VBF"
      self.kbkg_VBF = ROOT.RooRealVar(varname, varname, 1., 0., 2.)
      self.kbkg_VBF.setBins(200)

      varname = "fai1"
      self.fai1 = ROOT.RooRealVar(varname, varname, 0., -1., 1.)
      self.fai1.setBins(200)
      varname = "phiai1"
      self.phiai1 = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      self.phiai1.setBins(200)
      varname = "fai2"
      self.fai2 = ROOT.RooRealVar(varname, varname, 0., -1., 1.)
      self.fai2.setBins(200)
      varname = "phiai2"
      self.phiai2 = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      self.phiai2.setBins(200)

      self.phia1_gg = None # Could itself be a RooFormulaVar (e.g. phia1_gg = phia1+phi_SB_gg)
      self.phia1_VBF = None # Could itself be a RooFormulaVar (e.g. phia1_VBF = phia1+phi_SB_VBF/2)
      varname = "phia1"
      self.phia1 = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      self.phia1.setBins(200)
      varname = "phia1_SB_gg"
      self.phia1_SB_gg = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      self.phia1_SB_gg.setBins(200)
      varname = "phia1_SB_VBF"
      self.phia1_SB_VBF = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      self.phia1_SB_VBF.setBins(200)

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


   def makeRFVs_BSI(self):
   # Construct muF/V
      if self.GHmodel==1:
         self.muF = ROOT.RooFormulaVar("muF", "@0*@1*@2", ROOT.RooArgList(self.R,self.RF,self.GHratio))
         self.muV = ROOT.RooFormulaVar("muV", "@0*@1*@2", ROOT.RooArgList(self.R,self.RV,self.GHratio))
      elif self.GHmodel==-1:
         self.muF = ROOT.RooFormulaVar("muF", "@0*@1/@2", ROOT.RooArgList(self.R,self.RF,self.GHratio))
         self.muV = ROOT.RooFormulaVar("muV", "@0*@1/@2", ROOT.RooArgList(self.R,self.RV,self.GHratio))
      elif self.GHmodel==2:
         self.muF = ROOT.RooFormulaVar("muF", "@0*@1*@2*@3", ROOT.RooArgList(self.R,self.RF,self.GHratio,self.GH))
         self.muV = ROOT.RooFormulaVar("muV", "@0*@1*@2*@3", ROOT.RooArgList(self.R,self.RV,self.GHratio,self.GH))
      elif self.GHmodel==-2:
         self.muF = ROOT.RooFormulaVar("muF", "@0*@1/(@2*@3)", ROOT.RooArgList(self.R,self.RF,self.GHratio,self.GH))
         self.muV = ROOT.RooFormulaVar("muV", "@0*@1/(@2*@3)", ROOT.RooArgList(self.R,self.RV,self.GHratio,self.GH))
      else:
         self.muF = ROOT.RooFormulaVar("muF", "@0*@1", ROOT.RooArgList(self.R,self.RF))
         self.muV = ROOT.RooFormulaVar("muV", "@0*@1", ROOT.RooArgList(self.R,self.RV))

      self.phia1_gg = ROOT.RooFormulaVar("phia1_gg", "@0+@1", ROOT.RooArgList(self.phia1,self.phia1_SB_gg))
      self.phia1_VBF = ROOT.RooFormulaVar("phia1_VBF", "@0+@1/2.", ROOT.RooArgList(self.phia1,self.phia1_SB_VBF))

   # ggVV FORMULAE
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

   # vvVV FORMULAE
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


