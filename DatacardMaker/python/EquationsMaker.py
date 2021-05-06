#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array


# Container class for the equations pertaining to different datacard models
# Note: This class is independent of the input card, so luminosity and other variables common per-sqrts have to be created somewhere else.
class EquationsMaker:
   def __init__(self, options, theInputCard):
      self.mPOLE = options.mPOLE
      self.mLow = options.mLow
      self.mHigh = options.mHigh
      self.anomCoupl = options.anomCouplIndex
      self.GHmodel = options.GHmodel
      self.GHrefval = options.GHrefval
      self.sqrts = theInputCard.sqrts
      self.theSqrtsPeriod = theInputCard.theSqrtsPeriod
      self.lumi = theInputCard.lumi

      # Set the GHmodel name: Important to distinguish on-shell vs off-shell
      if self.GHmodel==1:
         self.GHmodelName="OnshellR_times_GGsm"
      elif self.GHmodel==-1:
         self.GHmodelName="OnshellR_over_GGsm"
      elif self.GHmodel==2:
         self.GHmodelName="OnshellR_times_GGsm_times_GHref"
      elif self.GHmodel==-2:
         self.GHmodelName="OnshellR_over_GGsm_over_GHref"
      else:
         self.GHmodelName="OnshellR_only"

      self.rrvars = dict()

      # MH
      var = ROOT.RooRealVar("MH", "MH", self.mPOLE, 0, self.sqrts*1000.)
      var.setConstant(True)
      self.rrvars["MH"]=var

      # LUMI
      var = ROOT.RooConstVar("LUMI_{0}".format(self.theSqrtsPeriod), "LUMI_{0}".format(self.theSqrtsPeriod), self.lumi)
      self.rrvars["lumi"]=var

   # Variables for the template dimensions
      varname = "mass"
      var = ROOT.RooRealVar(varname, varname, self.mLow, self.mHigh)
      var.setBins(69) # To be reset later
      self.rrvars["mass"]=var
      varname = "KD1"
      var = ROOT.RooRealVar(varname, varname, 0., 0., 1.)
      var.setBins(30) # To be reset later
      self.rrvars["KD1"]=var
      varname = "KD2"
      var = ROOT.RooRealVar(varname, varname, 0., 0., 1.)
      var.setBins(30) # To be reset later
      self.rrvars["KD2"]=var
      varname = "KD3"
      var = ROOT.RooRealVar(varname, varname, 0., 0., 1.)
      var.setBins(30) # To be reset later
      self.rrvars["KD3"]=var
      varname = "weight"
      var = ROOT.RooRealVar(varname, varname, 1., 0., 10.)
      var.removeMin()
      var.removeMax()
      self.rrvars["weight"]=var

      Rnames = [ "R", "RV", "RF", "R_{0:.0f}TeV".format(self.sqrts), "RV_{0:.0f}TeV".format(self.sqrts), "RF_{0:.0f}TeV".format(self.sqrts) ]
      Rlabels = [ "R", "RV", "RF", "Rsqrts", "RVsqrts", "RFsqrts" ]
      for varname,varlabel in zip(Rnames,Rlabels):
         var = ROOT.RooRealVar(varname, varname, 1., 0., 400.)
         var.setVal(1)
         var.setBins(100)
         self.rrvars[varlabel]=var

      varname = "GGsm"
      var = ROOT.RooRealVar(varname, varname, 1., 0., 50.)
      var.setBins(500)
      self.rrvars["GHratio"]=var
      varname = "CMS_zz4l_GHrefval"
      var = ROOT.RooConstVar(varname, varname, self.GHrefval)
      self.rrvars["GHrefval"]=var
      varname = "kbkg_gg"
      var = ROOT.RooRealVar(varname, varname, 1., 0., 2.)
      var.setBins(200)
      self.rrvars[varname]=var
      varname = "kbkg_VBF"
      var = ROOT.RooRealVar(varname, varname, 1., 0., 2.)
      var.setBins(200)
      self.rrvars[varname]=var

      varname = "CMS_zz4l_fai1"
      var = ROOT.RooRealVar(varname, varname, 0., -1., 1.)
      var.setBins(200)
      self.rrvars["fai1"]=var
      varname = "CMS_zz4l_phiai1"
      var = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      var.setBins(200)
      self.rrvars["phiai1"]=var
      varname = "CMS_zz4l_fai2"
      var = ROOT.RooRealVar(varname, varname, 0., -1., 1.)
      var.setBins(200)
      self.rrvars["fai2"]=var
      varname = "CMS_zz4l_phiai2"
      var = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      var.setBins(200)
      self.rrvars["phiai2"]=var

      varname = "phia1"
      var = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      var.setBins(200)
      self.rrvars[varname]=var
      varname = "phia1_SB_gg"
      var = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      var.setBins(200)
      self.rrvars[varname]=var
      varname = "phia1_SB_VBF"
      var = ROOT.RooRealVar(varname, varname, 0., -math.pi, math.pi)
      var.setBins(200)
      self.rrvars[varname]=var

   # BSI+AC FORMULAE
      # Bare formula strings
      self.ggSigFormula_list = []
      self.ggInterfFormula_list = []
      self.VBFSigFormula_list = []
      self.VBFInterfFormula_list = []

      # RooFormulaVars
      # Version with no mu
      self.ggSigRFV_noMu_list = []
      self.ggInterfRFV_noMu_list = []
      self.VBFSigRFV_noMu_list = []
      self.VBFInterfRFV_noMu_list = []
      # Version with mu
      self.ggSigRFV_list = []
      self.ggInterfRFV_list = []
      self.VBFSigRFV_list = []
      self.VBFInterfRFV_list = []
      # Version with forced on-shell mu
      self.ggSigRFV_onshellMu_list = []
      self.ggInterfRFV_onshellMu_list = []
      self.VBFSigRFV_onshellMu_list = []
      self.VBFInterfRFV_onshellMu_list = []

      self.makeRFVs_BSI()


   def makeRFVs_BSI(self):
      self.rrvars["one"] = ROOT.RooConstVar("VarOne","VarOne",1.0)

   # Construct muF/V
      if self.GHmodel==1:
         self.GHmodelName="OnshellR_times_GGsm"
         self.rrvars["muF"] = ROOT.RooFormulaVar("muF_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3*@4", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"], self.rrvars["GHratio"]))
         self.rrvars["muV"] = ROOT.RooFormulaVar("muV_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3*@4", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"], self.rrvars["GHratio"]))
         self.rrvars["muF_onshell"] = ROOT.RooFormulaVar("muF_OnshellR_only_{0:.0f}TeV".format(self.sqrts), "@0*@1*@2*@3", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"]))
         self.rrvars["muV_onshell"] = ROOT.RooFormulaVar("muV_OnshellR_only_{0:.0f}TeV".format(self.sqrts), "@0*@1*@2*@3", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"]))
      elif self.GHmodel==-1:
         self.GHmodelName="OnshellR_over_GGsm"
         self.rrvars["muF"] = ROOT.RooFormulaVar("muF_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3/@4", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"], self.rrvars["GHratio"]))
         self.rrvars["muV"] = ROOT.RooFormulaVar("muV_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3/@4", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"], self.rrvars["GHratio"]))
         self.rrvars["muF_onshell"] = self.rrvars["muF"]
         self.rrvars["muV_onshell"] = self.rrvars["muV"]
      elif self.GHmodel==2:
         self.GHmodelName="OnshellR_times_GGsm_times_GHref"
         self.rrvars["muF"] = ROOT.RooFormulaVar("muF_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3*@4*@5", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"], self.rrvars["GHratio"],self.rrvars["GHrefval"]))
         self.rrvars["muV"] = ROOT.RooFormulaVar("muV_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3*@4*@5", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"], self.rrvars["GHratio"],self.rrvars["GHrefval"]))
         self.rrvars["muF_onshell"] = ROOT.RooFormulaVar("muF_OnshellR_only_{0:.0f}TeV".format(self.sqrts), "@0*@1*@2*@3", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"]))
         self.rrvars["muV_onshell"] = ROOT.RooFormulaVar("muV_OnshellR_only_{0:.0f}TeV".format(self.sqrts), "@0*@1*@2*@3", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"]))
      elif self.GHmodel==-2:
         self.GHmodelName="OnshellR_over_GGsm_over_GHref"
         self.rrvars["muF"] = ROOT.RooFormulaVar("muF_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3/(@4*@5)", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"], self.rrvars["GHratio"],self.rrvars["GHrefval"]))
         self.rrvars["muV"] = ROOT.RooFormulaVar("muV_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3/(@4*@5)", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"], self.rrvars["GHratio"],self.rrvars["GHrefval"]))
         self.rrvars["muF_onshell"] = self.rrvars["muF"]
         self.rrvars["muV_onshell"] = self.rrvars["muV"]
      else:
         self.GHmodelName="OnshellR_only"
         self.rrvars["muF"] = ROOT.RooFormulaVar("muF_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RF"], self.rrvars["Rsqrts"],self.rrvars["RFsqrts"]))
         self.rrvars["muV"] = ROOT.RooFormulaVar("muV_{0}_{1:.0f}TeV".format(self.GHmodelName,self.sqrts), "@0*@1*@2*@3", ROOT.RooArgList(self.rrvars["R"],self.rrvars["RV"], self.rrvars["Rsqrts"],self.rrvars["RVsqrts"]))
         self.rrvars["muF_onshell"] = self.rrvars["muF"]
         self.rrvars["muV_onshell"] = self.rrvars["muV"]

      # Construct phia1 variables
      self.rrvars["phia1_gg"] = ROOT.RooFormulaVar("phia1_gg", "@0+@1", ROOT.RooArgList(self.rrvars["phia1"],self.rrvars["phia1_SB_gg"]))
      self.rrvars["phia1_VBF"] = ROOT.RooFormulaVar("phia1_VBF", "@0+@1/2.", ROOT.RooArgList(self.rrvars["phia1"],self.rrvars["phia1_SB_VBF"]))

      # ggVV FORMULAE
      # In signals and interferences, @0==muF/V; additionally in interferences, @1==kbkg_gg/VBF
      if self.anomCoupl == 1: # Full parameterization with fai1=[-1, 1], and phases phiai1 and phia1_gg, phia1_VBF (phia1 are different since the bkg phase could be different)
         # 0-3 are templates, @0==fai1, @1==phiai1
         self.ggSigFormula_list.append("(1-abs(@0))")
         self.ggSigFormula_list.append("sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*cos(@1)")
         self.ggSigFormula_list.append("sign(@0)*sqrt(abs(@0)*(1-abs(@0)))*sin(@1)")
         self.ggSigFormula_list.append("abs(@0)")
         # 0-3 are templates, @0==kbkg_gg, @1==fai1, @2==phiai1, @3==phia1 (gg)
         self.ggInterfFormula_list.append("sqrt(@0)*sqrt(1-abs(@1))*cos(@3)")
         self.ggInterfFormula_list.append("sqrt(@0)*sqrt(1-abs(@1))*sin(@3)")
         self.ggInterfFormula_list.append("sqrt(@0)*sign(@1)*sqrt(abs(@1))*cos(@2+@3)")
         self.ggInterfFormula_list.append("sqrt(@0)*sign(@1)*sqrt(abs(@1))*sin(@2+@3)")
      elif self.anomCoupl == 2: # No phases, just fai1=[-1, 1]
         # 0-2 are templates, @0==fai1
         self.ggSigFormula_list.append("(1-abs(@0))")
         self.ggSigFormula_list.append("sign(@0)*sqrt(abs(@0)*(1-abs(@0)))")
         self.ggSigFormula_list.append("abs(@0)")
         # 0-1 are templates, @0==kbkg_gg, @1==fai1
         self.ggInterfFormula_list.append("sqrt(@0)*sqrt(1-abs(@1))")
         self.ggInterfFormula_list.append("sqrt(@0)*sign(@1)*sqrt(abs(@1))")
      else: # No ai1 dependence
         self.ggSigFormula_list.append("1")
         self.ggInterfFormula_list.append("sqrt(@0)") # @0==kbkg_gg

      for irfv in range(0,len(self.ggSigFormula_list)):
         rfvname_noMu = "HVV_Sig_ai1_{0:.0f}_noMu_Coef".format(irfv)
         rfvname = "HVV_Sig_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muF"].GetName())
         rfvname_onshellMu = "HVV_Sig_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muF_onshell"].GetName())
         rfvargs = ROOT.RooArgList()
         if self.anomCoupl == 1:
            rfvargs.add(self.rrvars["fai1"])
            rfvargs.add(self.rrvars["phiai1"])
         elif self.anomCoupl == 2:
            rfvargs.add(self.rrvars["fai1"])
         seg_rfv = None
         if len(self.ggSigFormula_list)>1:
            seg_rfv = ROOT.RooFormulaVar( rfvname_noMu , self.ggSigFormula_list[irfv] , rfvargs )
         else:
            seg_rfv = self.rrvars["one"]
         seg_rfv2 = ROOT.RooFormulaVar( rfvname , "abs(@0)*@1" , ROOT.RooArgList(self.rrvars["muF"] , seg_rfv) )
         seg_rfv2_onshellMu = None
         if rfvname == rfvname_onshellMu:
            seg_rfv2_onshellMu = seg_rfv2
         else:
            seg_rfv2_onshellMu = ROOT.RooFormulaVar( rfvname_onshellMu , "abs(@0)*@1" , ROOT.RooArgList(self.rrvars["muF_onshell"] , seg_rfv) )
         self.ggSigRFV_noMu_list.append(seg_rfv)
         self.ggSigRFV_list.append(seg_rfv2)
         self.ggSigRFV_onshellMu_list.append(seg_rfv2_onshellMu)

      for irfv in range(0,len(self.ggInterfFormula_list)):
         rfvname_noMu = "HVV_Interf_ai1_{0:.0f}_noMu_Coef".format(irfv)
         rfvname = "HVV_Interf_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muF"].GetName())
         rfvname_onshellMu = "HVV_Interf_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muF_onshell"].GetName())
         rfvargs = ROOT.RooArgList()
         rfvargs.add(self.rrvars["kbkg_gg"])
         if self.anomCoupl == 1:
            rfvargs.add(self.rrvars["fai1"])
            rfvargs.add(self.rrvars["phiai1"])
            rfvargs.add(self.rrvars["phia1_gg"])
         elif self.anomCoupl == 2:
            rfvargs.add(self.rrvars["fai1"])
         seg_rfv = ROOT.RooFormulaVar( rfvname_noMu , self.ggInterfFormula_list[irfv] , rfvargs )
         seg_rfv2 = ROOT.RooFormulaVar( rfvname , "sign(@0)*sqrt(abs(@0))*@1" , ROOT.RooArgList(self.rrvars["muF"] , seg_rfv) )
         seg_rfv2_onshellMu = None
         if rfvname == rfvname_onshellMu:
            seg_rfv2_onshellMu = seg_rfv2
         else:
            seg_rfv2_onshellMu = ROOT.RooFormulaVar( rfvname_onshellMu , "sign(@0)*sqrt(abs(@0))*@1" , ROOT.RooArgList(self.rrvars["muF_onshell"] , seg_rfv) )
         self.ggInterfRFV_noMu_list.append(seg_rfv)
         self.ggInterfRFV_list.append(seg_rfv2)
         self.ggInterfRFV_onshellMu_list.append(seg_rfv2_onshellMu)


      # vvVV FORMULAE
      # In signals and interferences, @0==muV; additionally in interferences, @0==kbkg_VBF
      if self.anomCoupl == 1: # Full parameterization with fai1=[-1, 1], and phases phiai1 and phia1_gg, phia1_VBF (phia1 are different since the bkg phase could be different)
         # 0-8 are templates, @0==fai1, @1==phiai1
         self.VBFSigFormula_list.append("pow((1-abs(@0)),2)")
         self.VBFSigFormula_list.append("sign(@0)*sqrt(abs(@0))*pow(sqrt(1-abs(@0)),3)*cos(@1)")
         self.VBFSigFormula_list.append("sign(@0)*sqrt(abs(@0))*pow(sqrt(1-abs(@0)),3)*sin(@1)")
         self.VBFSigFormula_list.append("abs(@0)*(1-abs(@0))")
         self.VBFSigFormula_list.append("abs(@0)*(1-abs(@0))*cos(2*@1)")
         self.VBFSigFormula_list.append("abs(@0)*(1-abs(@0))*sin(2*@1)")
         self.VBFSigFormula_list.append("sign(@0)*pow(sqrt(abs(@0)),3)*sqrt(1-abs(@0))*cos(@1)")
         self.VBFSigFormula_list.append("sign(@0)*pow(sqrt(abs(@0)),3)*sqrt(1-abs(@0))*sin(@1)")
         self.VBFSigFormula_list.append("pow(@0,2)")
         # 0-5 are templates, @0==kbkg_VBF, @1==fai1, @2==phiai1, @3==phia1 (VBF)
         self.VBFInterfFormula_list.append("sqrt(@0)*(1-abs(@1))*cos(2*@3)")
         self.VBFInterfFormula_list.append("sqrt(@0)*(1-abs(@1))*sin(2*@3)")
         self.VBFInterfFormula_list.append("sqrt(@0)*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))*cos(@2+2*@3)")
         self.VBFInterfFormula_list.append("sqrt(@0)*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))*sin(@2+2*@3)")
         self.VBFInterfFormula_list.append("sqrt(@0)*abs(@1)*cos(2*(@2+@3))")
         self.VBFInterfFormula_list.append("sqrt(@0)*abs(@1)*sin(2*(@2+@3))")
      elif self.anomCoupl == 2: # No phases, just fai1=[-1, 1]
         # 0-4 are templates, @0==fai1
         self.VBFSigFormula_list.append("pow((1-abs(@0)),2)")
         self.VBFSigFormula_list.append("sign(@0)*sqrt(abs(@0))*pow(sqrt(1-abs(@0)),3)")
         self.VBFSigFormula_list.append("abs(@0)*(1-abs(@0))")
         self.VBFSigFormula_list.append("sign(@0)*pow(sqrt(abs(@0)),3)*sqrt(1-abs(@0))")
         self.VBFSigFormula_list.append("pow(@0,2)")
         # 0-2 are templates, @0==kbkg_VBF, @1==fai1
         self.VBFInterfFormula_list.append("sqrt(@0)*(1-abs(@1))")
         self.VBFInterfFormula_list.append("sqrt(@0)*sign(@1)*sqrt(abs(@1)*(1-abs(@1)))")
         self.VBFInterfFormula_list.append("sqrt(@0)*abs(@1)")
      else: # No ai1 dependence
         self.VBFSigFormula_list.append("1")
         self.VBFInterfFormula_list.append("sqrt(@0)") # @0==kbkg_VBF

      for irfv in range(0,len(self.VBFSigFormula_list)):
         rfvname_noMu = "vvHVV_Sig_ai1_{0:.0f}_noMu_Coef".format(irfv)
         rfvname = "vvHVV_Sig_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muV"].GetName())
         rfvname_onshellMu = "vvHVV_Sig_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muV_onshell"].GetName())
         rfvargs = ROOT.RooArgList()
         if self.anomCoupl == 1:
            rfvargs.add(self.rrvars["fai1"])
            rfvargs.add(self.rrvars["phiai1"])
         elif self.anomCoupl == 2:
            rfvargs.add(self.rrvars["fai1"])
         if len(self.ggSigFormula_list)>1:
            seg_rfv = ROOT.RooFormulaVar( rfvname_noMu , self.VBFSigFormula_list[irfv] , rfvargs )
         else:
            seg_rfv = self.rrvars["one"]
         seg_rfv = ROOT.RooFormulaVar( rfvname_noMu , self.VBFSigFormula_list[irfv] , rfvargs )
         seg_rfv2 = ROOT.RooFormulaVar( rfvname , "abs(@0)*@1" , ROOT.RooArgList(self.rrvars["muV"] , seg_rfv) )
         seg_rfv2_onshellMu = None
         if rfvname == rfvname_onshellMu:
            seg_rfv2_onshellMu = seg_rfv2
         else:
            seg_rfv2_onshellMu = ROOT.RooFormulaVar( rfvname_onshellMu , "abs(@0)*@1" , ROOT.RooArgList(self.rrvars["muV_onshell"] , seg_rfv) )
         self.VBFSigRFV_noMu_list.append(seg_rfv)
         self.VBFSigRFV_list.append(seg_rfv2)
         self.VBFSigRFV_onshellMu_list.append(seg_rfv2_onshellMu)

      for irfv in range(0,len(self.VBFInterfFormula_list)):
         rfvname_noMu = "vvHVV_Interf_ai1_{0:.0f}_noMu_Coef".format(irfv)
         rfvname = "vvHVV_Interf_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muV"].GetName())
         rfvname_onshellMu = "vvHVV_Interf_ai1_{0:.0f}_{1}_Coef".format(irfv,self.rrvars["muV_onshell"].GetName())
         rfvargs = ROOT.RooArgList()
         rfvargs.add(self.rrvars["kbkg_VBF"])
         if self.anomCoupl == 1:
            rfvargs.add(self.rrvars["fai1"])
            rfvargs.add(self.rrvars["phiai1"])
            rfvargs.add(self.rrvars["phia1_VBF"])
         elif self.anomCoupl == 2:
            rfvargs.add(self.rrvars["fai1"])
         seg_rfv = ROOT.RooFormulaVar( rfvname_noMu , self.VBFInterfFormula_list[irfv] , rfvargs )
         seg_rfv2 = ROOT.RooFormulaVar( rfvname , "sign(@0)*sqrt(abs(@0))*@1" , ROOT.RooArgList(self.rrvars["muV"] , seg_rfv) )
         seg_rfv2_onshellMu = None
         if rfvname == rfvname_onshellMu:
            seg_rfv2_onshellMu = seg_rfv2
         else:
            seg_rfv2_onshellMu = ROOT.RooFormulaVar( rfvname_onshellMu , "sign(@0)*sqrt(abs(@0))*@1" , ROOT.RooArgList(self.rrvars["muV_onshell"] , seg_rfv) )
         self.VBFInterfRFV_noMu_list.append(seg_rfv)
         self.VBFInterfRFV_list.append(seg_rfv2)
         self.VBFInterfRFV_onshellMu_list.append(seg_rfv2_onshellMu)
