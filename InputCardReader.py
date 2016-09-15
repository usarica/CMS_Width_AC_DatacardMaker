#!/usr/bin/python
import os
import re
import math
import collections
#from ROOT import *
from array import array

## ---------------------------------------------------------------
## card reader class
## ---------------------------------------------------------------

class InputCardReader:

   def __init__(self, inputPath):

      if not os.path.exists(inputPath):
         raise RuntimeError, "File {0} does not exist!!!".format(inputPath)

      # input file
      self.theInput = inputPath
      # decay channel
      self.decayChan = None
      # lumi
      self.lumi = None
      # sqrts
      self.sqrts = None

      # list of
      # [ channel name, rate, lumi, (int) iBkg ]
      self.channels = []

      # list of
      # [ parameter name, value ]
      self.parameters = []

      # list of either
      # [ systematics name, systematics type=lnN, [ [channel, 1+sigma] ] ]
      # or
      # [ systematics name, systematics type=lnN, [ [channel, 1+sigma, 1-sigma] ] ]
      # or
      # [ systematics name, systematics type=param, [central value, 1 sigma error, (optional) parameter minimum, (optional) parameter maximum] ]
      # or
      # [ systematics name, systematics type=template, [ [ channel, systematics up appendix name, systematics dn appendix name ] ]
      self.systematics = []


   def parseBoolString(self,theString):
      return theString[0].upper()=='T'

   def isGoodEntry(self, var):
      if (var is None):
         return False
      elif (var == []):
         return False
      else:
         return True

   def readInputs(self):
      for line in open(self.theInput,'r'):
         f = line.split()

         if len(f) < 1: continue
         if f[0].startswith("#"): continue

         if f[0].lower().startswith("lumi"):
            self.lumi = float(f[1])

         if f[0].lower().startswith("sqrts"):
            self.sqrts = float(f[1])

         if f[0].lower().startswith("decay"):
            if f[1] == "4mu": self.decayChan = 1
            elif f[1] == "4e": self.decayChan = 2
            elif f[1] == "2e2mu": self.decayChan = 3
            elif f[1] == "2mu2e": self.decayChan = 4
            else : raise RuntimeError, "Unknown decay channel {0}, choices are 4mu, 4e, 2e2mu or 2mu2e".format(f[1])

         if f[0].lower().startswith("channel"):
            channame = f[1]
            chanrate = 1.0
            chanlumi = -1.0
            iBkg = 0
            if len(f)>2:
               chanrate = float(f[2])
               if len(f)>3:
                  chanlumi = float(f[3])
            chanlist = [ channame, chanrate, chanlumi, iBkg ]
            self.channels.append(chanlist)

         if f[0].lower().startswith("parameter"):
            parname = f[1]
            parvalue = float(f[2])
            parlist = [ parname, parvalue ]
            self.parameters.append(parlist)

         if f[0].lower().startswith("systematics"):
            parname = f[1]
            partype = f[2]
            parconfig = []
            if partype == "lnN":
               if len(f)<4:
                  raise RuntimeError, "{0} uncertainty for systematic {1} is not given any process!".format(partype, parname)
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":")
                  for itc in range(1,len(tmpconfig)): # convert strings to floats
                     tmpconfig[itc] = float(tmpconfig[itc])
                  parconfig.append(tmpconfig)
            elif partype == "param":
               if len(f)!=4:
                  raise RuntimeError, "{0} uncertainty for systematic {1} has to consist of 4 whitespace-separated string!".format(partype, parname)
               tmpconfig = f[4].split(":")
               for itc in range(0,len(tmpconfig)): # convert strings to floats
                  tmpconfig[itc] = float(tmpconfig[itc])
               parconfig = tmpconfig
            elif partype.lower() == "template":
               if len(f)<4:
                  raise RuntimeError, "{0} uncertainty name strings for systematic {1} is not given any process!".format(partype, parname)
               for ic in range(3,len(f)):
                  tmpconfig = f[ic].split(":") # Leave strings as strings
                  if ((tmpconfig[0] == "") or (tmpconfig[1] == "") or (tmpconfig[2] == "")):
                     raise RuntimeError, "{0} uncertainty does not specify any process or template appendix names!".format(parname, tmpconfig[0])
                  parconfig.append(tmpconfig)

            parlist = [ parname, partype, parconfig ]
            self.parameters.append(parlist)

    def getInputs(self):

        dict = {}

        if not self.isGoodEntry(self.sqrts): raise RuntimeError, "{0} is not set. Check inputs!".format("sqrts")
        if not self.isGoodEntry(self.lumi): raise RuntimeError, "{0} is not set. Check inputs!".format("lumi")
        if not self.isGoodEntry(self.channels): raise RuntimeError, "{0} is empty. Check inputs!".format("channels")
        if not self.isGoodEntry(self.parameters): raise RuntimeError, "{0} is empty. Check inputs!".format("parameters")
        if not self.isGoodEntry(self.systematics): raise RuntimeError, "{0} is empty. Check inputs!".format("systematics")



        dict['decayChannel'] = int(self.decayChan)
        dict['model'] = str(self.model)
        dict['lumi'] = float(self.lumi)
        dict['sqrts'] = float(self.sqrts)

        dict['all'] = self.all_chan
        dict['ggH'] = self.ggH_chan
        dict['qqH'] = self.qqH_chan
        dict['WH'] = self.WH_chan
        dict['ZH'] = self.ZH_chan
        dict['ttH'] = self.ttH_chan
        dict['ggH_SM'] = self.ggH_chan
        dict['qqH_SM'] = self.qqH_chan
        dict['WH_SM'] = self.WH_chan
        dict['ZH_SM'] = self.ZH_chan
        dict['ttH_SM'] = self.ttH_chan
        dict['qqZZ'] = self.qqZZ_chan
        dict['ggZZ'] = self.ggZZ_chan
        dict['ggZZ_signal'] = self.ggZZ_signal_chan
        dict['ggZZbkg'] = self.ggZZ_bkg_chan
        dict['ggZZ_interf'] = self.ggZZ_interf_chan
        dict['VBF_offshell'] = self.VBF_offshell_chan
        dict['zjets'] = self.zjets_chan
        dict['ttbar'] = self.ttbar_chan
        dict['zbb'] = self.zbb_chan

        dict['qqZZ_rate'] = self.qqZZ_rate
        dict['ggZZ_rate'] = self.ggZZ_rate
        dict['ggZZ_signal_rate'] = self.ggZZ_signal_rate
        dict['ggZZ_bkg_rate'] = self.ggZZ_bkg_rate
        dict['ggZZ_interf_rate'] = self.ggZZ_interf_rate
        dict['zjets_rate'] = self.zjets_rate
        dict['ttbar_rate'] = self.ttbar_rate
        dict['zbb_rate'] = self.zbb_rate
        dict['catscale_ggzz'] = self.djetscale_ggzz
        dict['catscale_vbf_offshell'] = self.djetscale_vbf_offshell
        dict['catscale_bkg_qqzz'] = self.djetscale_bkg_qqzz
        dict['catscale_bkg_zjets'] = self.djetscale_bkg_zjets

        dict['qqZZ_lumi'] = float(self.qqZZ_lumi)
        dict['ggZZ_lumi'] = float(self.ggZZ_lumi)
        dict['zjets_lumi'] = float(self.zjets_lumi)
        dict['ttbar_lumi'] = float(self.ttbar_lumi)
        dict['zbb_lumi'] = float(self.zbb_lumi)

        dict['n_CB_shape'] = self.n_CB_shape
        dict['alpha_CB_shape'] = self.alpha_CB_shape
        dict['n2_CB_shape'] = self.n2_CB_shape
        dict['alpha2_CB_shape'] = self.alpha2_CB_shape
        dict['mean_CB_shape'] = self.mean_CB_shape
        dict['sigma_CB_shape'] = self.sigma_CB_shape

        dict['useHighMassReweightedShapes'] = self.useHighMassReweightedShapes
        dict['n_CB_shape_HM'] = self.n_CB_shape_HM
        dict['alpha_CB_shape_HM'] = self.alpha_CB_shape_HM
        dict['n2_CB_shape_HM'] = self.n2_CB_shape_HM
        dict['alpha2_CB_shape_HM'] = self.alpha2_CB_shape_HM
        dict['mean_CB_shape_HM'] = self.mean_CB_shape_HM
        dict['sigma_CB_shape_HM'] = self.sigma_CB_shape_HM
        dict['gamma_BW_shape_HM'] = self.gamma_BW_shape_HM

        dict['sigEff_a1'] = float(self.sigeff_a1)
        dict['sigEff_a2'] = float(self.sigeff_a2)
        dict['sigEff_a3'] = float(self.sigeff_a3)
        dict['sigEff_a4'] = float(self.sigeff_a4)
        dict['sigEff_b1'] = float(self.sigeff_b1)
        dict['sigEff_b2'] = float(self.sigeff_b2)
        dict['sigEff_b3'] = float(self.sigeff_b3)
        dict['sigEff_g1'] = float(self.sigeff_g1)
        dict['sigEff_g2'] = float(self.sigeff_g2)
        dict['sigEff_g3'] = float(self.sigeff_g3)

        dict['sigEff_qqHa1'] = float(self.sigeff_qqHa1)
        dict['sigEff_qqHa2'] = float(self.sigeff_qqHa2)
        dict['sigEff_qqHa3'] = float(self.sigeff_qqHa3)
        dict['sigEff_qqHa4'] = float(self.sigeff_qqHa4)
        dict['sigEff_qqHb1'] = float(self.sigeff_qqHb1)
        dict['sigEff_qqHb2'] = float(self.sigeff_qqHb2)
        dict['sigEff_qqHb3'] = float(self.sigeff_qqHb3)
        dict['sigEff_qqHg1'] = float(self.sigeff_qqHg1)
        dict['sigEff_qqHg2'] = float(self.sigeff_qqHg2)
        dict['sigEff_qqHg3'] = float(self.sigeff_qqHg3)

        dict['sigEff_WHa1'] = float(self.sigeff_WHa1)
        dict['sigEff_WHa2'] = float(self.sigeff_WHa2)
        dict['sigEff_WHa3'] = float(self.sigeff_WHa3)
        dict['sigEff_WHa4'] = float(self.sigeff_WHa4)
        dict['sigEff_WHb1'] = float(self.sigeff_WHb1)
        dict['sigEff_WHb2'] = float(self.sigeff_WHb2)
        dict['sigEff_WHb3'] = float(self.sigeff_WHb3)
        dict['sigEff_WHg1'] = float(self.sigeff_WHg1)
        dict['sigEff_WHg2'] = float(self.sigeff_WHg2)
        dict['sigEff_WHg3'] = float(self.sigeff_WHg3)

        dict['sigEff_ZHa1'] = float(self.sigeff_ZHa1)
        dict['sigEff_ZHa2'] = float(self.sigeff_ZHa2)
        dict['sigEff_ZHa3'] = float(self.sigeff_ZHa3)
        dict['sigEff_ZHa4'] = float(self.sigeff_ZHa4)
        dict['sigEff_ZHb1'] = float(self.sigeff_ZHb1)
        dict['sigEff_ZHb2'] = float(self.sigeff_ZHb2)
        dict['sigEff_ZHb3'] = float(self.sigeff_ZHb3)
        dict['sigEff_ZHg1'] = float(self.sigeff_ZHg1)
        dict['sigEff_ZHg2'] = float(self.sigeff_ZHg2)
        dict['sigEff_ZHg3'] = float(self.sigeff_ZHg3)

        dict['sigEff_ttHa1'] = float(self.sigeff_ttHa1)
        dict['sigEff_ttHa2'] = float(self.sigeff_ttHa2)
        dict['sigEff_ttHa3'] = float(self.sigeff_ttHa3)
        dict['sigEff_ttHa4'] = float(self.sigeff_ttHa4)
        dict['sigEff_ttHb1'] = float(self.sigeff_ttHb1)
        dict['sigEff_ttHb2'] = float(self.sigeff_ttHb2)
        dict['sigEff_ttHb3'] = float(self.sigeff_ttHb3)
        dict['sigEff_ttHg1'] = float(self.sigeff_ttHg1)
        dict['sigEff_ttHg2'] = float(self.sigeff_ttHg2)
        dict['sigEff_ttHg3'] = float(self.sigeff_ttHg3)

        dict['tagged_ggH_ratio'] = self.tagged_ggH_ratio
        dict['tagged_qqH_ratio'] = self.tagged_qqH_ratio
        dict['tagged_WH_ratio'] = self.tagged_WH_ratio
        dict['tagged_ZH_ratio'] = self.tagged_ZH_ratio
        dict['tagged_ttH_ratio'] = self.tagged_ttH_ratio
        dict['QCD_scale_ggH_2j_sys'] = float(self.QCD_scale_ggH_2j_sys)
        dict['QCD_scale_qqH_2j_sys'] = float(self.QCD_scale_qqH_2j_sys)
        #dict['QCD_scale_qqZZ_2j_sys'] = float(self.QCD_scale_qqZZ_2j_sys)

        dict['dijetRatio'] = self.dijetRatio

        dict['qqZZshape_a0'] = float(self.qqZZshape_a0)
        dict['qqZZshape_a1'] = float(self.qqZZshape_a1)
        dict['qqZZshape_a2'] = float(self.qqZZshape_a2)
        dict['qqZZshape_a3'] = float(self.qqZZshape_a3)
        dict['qqZZshape_a4'] = float(self.qqZZshape_a4)
        dict['qqZZshape_a5'] = float(self.qqZZshape_a5)
        dict['qqZZshape_a6'] = float(self.qqZZshape_a6)
        dict['qqZZshape_a7'] = float(self.qqZZshape_a7)
        dict['qqZZshape_a8'] = float(self.qqZZshape_a8)
        dict['qqZZshape_a9'] = float(self.qqZZshape_a9)
        dict['qqZZshape_a10'] = float(self.qqZZshape_a10)
        dict['qqZZshape_a11'] = float(self.qqZZshape_a11)
        dict['qqZZshape_a12'] = float(self.qqZZshape_a12)
        dict['qqZZshape_a13'] = float(self.qqZZshape_a13)

        dict['ggZZshape_a0'] = float(self.ggZZshape_a0)
        dict['ggZZshape_a1'] = float(self.ggZZshape_a1)
        dict['ggZZshape_a2'] = float(self.ggZZshape_a2)
        dict['ggZZshape_a3'] = float(self.ggZZshape_a3)
        dict['ggZZshape_a4'] = float(self.ggZZshape_a4)
        dict['ggZZshape_a5'] = float(self.ggZZshape_a5)
        dict['ggZZshape_a6'] = float(self.ggZZshape_a6)
        dict['ggZZshape_a7'] = float(self.ggZZshape_a7)
        dict['ggZZshape_a8'] = float(self.ggZZshape_a8)
        dict['ggZZshape_a9'] = float(self.ggZZshape_a9)

        dict['zjetsShape_mean_3P1F'] = float(self.zjetsShape_mean_3P1F)
        dict['zjetsShape_sigma_3P1F'] = float(self.zjetsShape_sigma_3P1F)
        dict['zjetsShape_norm_3P1F'] = float(self.zjetsShape_norm_3P1F)

        dict['zjetsShape_mean_2P2F'] = float(self.zjetsShape_mean_2P2F)
        dict['zjetsShape_sigma_2P2F'] = float(self.zjetsShape_sigma_2P2F)
        dict['zjetsShape_norm_2P2F'] = float(self.zjetsShape_norm_2P2F)
        dict['zjetsShape_pol0_2P2F'] = float(self.zjetsShape_pol0_2P2F)
        dict['zjetsShape_pol1_2P2F'] = float(self.zjetsShape_pol1_2P2F)

        dict['zjetsShape_mean_2P2F_2e2mu'] = float(self.zjetsShape_mean_2P2F_2e2mu)
        dict['zjetsShape_sigma_2P2F_2e2mu'] = float(self.zjetsShape_sigma_2P2F_2e2mu)
        dict['zjetsShape_norm_2P2F_2e2mu'] = float(self.zjetsShape_norm_2P2F_2e2mu)

        dict['zjetsKappaLow'] = float(self.zjetsKappaLow)
        dict['zjetsKappaHigh'] = float(self.zjetsKappaHigh)


        dict['lumiUnc'] = self.lumiUnc
        dict['muonFullUnc'] = float(self.muonFullUnc)
        dict['muonFullUnc_HM'] = float(self.muonFullUnc_HM)
        dict['muonFullCutoff'] = float(self.muonFullCutoff)
        dict['elecFullUnc'] = float(self.elecFullUnc)
        dict['elecFullUnc_HM'] = float(self.elecFullUnc_HM)
        dict['elecFullCutoff'] = float(self.elecFullCutoff)

        dict['muonTrigUnc'] = float(self.muonTrigUnc)
        dict['muonTrigUnc_HM'] = float(self.muonTrigUnc_HM)
        dict['muonTrigCutoff'] = float(self.muonTrigCutoff)
        dict['elecTrigUnc'] = float(self.elecTrigUnc)
        dict['elecTrigUnc_HM'] = float(self.elecTrigUnc_HM)
        dict['elecTrigCutoff'] = float(self.elecTrigCutoff)

        dict['useLumiUnc'] = self.useLumiUnc
        dict['usePdf_gg'] = self.usePdf_gg
        dict['usePdf_qqbar'] = self.usePdf_qqbar
        dict['usePdf_hzz4l_accept'] = self.usePdf_hzz4l_accept
        dict['useQCDscale_ggH'] = self.useQCDscale_ggH
        dict['useQCDscale_qqH'] = self.useQCDscale_qqH
        dict['useQCDscale_VH'] = self.useQCDscale_VH
        dict['useQCDscale_ttH'] = self.useQCDscale_ttH
        dict['useTheoryUncXS_HighMH'] = self.useTheoryUncXS_HighMH
        dict['useQCDscale_ggVV'] = self.useQCDscale_ggVV
        dict['useQCDscale_VV'] = self.useQCDscale_VV
        dict['useBRhiggs_hzz4l'] = self.useBRhiggs_hzz4l
        dict['useCMS_eff'] = self.useCMS_eff
        dict['useCMS_hzz4l_Zjets'] = self.useCMS_hzz4l_Zjets
        dict['useCMS_zz4l_bkgMELA'] = self.useCMS_zz4l_bkgMELA
        dict['useCMS_zz4l_sigMELA'] = self.useCMS_zz4l_sigMELA
        dict['useCMS_zz4l_mean'] = self.useCMS_zz4l_mean
        dict['useCMS_zz4l_sigma'] = self.useCMS_zz4l_sigma
        dict['useCMS_zz4l_n'] = self.useCMS_zz4l_n
        dict['useCMS_zz4l_gamma'] = self.useCMS_zz4l_gamma

        dict['CMS_zz4l_mean_m_sig'] = float(self.CMS_zz4l_mean_m_sig)
        dict['CMS_zz4l_sigma_m_sig'] = float(self.CMS_zz4l_sigma_m_sig)
        dict['CMS_zz4l_mean_e_sig'] = float(self.CMS_zz4l_mean_e_sig)
        dict['CMS_zz4l_sigma_e_sig'] = float(self.CMS_zz4l_sigma_e_sig)
        dict['CMS_zz4l_n_sig'] = float(self.CMS_zz4l_n_sig)
        dict['CMS_zz4l_gamma_sig'] = float(self.CMS_zz4l_gamma_sig)

        dict['doHypTest'] = self.doHypTest
        dict['altHypLabel'] = str(self.altHypLabel)

        dict['useCMS_zz4l_doVBFtest'] = self.useCMS_zz4l_doVBFtest
        dict['useCMS_zz4l_Fisher_sys'] = self.useCMS_zz4l_Fisher_sys
        dict['useCMS_zz4l_Pt_sys'] = self.useCMS_zz4l_Pt_sys


	dict['mekd_sig_a0_shape'] = self.mekd_sig_a0_shape
	dict['mekd_sig_a1_shape'] = self.mekd_sig_a1_shape
	dict['mekd_sig_a2_shape'] = self.mekd_sig_a2_shape
	dict['mekd_sig_a3_shape'] = self.mekd_sig_a3_shape
	dict['mekd_sig_a4_shape'] = self.mekd_sig_a4_shape
	dict['mekd_qqZZ_a0_shape'] = self.mekd_qqZZ_a0_shape
	dict['mekd_qqZZ_a1_shape'] = self.mekd_qqZZ_a1_shape
	dict['mekd_qqZZ_a2_shape'] = self.mekd_qqZZ_a2_shape
	dict['mekd_qqZZ_a3_shape'] = self.mekd_qqZZ_a3_shape
	dict['mekd_qqZZ_a4_shape'] = self.mekd_qqZZ_a4_shape

	dict['relerr_ggH_ld_frac'] = self.relerr_ggH_ld_frac
	dict['relerr_ggH_ld_mean'] = self.relerr_ggH_ld_mean
	dict['relerr_ggH_ld_sigma'] = self.relerr_ggH_ld_sigma
	dict['relerr_ggH_gs_mean'] = self.relerr_ggH_gs_mean
	dict['relerr_ggH_gs_sigma'] = self.relerr_ggH_gs_sigma
	dict['relerr_qqzz_ld_frac'] = self.relerr_qqzz_ld_frac
	dict['relerr_qqzz_ld_mean'] = self.relerr_qqzz_ld_mean
	dict['relerr_qqzz_ld_sigma'] = self.relerr_qqzz_ld_sigma
	dict['relerr_qqzz_gs_mean'] = self.relerr_qqzz_gs_mean
	dict['relerr_qqzz_gs_sigma'] = self.relerr_qqzz_gs_sigma
	dict['relerr_zx_ld_frac'] = self.relerr_zx_ld_frac
	dict['relerr_zx_ld_mean'] = self.relerr_zx_ld_mean
	dict['relerr_zx_ld_sigma'] = self.relerr_zx_ld_sigma
	dict['relerr_zx_gs_mean'] = self.relerr_zx_gs_mean
	dict['relerr_zx_gs_sigma'] = self.relerr_zx_gs_sigma

        return dict
