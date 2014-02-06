#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClass import *
from inputReader import *

## ------------------------------------
##  card and workspace class
## ------------------------------------

class datacardClass:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True

    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")


    # cs x br function 
    def makeXsBrFunction(self,signalProc,rrvMH):
            
        procName = "ggH"
        if(signalProc == 0): procName = "ggH" #dummy, when you sum up all the 5 chans
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"

        
        
        
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)" 

     
        
        myCSWrhf = HiggsCSandWidth()
        
        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 8905, 109.55, 1000.05)
                
        for i in range(1,8906):
            
            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0 
            if (self.channel == self.ID_2e2mu):
                BR = myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = myCSWrhf.HiggsBR(12,mHVal)

            if (signalProc == 3 or signalProc == 4 or signalProc == 5):
                #overwrite BR if VH,ttH sample
                #these samples have inclusive Z decay
                BR = myCSWrhf.HiggsBR(11,mHVal)

            if (signalProc==0):
                totXs=0
                for ch in range(1,6):
                    totXs+=myCSWrhf.HiggsCS(ch, mHVal, self.sqrts)
                histXsBr.SetBinContent(i, totXs * BR)
            else:
                histXsBr.SetBinContent(i, myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            #print '\nmakeXsBrFunction : procName=',procName,'   signalProc=',signalProc,'  mH (input)=',rrvMH.getVal(),
            #print '   CS=',myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts),'   BR=',BR
            
        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)  
        
        return rdhXsBr
    
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
    
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theis2D, theOutputDir, theInputs,theTemplateDir="templates2D", theIncludingError=False, theMEKD=False, theVBFcat=False, theUse3D=False):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = theMH
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.is2D = theis2D
        self.outputDir = theOutputDir
        self.sigMorph = theInputs['useCMS_zz4l_sigMELA']
        self.bkgMorph = theInputs['useCMS_zz4l_bkgMELA']
        self.templateDir = theTemplateDir
	self.bIncludingError=theIncludingError
	self.bMEKD = False
        self.useMEKDTemplates = theMEKD
        self.bVBF = False
        if(theInputs['useCMS_zz4l_doVBFtest']):
            self.bVBF = True
        self.VBFcat = theVBFcat
        self.Use3D = theUse3D
        self.FisherMorph = theInputs['useCMS_zz4l_Fisher_sys']
        self.PtMorph = theInputs['useCMS_zz4l_Pt_sys']
        
        FactorizedShapes = False

        self.all_chan = theInputs['all']
        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']
        self.ttbar_chan = theInputs['ttbar']
        self.zbb_chan = theInputs['zbb']
        
        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        ## ---------------- VARIABLES FOR LATER --------------- ##
        self.bUseCBnoConvolution = False
        ForXSxBR = False

        myCSW = HiggsCSandWidth()
                
        ## ----------------- WIDTH AND RANGES ----------------- ##
        self.widthHVal =  myCSW.HiggsWidth(0,self.mH)
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
        self.isHighMass = False
        if self.mH >= 390:
            if theInputs['useHighMassReweightedShapes']:
                self.isHighMass = True
            else: print "useHighMassReweightedShapes set to FALSE, using non-reweighted shapes!"

            
        print "width: ",self.widthHVal
        
        self.windowVal = max( self.widthHVal, 1.0)
        lowside = 100.0
        highside = 1000.0
        
        if (self.mH >= 275):
            lowside = 180.0
            highside = 650.0
        if (self.mH >= 350):
            lowside = 200.0
            highside = 900.0
        if (self.mH >= 500):
            lowside = 250.0
            highside = 1000.0
        if (self.mH >= 700):
            lowside = 350.0
            highside = 1400.0
        
        self.low_M = max( (self.mH - 20.*self.windowVal), lowside)
        self.high_M = min( (self.mH + 15.*self.windowVal), highside)

        #self.low_M = 100.0
        #self.high_M = 800.0
       
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"

        self.isAltSig = False
        if (theInputs['doHypTest']):
            self.isAltSig = True
            
        if self.isAltSig and not self.all_chan :
            raise RuntimeError, "You asked to prepare DC and WS for Hyp Test but you did not want to sum over all signal channels. This is forbidden. Check inputs ! (it should have already send you this error message, strange that  you are here...)"

        if (self.isAltSig and not (self.is2D==1)):
            raise RunTimeError, "Cannot perform hypothesis testing without a 2D analysis, feature not supported yet. Exiting."

        if (self.bVBF and not (self.is2D==1)):
            raise RunTimeError, "Cannot perform jet catagory separation without a 2D analysis, feature not supported. Exiting."
        if (self.bVBF and self.bIncludingError):
            raise RunTimeError, "Cannot perform jet catagory separation with EbE Errors, feature not tested yet. Exiting."
        if (self.bVBF and self.isAltSig):
            raise RunTimeError, "Cannot perform jet catagory separation with alt. hpyth., feature not supported yet. Exiting."
        if (self.bVBF and self.bMEKD):
            raise RunTimeError, "Cannot perform jet catagory separation with MEKD, feature not tested yet. Exiting."
        

        self.appendHypType = theInputs['altHypLabel']
        if self.isAltSig and self.appendHypType=="" :
            self.appendHypType = "_ALT"
            
        
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)
        systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        bins = 1000
        if(self.bUseCBnoConvolution): bins = 200
        if(self.bIncludingError): bins = 1000

        CMS_zz4l_mass_name = "CMS_zz4l_mass"
            
        CMS_zz4l_mass = ROOT.RooRealVar(CMS_zz4l_mass_name,CMS_zz4l_mass_name,self.low_M,self.high_M)    
        CMS_zz4l_mass.setBins(bins)

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)

	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB

        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        gamma_BW_d = 0.0
        
        if(self.all_chan):
            rdhXsBrFuncV_1 = self.makeXsBrFunction(0,self.MH)
        else:
            rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        if not self.bVBF:
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.sqrts)
            rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)
            
            rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.sqrts)
            rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)
            
            rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.sqrts)
            rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)
            
            rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.sqrts)
            rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)
            
            rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.sqrts)
            rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)
        else:
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}_{3}".format("ggH",self.channel,self.sqrts,self.VBFcat)
            rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)
            
            rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}_{3}".format("VBF",self.channel,self.sqrts,self.VBFcat)
            rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)
            
            rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}_{3}".format("WH",self.channel,self.sqrts,self.VBFcat)
            rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)
            
            rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}_{3}".format("ZH",self.channel,self.sqrts,self.VBFcat)
            rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)
            
            rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
            rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}_{3}".format("ttH",self.channel,self.sqrts,self.VBFcat)
            rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)
    
        ## -------- Variable Definitions -------- ##
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_e_err_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_e_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_err",float(theInputs['CMS_zz4l_mean_e_sig']),-0.99,0.99)
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_e_sig",3.0,0.0,30.0)
        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_m_err_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_m_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_m_err",float(theInputs['CMS_zz4l_mean_m_sig']),-0.99,0.99)
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
            
        
        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",1.,-10.,10.)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",2.,-10.,10.)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.low_M,self.high_M)
        name = "interf_ggH"
        #name = "CMS_zz4l_gamma_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_gamma = ROOT.RooRealVar(name,"CMS_zz4l_gamma",10.,0.001,1000.)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)
            
        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)
    
        CMS_zz4l_mean_BW.setVal( mean_BW_d )
        CMS_zz4l_gamma.setVal(0)
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_mean_e_err.setConstant(kTRUE)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_mean_m_err.setConstant(kTRUE)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
        CMS_zz4l_alpha2.setVal(0)
        CMS_zz4l_n2.setVal(0)
    
        CMS_zz4l_widthScale.setConstant(True)
        #CMS_zz4l_alpha.setConstant(True)  # also read from input file
        CMS_zz4l_mean_BW.setConstant(True)
        #CMS_zz4l_gamma_BW.setConstant(True)

        print "HEEERRRRRRRRRRRRRRRRREEEEEEE"

        print "mean_BW ", CMS_zz4l_mean_BW.getVal()
        print "gamma_BW ", CMS_zz4l_gamma.getVal()
        print "mean_e_sig ", CMS_zz4l_mean_e_sig.getVal()
        print "mean_e_err ", CMS_zz4l_mean_e_err.getVal()
        print "sigma_e ", CMS_zz4l_sigma_e_sig.getVal()
        print "mean_m_sig ",CMS_zz4l_mean_m_sig.getVal()
        print "mean_m_err ", CMS_zz4l_mean_m_err.getVal()
        print "sigma_m ", CMS_zz4l_sigma_m_sig.getVal()
        print "alpha ", CMS_zz4l_alpha.getVal()
        print "n ", CMS_zz4l_n.getVal()
        print "alpha2 ", CMS_zz4l_alpha2.getVal()
        print "n2 ", CMS_zz4l_n2.getVal()

                                                                


        ## -------------------- RooFormulaVar's -------------------- ##
        rfv_n_CB = ROOT.RooFormulaVar()
        rfv_alpha_CB = ROOT.RooFormulaVar()
        rfv_n2_CB = ROOT.RooFormulaVar()
        rfv_alpha2_CB = ROOT.RooFormulaVar()
        rfv_mean_CB = ROOT.RooFormulaVar()
        rfv_sigma_CB = ROOT.RooFormulaVar()

        if not self.bVBF:
            name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
            if self.isHighMass : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
            else : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
            
            name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
            if self.isHighMass : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape_HM'], ROOT.RooArgList(self.MH))
            else : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))
            
            name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
            #if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
            #else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
            if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")",ROOT.RooArgList(self.MH))
            else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))

            name = "CMS_zz4l_alpha2_{0:.0f}_centralValue".format(self.channel)
            if self.isHighMass : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape_HM'], ROOT.RooArgList(self.MH))
            else : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))
            
            name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        else:
            name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_{2}_centralValue".format(self.channel,self.sqrts,self.VBFcat)
            if self.isHighMass : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
            else : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
            
            name = "CMS_zz4l_alpha_{0:.0f}_{1}_centralValue".format(self.channel,self.VBFcat)
            if self.isHighMass : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape_HM'], ROOT.RooArgList(self.MH))
            else : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))
            
            name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_{2}_centralValue".format(self.channel,self.sqrts,self.VBFcat)
            #if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
            #else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
            if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")",ROOT.RooArgList(self.MH))
            else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))

            name = "CMS_zz4l_alpha2_{0:.0f}_{1}_centralValue".format(self.channel,self.VBFcat)
            if self.isHighMass : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape_HM'], ROOT.RooArgList(self.MH))
            else : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))
            
            name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_{2}_centralValue".format(self.channel,self.sqrts,self.VBFcat)
            
        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+ (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+ (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))
        

        if not self.bVBF:
            name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        else:
            name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_{2}_centralValue".format(self.channel,self.sqrts,self.VBFcat)
        
        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))

        if not self.bVBF:
            name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        else:
            name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_{2}_centralValue".format(self.channel,self.sqrts,self.VBFcat)
        rfv_gamma_BW = ROOT.RooFormulaVar(name,"("+theInputs['gamma_BW_shape_HM']+")"+"*(1+@1*0.05)",ROOT.RooArgList(self.MH,CMS_zz4l_gamma))

        if (DEBUG): print " DEBUG *********  ", theInputs['sigma_CB_shape'] 

        print "n_CB ", rfv_n_CB.getVal()
        print "alpha_CB ", rfv_alpha_CB.getVal()
        print "n2_CB ", rfv_n2_CB.getVal()
        print "alpha2_CB ", rfv_alpha2_CB.getVal()
        print "mean_CB ", rfv_mean_CB.getVal()
        print "sigma_CB ", rfv_sigma_CB.getVal()
        print "gamma_BW ", rfv_gamma_BW.getVal()    

        if not self.bVBF:
            CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))
        else:
            CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))

        print "mean_sig_NoConv ", CMS_zz4l_mean_sig_NoConv.getVal()

        # bIncludingError variables
        mrelerrVarName = "CMS_zz4l_massRelErr"
        RelErr = ROOT.RooRealVar(mrelerrVarName,mrelerrVarName,0.002,0.2)
        RelErr.setBins(30)
        CMS_zz4l_massErr = ROOT.RooFormulaVar()
        if (self.channel == self.ID_4mu) :
            CMS_zz4l_massErr = ROOT.RooFormulaVar("CMS_zz4l_massErr","@0*@1*(1+@2)",ROOT.RooArgList(CMS_zz4l_mass,RelErr,CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            CMS_zz4l_massErr = ROOT.RooFormulaVar("CMS_zz4l_massErr","@0*@1*(1+@2)",ROOT.RooArgList(CMS_zz4l_mass,RelErr,CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            CMS_zz4l_massErr = ROOT.RooFormulaVar("CMS_zz4l_massErr","@0*@1*TMath::Sqrt((1+@2)*(1+@3))",ROOT.RooArgList(CMS_zz4l_mass,RelErr,CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))
        # end bIncludingError
        
        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
    
        signalCB_ggH = ROOT.RooDoubleCB("signalCB_ggH","signalCB_ggH",CMS_zz4l_mass, self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB, self.bUseCBnoConvolution) , self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH,signalCB_ggH, 2)
        #High mass pdf
        signalBW_ggH_HM = ROOT.RooRelBWHighMass("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ggH_HM =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH_HM,signalCB_ggH, 2)
  
        
        signalCB_VBF = ROOT.RooDoubleCB("signalCB_VBF","signalCB_VBF",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_VBF = ROOT.RooRelBWUFParam("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_VBF = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF,signalCB_VBF, 2)
        #High mass pdf
        signalBW_VBF_HM = ROOT.RooRelBWHighMass("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_VBF_HM = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF_HM,signalCB_VBF, 2)
                       
        
        signalCB_WH = ROOT.RooDoubleCB("signalCB_WH","signalCB_WH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_WH = ROOT.RooRelBWUFParam("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_WH = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH,signalCB_WH, 2)
        #High mass pdf
        signalBW_WH_HM = ROOT.RooRelBWHighMass("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_WH_HM = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH_HM,signalCB_WH, 2)

        
        signalCB_ZH = ROOT.RooDoubleCB("signalCB_ZH","signalCB_ZH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ZH = ROOT.RooRelBWUFParam("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ZH = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH,signalCB_ZH, 2)
        #High mass pdf
        signalBW_ZH_HM = ROOT.RooRelBWHighMass("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ZH_HM = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH_HM,signalCB_ZH, 2)

        
        signalCB_ttH = ROOT.RooDoubleCB("signalCB_ttH","signalCB_ttH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ttH = ROOT.RooRelBWUFParam("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ttH = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH,signalCB_ttH, 2) 
        #High mass pdf
        signalBW_ttH_HM = ROOT.RooRelBWHighMass("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ttH_HM = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH_HM,signalCB_ttH, 2)
        
        
        ## Buffer fraction for cyclical behavior
        sig_ggH.setBufferFraction(0.2)
        sig_VBF.setBufferFraction(0.2)
        sig_WH.setBufferFraction(0.2)
        sig_ZH.setBufferFraction(0.2)
        sig_ttH.setBufferFraction(0.2)
        
        sig_ggH_HM.setBufferFraction(0.2)
        sig_VBF_HM.setBufferFraction(0.2)
        sig_WH_HM.setBufferFraction(0.2)
        sig_ZH_HM.setBufferFraction(0.2)
        sig_ttH_HM.setBufferFraction(0.2)


	#------------------------------------------------begin  bIncludingError 

	name = "CMS_zz4l_massErrS_ln_kappa_{0:.0f}".format(self.channel)
	rfv_EBE_sig_ln_kappa = ROOT.RooFormulaVar(name, "("+theInputs['relerr_ggH_gs_sigma']+")", ROOT.RooArgList(self.MH))
	name = "CMS_zz4l_massErrS_ln_mean_{0:.0f}".format(self.channel)
	rfv_EBE_sig_ln_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_ggH_gs_mean']+")", ROOT.RooArgList(self.MH))
	EBE_sig_ln = ROOT.RooLognormal("errLN_ggH","errLN_ggH", RelErr, rfv_EBE_sig_ln_mean, rfv_EBE_sig_ln_kappa)
	if self.channel!=1: EBE_sig_ln = ROOT.RooGaussian("errGaus_ggH","errGaus_ggH", RelErr, rfv_EBE_sig_ln_mean, rfv_EBE_sig_ln_kappa)
	name = "CMS_zz4l_massErrS_ld_sigma_{0:.0f}".format(self.channel)
	rfv_EBE_sig_ld_sigma = ROOT.RooFormulaVar(name, "("+theInputs['relerr_ggH_ld_sigma']+")", ROOT.RooArgList(self.MH))
	name = "CMS_zz4l_massErrS_ld_mean_{0:.0f}".format(self.channel)
	rfv_EBE_sig_ld_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_ggH_ld_mean']+")", ROOT.RooArgList(self.MH))
	EBE_sig_ld = ROOT.RooLandau("errLD_ggH","errLD_ggH", RelErr, rfv_EBE_sig_ld_mean, rfv_EBE_sig_ld_sigma)
	name = "CMS_zz4l_massErrS_ld_frac_{0:.0f}".format(self.channel)
	rfv_EBE_sig_frac = ROOT.RooFormulaVar(name, "("+theInputs['relerr_ggH_ld_frac']+")", ROOT.RooArgList(self.MH))
	pdfErrS = ROOT.RooAddPdf("pdfErrS","pdfErrS", EBE_sig_ld, EBE_sig_ln, rfv_EBE_sig_frac)

	name = "CMS_zz4l_massErrZZ_ln_kappa_{0:.0f}".format(self.channel)
	rfv_EBE_zz_ln_kappa = ROOT.RooFormulaVar(name, "("+theInputs['relerr_qqzz_gs_sigma']+")", ROOT.RooArgList(CMS_zz4l_mass))
	name = "CMS_zz4l_massErrZZ_ln_mean_{0:.0f}".format(self.channel)
	rfv_EBE_zz_ln_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_qqzz_gs_mean']+")", ROOT.RooArgList(CMS_zz4l_mass))
	EBE_zz_ln = ROOT.RooLognormal("errLN_qqzz","errLN_qqzz", RelErr, rfv_EBE_zz_ln_mean, rfv_EBE_zz_ln_kappa)
	if self.channel!=1: EBE_zz_ln = ROOT.RooGaussian("errGaus_qqzz","errGaus_qqzz", RelErr, rfv_EBE_zz_ln_mean, rfv_EBE_zz_ln_kappa)	
	name = "CMS_zz4l_massErrZZ_ld_sigma_{0:.0f}".format(self.channel)
	rfv_EBE_zz_ld_sigma = ROOT.RooFormulaVar(name, "("+theInputs['relerr_qqzz_ld_sigma']+")", ROOT.RooArgList(CMS_zz4l_mass))
	name = "CMS_zz4l_massErrZZ_ld_mean_{0:.0f}".format(self.channel)
	rfv_EBE_zz_ld_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_qqzz_ld_mean']+")", ROOT.RooArgList(CMS_zz4l_mass)) 
	EBE_zz_ld = ROOT.RooLandau("errLD_qqzz","errLD_qqzz", RelErr, rfv_EBE_zz_ld_mean, rfv_EBE_zz_ld_sigma)
	name = "CMS_zz4l_massErrZZ_ld_frac_{0:.0f}".format(self.channel)
	rfv_EBE_zz_frac = ROOT.RooFormulaVar(name, "("+theInputs['relerr_qqzz_ld_frac']+")", ROOT.RooArgList(self.MH)) 
	pdfErrZZ = ROOT.RooAddPdf("pdfErrZZ","pdfErrZZ", EBE_zz_ld, EBE_zz_ln, rfv_EBE_zz_frac)

	name = "CMS_zz4l_massErrZX_ln_kappa_{0:.0f}".format(self.channel)
	rfv_EBE_zx_ln_kappa = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_gs_sigma']+")", ROOT.RooArgList(CMS_zz4l_mass)) 
	name = "CMS_zz4l_massErrZX_ln_mean_{0:.0f}".format(self.channel)
	rfv_EBE_zx_ln_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_gs_mean']+")", ROOT.RooArgList(CMS_zz4l_mass)) 
	EBE_zx_ln = ROOT.RooLognormal("errLN_zx","errLN_zx", RelErr, rfv_EBE_zx_ln_mean, rfv_EBE_zx_ln_kappa)
	if self.channel!=1: EBE_zx_ln = ROOT.RooGaussian("errGaus_zx","errGaus_zx", RelErr, rfv_EBE_zx_ln_mean, rfv_EBE_zx_ln_kappa)	
	name = "CMS_zz4l_massErrZX_ld_sigma_{0:.0f}".format(self.channel)
	rfv_EBE_zx_ld_sigma = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_ld_sigma']+")", ROOT.RooArgList(CMS_zz4l_mass)) 
	name = "CMS_zz4l_massErrZX_ld_mean_{0:.0f}".format(self.channel)
	rfv_EBE_zx_ld_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_ld_mean']+")", ROOT.RooArgList(CMS_zz4l_mass)) 
	EBE_zx_ld = ROOT.RooLandau("errLD_zx","errLD_zx", RelErr, rfv_EBE_zx_ld_mean, rfv_EBE_zx_ld_sigma)	
	name = "CMS_zz4l_massErrZX_ld_frac_{0:.0f}".format(self.channel)
	rfv_EBE_zx_frac = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_ld_frac']+")", ROOT.RooArgList(self.MH)) 
	pdfErrZX = ROOT.RooAddPdf("pdfErrZX","pdfErrZX", EBE_zx_ld, EBE_zx_ln, rfv_EBE_zx_frac)


	sig_ggHErr = ROOT.RooProdPdf("sig_ggHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(signalCB_ggH), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrS), ROOT.RooArgSet(RelErr)))
	sig_VBFErr = ROOT.RooProdPdf("sig_VBFErr","BW (X) CB * pdfErr", ROOT.RooArgSet(signalCB_VBF), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrS), ROOT.RooArgSet(RelErr)))
	sig_WHErr = ROOT.RooProdPdf("sig_WHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(signalCB_WH), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrS), ROOT.RooArgSet(RelErr)))
	sig_ZHErr = ROOT.RooProdPdf("sig_ZHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(signalCB_ZH), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrS), ROOT.RooArgSet(RelErr)))
	sig_ttHErr = ROOT.RooProdPdf("sig_ttHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(signalCB_ttH), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrS), ROOT.RooArgSet(RelErr)))
	#------------------------------------------------end bIncludingError

        #------------------------------------------------begin bUseVBF
        if(self.bVBF):

            VD_name = "CMS_zz4l_Fisher"
            VD = ROOT.RooRealVar(VD_name,VD_name,0,2.)
            VD.setBins(50)
            pt_name = "CMS_zz4l_Pt"
            pt = ROOT.RooRealVar(pt_name,pt_name,0,200.)
            pt.setBins(50)

            sig2d_ggH_Fisher = ROOT.RooProdPdf()
            sig2d_qqH_Fisher = ROOT.RooProdPdf()
            sig2d_WH_Fisher = ROOT.RooProdPdf()
            sig2d_ZH_Fisher = ROOT.RooProdPdf()
            sig2d_ttH_Fisher = ROOT.RooProdPdf()
            sig2d_ggH_Pt = ROOT.RooProdPdf()
            sig2d_qqH_Pt = ROOT.RooProdPdf()
            sig2d_WH_Pt = ROOT.RooProdPdf()
            sig2d_ZH_Pt = ROOT.RooProdPdf()
            sig2d_ttH_Pt = ROOT.RooProdPdf()
            
            sigCB2d_ggH_Fisher = ROOT.RooProdPdf()
            sigCB2d_qqH_Fisher = ROOT.RooProdPdf()
            sigCB2d_WH_Fisher = ROOT.RooProdPdf()
            sigCB2d_ZH_Fisher = ROOT.RooProdPdf()
            sigCB2d_ttH_Fisher = ROOT.RooProdPdf()
            sigCB2d_ggH_Pt = ROOT.RooProdPdf()
            sigCB2d_qqH_Pt = ROOT.RooProdPdf()
            sigCB2d_WH_Pt = ROOT.RooProdPdf()
            sigCB2d_ZH_Pt = ROOT.RooProdPdf()
            sigCB2d_ttH_Pt = ROOT.RooProdPdf()
        if(self.bVBF and self.VBFcat):
            
            ggHtempFileName = "{0}/ggH_fisher.root".format(self.templateDir)
            ggHtempFile = ROOT.TFile(ggHtempFileName)
            ggHtemplate = ggHtempFile.Get("h_Fisher")
            ggHtemplate_Up = ggHtempFile.Get("h_Fisher_up")
            ggHtemplate_Dn = ggHtempFile.Get("h_Fisher_dn")

            Fisher_ggH_dataHist = ROOT.RooDataHist("temp_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ggHtemplate)
            Fisher_ggH_dataHist_Up = ROOT.RooDataHist("temp_ggH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ggH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ggHtemplate_Up)
            Fisher_ggH_dataHist_Dn = ROOT.RooDataHist("temp_ggH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ggH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ggHtemplate_Dn)
            
            TemplateName = "FisherTempDataHist_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ggH_dataHist)
            TemplateName = "FisherTempDataHist_ggH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ggH_dataHist_Up)
            TemplateName = "FisherTempDataHist_ggH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ggH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ggH_dataHist_Dn)

            qqHtempFileName = "{0}/qqH_fisher.root".format(self.templateDir)
            qqHtempFile = ROOT.TFile(qqHtempFileName)
            qqHtemplate = qqHtempFile.Get("h_Fisher")
            qqHtemplate_Up = qqHtempFile.Get("h_Fisher_up")
            qqHtemplate_Dn = qqHtempFile.Get("h_Fisher_dn")
            
            Fisher_qqH_dataHist = ROOT.RooDataHist("temp_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),qqHtemplate)
            Fisher_qqH_dataHist_Up = ROOT.RooDataHist("temp_qqH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_qqH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),qqHtemplate_Up)
            Fisher_qqH_dataHist_Dn = ROOT.RooDataHist("temp_qqH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_qqH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),qqHtemplate_Dn)
            
            TemplateName = "FisherTempDataHist_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_qqH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_qqH_dataHist)
            TemplateName = "FisherTempDataHist_qqH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_qqH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_qqH_dataHist_Up)
            TemplateName = "FisherTempDataHist_qqH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_qqH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_qqH_dataHist_Dn)

            WHtempFileName = "{0}/WH_fisher.root".format(self.templateDir)
            WHtempFile = ROOT.TFile(WHtempFileName)
            WHtemplate = WHtempFile.Get("h_Fisher")
            
            
            Fisher_WH_dataHist = ROOT.RooDataHist("temp_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),WHtemplate)
            
            TemplateName = "FisherTempDataHist_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_WH_dataHist)

            ZHtempFileName = "{0}/ZH_fisher.root".format(self.templateDir)
            ZHtempFile = ROOT.TFile(ZHtempFileName)
            ZHtemplate = ZHtempFile.Get("h_Fisher")
            
            Fisher_ZH_dataHist = ROOT.RooDataHist("temp_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ZHtemplate)
            
            TemplateName = "FisherTempDataHist_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ZH_dataHist)

            ttHtempFileName = "{0}/ttH_fisher.root".format(self.templateDir)
            ttHtempFile = ROOT.TFile(ttHtempFileName)
            ttHtemplate = ttHtempFile.Get("h_Fisher")
            
            Fisher_ttH_dataHist = ROOT.RooDataHist("temp_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ttHtemplate)
            
            TemplateName = "FisherTempDataHist_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ttH_dataHist)

            FisherList_ggH = ROOT.RooArgList()  
            FisherList_qqH = ROOT.RooArgList()
            FisherList_WH  = ROOT.RooArgList()
            FisherList_ZH  = ROOT.RooArgList()
            FisherList_ttH = ROOT.RooArgList()

            if(self.FisherMorph):
                
                FisherList_ggH.add(FisherTemplatePdf_ggH)
                FisherList_ggH.add(FisherTemplatePdf_ggH_Up)
                FisherList_ggH.add(FisherTemplatePdf_ggH_Dn)  
                
                FisherList_qqH.add(FisherTemplatePdf_qqH)
                FisherList_qqH.add(FisherTemplatePdf_qqH_Up)
                FisherList_qqH.add(FisherTemplatePdf_qqH_Dn)
                
                FisherList_WH.add(FisherTemplatePdf_WH) 
                
                FisherList_ZH.add(FisherTemplatePdf_ZH) 
                
                FisherList_ttH.add(FisherTemplatePdf_ttH)
            else:
            
                FisherList_ggH.add(FisherTemplatePdf_ggH)
                FisherList_qqH.add(FisherTemplatePdf_qqH)
                FisherList_WH.add(FisherTemplatePdf_WH)
                FisherList_ZH.add(FisherTemplatePdf_ZH)
                FisherList_ttH.add(FisherTemplatePdf_ttH)

            morphFisherVarName = "CMS_zz4l_ggH_Fisher_sys"
            alphaMorphFisher_ggH = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            morphFisherVarName = "CMS_zz4l_qqH_Fisher_sys"
            alphaMorphFisher_qqH = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            morphFisherVarName = "CMS_zz4l_WH_Fisher_sys"
            alphaMorphFisher_WH = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            morphFisherVarName = "CMS_zz4l_ZH_Fisher_sys"
            alphaMorphFisher_ZH = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            morphFisherVarName = "CMS_zz4l_ttH_Fisher_sys"
            alphaMorphFisher_ttH = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            
            if(self.FisherMorph):
                alphaMorphFisher_ggH.setConstant(False)
                alphaMorphFisher_qqH.setConstant(False)
            else:
                alphaMorphFisher_ggH.setConstant(True)
                alphaMorphFisher_qqH.setConstant(True)
            alphaMorphFisher_WH.setConstant(True)
            alphaMorphFisher_ZH.setConstant(True)
            alphaMorphFisher_ttH.setConstant(True)
        
            morphVarListFisher_ggH = ROOT.RooArgList()
            morphVarListFisher_qqH = ROOT.RooArgList()
            morphVarListFisher_WH = ROOT.RooArgList()
            morphVarListFisher_ZH = ROOT.RooArgList()
            morphVarListFisher_ttH = ROOT.RooArgList()
                
            if(self.FisherMorph):
                morphVarListFisher_ggH.add(alphaMorphFisher_ggH)
                morphVarListFisher_qqH.add(alphaMorphFisher_qqH)
        
            TemplateName = "FisherTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_ggH,morphVarListFisher_ggH,1.0,1)
        
            TemplateName = "FisherTemplateMorphPdf_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_qqH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_qqH,morphVarListFisher_qqH,1.0,1)
        
            TemplateName = "FisherTemplateMorphPdf_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_WH,morphVarListFisher_WH,1.0,1)
        
            TemplateName = "FisherTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_ZH,morphVarListFisher_ZH,1.0,1)
            
            TemplateName = "FisherTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_ttH,morphVarListFisher_ttH,1.0,1)

            sig2d_ggH_Fisher = ROOT.RooProdPdf("sig2d_ggH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_ggH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ggH_HM,sig_ggH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ggH),ROOT.RooArgSet(VD)))
            sig2d_qqH_Fisher = ROOT.RooProdPdf("sig2d_qqH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_qqH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_VBF_HM,sig_VBF,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_qqH),ROOT.RooArgSet(VD)))
            sig2d_WH_Fisher = ROOT.RooProdPdf("sig2d_WH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_WH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_WH_HM,sig_WH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_WH),ROOT.RooArgSet(VD)))
            sig2d_ZH_Fisher = ROOT.RooProdPdf("sig2d_ZH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_ZH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ZH_HM,sig_ZH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ZH),ROOT.RooArgSet(VD)))
            sig2d_ttH_Fisher = ROOT.RooProdPdf("sig2d_ttH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_ttH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ttH_HM,sig_ttH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ttH),ROOT.RooArgSet(VD)))
        
            sigCB2d_ggH_Fisher = ROOT.RooProdPdf("sigCB2d_ggH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_ggH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ggHErr,signalCB_ggH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ggH),ROOT.RooArgSet(VD)))
            sigCB2d_qqH_Fisher = ROOT.RooProdPdf("sigCB2d_qqH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_qqH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_VBFErr,signalCB_VBF,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_qqH),ROOT.RooArgSet(VD)))
            sigCB2d_WH_Fisher = ROOT.RooProdPdf("sigCB2d_WH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_ZH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_WHErr,signalCB_WH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_WH),ROOT.RooArgSet(VD)))
            sigCB2d_ZH_Fisher = ROOT.RooProdPdf("sigCB2d_ZH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_WH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ZHErr,signalCB_ZH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ZH),ROOT.RooArgSet(VD)))
            sigCB2d_ttH_Fisher = ROOT.RooProdPdf("sigCB2d_ttH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_ttH_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ttHErr,signalCB_ttH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ttH),ROOT.RooArgSet(VD)))

        if(self.bVBF and not self.VBFcat):
            
            ggHtempPtFileName = "{0}/Pt_mZZ_gg_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            ggHtempPtFile = ROOT.TFile(ggHtempPtFileName)
            ggHtemplatePt = ggHtempPtFile.Get("h_Ptmzz_mzz")
            ggHtemplatePt_Up = ggHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            ggHtemplatePt_Dn = ggHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_ggH_dataHist = ROOT.RooDataHist("tempPt_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ggHtemplatePt)
            Pt_ggH_dataHist_Up = ROOT.RooDataHist("tempPt_ggH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ggH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ggHtemplatePt_Up)
            Pt_ggH_dataHist_Dn = ROOT.RooDataHist("tempPt_ggH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ggH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ggHtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ggH_dataHist)
            TemplateName = "PtTempDataHist_ggH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ggH_dataHist_Up)
            TemplateName = "PtTempDataHist_ggH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ggH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ggH_dataHist_Dn)

            qqHtempPtFileName = "{0}/Pt_mZZ_vbf_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            qqHtempPtFile = ROOT.TFile(qqHtempPtFileName)
            qqHtemplatePt = qqHtempPtFile.Get("h_Ptmzz_mzz")
            qqHtemplatePt_Up = qqHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            qqHtemplatePt_Dn = qqHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_qqH_dataHist = ROOT.RooDataHist("tempPt_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),qqHtemplatePt)
            Pt_qqH_dataHist_Up = ROOT.RooDataHist("tempPt_qqH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_qqH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),qqHtemplatePt_Up)
            Pt_qqH_dataHist_Dn = ROOT.RooDataHist("tempPt_qqH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_qqH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),qqHtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_qqH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_qqH_dataHist)
            TemplateName = "PtTempDataHist_qqH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_qqH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_qqH_dataHist_Up)
            TemplateName = "PtTempDataHist_qqH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_qqH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_qqH_dataHist_Dn)

            ttHtempPtFileName = "{0}/Pt_mZZ_tth_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            ttHtempPtFile = ROOT.TFile(ttHtempPtFileName)
            ttHtemplatePt = ttHtempPtFile.Get("h_Ptmzz_mzz")
            ttHtemplatePt_Up = ttHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            ttHtemplatePt_Dn = ttHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_ttH_dataHist = ROOT.RooDataHist("tempPt_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ttHtemplatePt)
            Pt_ttH_dataHist_Up = ROOT.RooDataHist("tempPt_ttH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ttH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ttHtemplatePt_Up)
            Pt_ttH_dataHist_Dn = ROOT.RooDataHist("tempPt_ttH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ttH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ttHtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ttH_dataHist)
            TemplateName = "PtTempDataHist_ttH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ttH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ttH_dataHist_Up)
            TemplateName = "PtTempDataHist_ttH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ttH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ttH_dataHist_Dn)

            WHtempPtFileName = "{0}/Pt_mZZ_wh_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            WHtempPtFile = ROOT.TFile(WHtempPtFileName)
            WHtemplatePt = WHtempPtFile.Get("h_Ptmzz_mzz")
            WHtemplatePt_Up = WHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            WHtemplatePt_Dn = WHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_WH_dataHist = ROOT.RooDataHist("tempPt_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),WHtemplatePt)
            Pt_WH_dataHist_Up = ROOT.RooDataHist("tempPt_WH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_WH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),WHtemplatePt_Up)
            Pt_WH_dataHist_Dn = ROOT.RooDataHist("tempPt_WH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_WH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),WHtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_WH_dataHist)
            TemplateName = "PtTempDataHist_WH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_WH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_WH_dataHist_Up)
            TemplateName = "PtTempDataHist_WH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_WH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_WH_dataHist_Dn)

            ZHtempPtFileName = "{0}/Pt_mZZ_zh_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            ZHtempPtFile = ROOT.TFile(ZHtempPtFileName)
            ZHtemplatePt = ZHtempPtFile.Get("h_Ptmzz_mzz")
            ZHtemplatePt_Up = ZHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            ZHtemplatePt_Dn = ZHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_ZH_dataHist = ROOT.RooDataHist("tempPt_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ZHtemplatePt)
            Pt_ZH_dataHist_Up = ROOT.RooDataHist("tempPt_ZH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ZH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ZHtemplatePt_Up)
            Pt_ZH_dataHist_Dn = ROOT.RooDataHist("tempPt_ZH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ZH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ZHtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ZH_dataHist)
            TemplateName = "PtTempDataHist_ZH_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ZH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ZH_dataHist_Up)
            TemplateName = "PtTempDataHist_ZH_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ZH_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ZH_dataHist_Dn)

            PtList_ggH = ROOT.RooArgList()  
            PtList_qqH = ROOT.RooArgList()
            PtList_WH  = ROOT.RooArgList()
            PtList_ZH  = ROOT.RooArgList()
            PtList_ttH = ROOT.RooArgList()

            if(self.PtMorph):
                PtList_ggH.add(PtTemplatePdf_ggH)
                PtList_ggH.add(PtTemplatePdf_ggH_Up)
                PtList_ggH.add(PtTemplatePdf_ggH_Dn)  
                
                PtList_qqH.add(PtTemplatePdf_qqH)
                PtList_qqH.add(PtTemplatePdf_qqH_Up)
                PtList_qqH.add(PtTemplatePdf_qqH_Dn)  
                
                PtList_WH.add(PtTemplatePdf_WH)
                PtList_WH.add(PtTemplatePdf_WH_Up)
                PtList_WH.add(PtTemplatePdf_WH_Dn)  
                
                PtList_ZH.add(PtTemplatePdf_ZH)
                PtList_ZH.add(PtTemplatePdf_ZH_Up)
                PtList_ZH.add(PtTemplatePdf_ZH_Dn)  
                
                PtList_ttH.add(PtTemplatePdf_ttH)
                PtList_ttH.add(PtTemplatePdf_ttH_Up)
                PtList_ttH.add(PtTemplatePdf_ttH_Dn)
            else:
            
                PtList_ggH.add(PtTemplatePdf_ggH)
                PtList_qqH.add(PtTemplatePdf_qqH)
                PtList_WH.add(PtTemplatePdf_WH)
                PtList_ZH.add(PtTemplatePdf_ZH)
                PtList_ttH.add(PtTemplatePdf_ttH)

            morphPtVarName = "CMS_zz4l_ggH_Pt_sys"
            alphaMorphPt_ggH = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            morphPtVarName = "CMS_zz4l_qqH_Pt_sys"
            alphaMorphPt_qqH = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            morphPtVarName = "CMS_zz4l_VH_Pt_sys"
            alphaMorphPt_VH = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            morphPtVarName = "CMS_zz4l_ttH_Pt_sys"
            alphaMorphPt_ttH = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            if(self.PtMorph):
                alphaMorphPt_ggH.setConstant(False)
                alphaMorphPt_qqH.setConstant(False)
                alphaMorphPt_VH.setConstant(False)
                alphaMorphPt_ttH.setConstant(False)
            else:
                alphaMorphPt_ggH.setConstant(True)
                alphaMorphPt_qqH.setConstant(True)
                alphaMorphPt_VH.setConstant(True)
                alphaMorphPt_ttH.setConstant(True)
        
            morphVarListPt_ggH = ROOT.RooArgList()
            morphVarListPt_qqH = ROOT.RooArgList()
            morphVarListPt_VH = ROOT.RooArgList()
            morphVarListPt_ttH = ROOT.RooArgList()
                
            if(self.PtMorph):
                morphVarListPt_ggH.add(alphaMorphPt_ggH)
                morphVarListPt_qqH.add(alphaMorphPt_qqH)
                morphVarListPt_VH.add(alphaMorphPt_VH)
                morphVarListPt_ttH.add(alphaMorphPt_ttH)
        
            PtTemplateName = "PtTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_ggH,morphVarListPt_ggH,1.0,1)
        
            PtTemplateName = "PtTemplateMorphPdf_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_qqH = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_qqH,morphVarListPt_qqH,1.0,1)
        
            PtTemplateName = "PtTemplateMorphPdf_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_WH,morphVarListPt_VH,1.0,1)
        
            PtTemplateName = "PtTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_ZH,morphVarListPt_VH,1.0,1)
            
            PtTemplateName = "PtTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_ttH,morphVarListPt_ttH,1.0,1)

            sig2d_ggH_Pt = ROOT.RooProdPdf("sig2d_ggH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_ggH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ggH_HM,sig_ggH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ggH),ROOT.RooArgSet(pt)))
            sig2d_qqH_Pt = ROOT.RooProdPdf("sig2d_VBF_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_VBF_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_VBF_HM,sig_VBF,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_qqH),ROOT.RooArgSet(pt)))
            sig2d_WH_Pt = ROOT.RooProdPdf("sig2d_WH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_WH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_WH_HM,sig_WH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_WH),ROOT.RooArgSet(pt)))
            sig2d_ZH_Pt = ROOT.RooProdPdf("sig2d_ZH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_ZH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ZH_HM,sig_ZH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ZH),ROOT.RooArgSet(pt)))
            sig2d_ttH_Pt = ROOT.RooProdPdf("sig2d_ttH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sig2d_ttH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ttH_HM,sig_ttH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ttH),ROOT.RooArgSet(pt)))

            sigCB2d_ggH_Pt = ROOT.RooProdPdf("sigCB2d_ggH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_ggH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ggHErr,signalCB_ggH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ggH),ROOT.RooArgSet(pt)))
            sigCB2d_qqH_Pt = ROOT.RooProdPdf("sigCB2d_qqH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_qqH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_VBFErr,signalCB_VBF,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_qqH),ROOT.RooArgSet(pt)))
            sigCB2d_WH_Pt = ROOT.RooProdPdf("sigCB2d_WH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_ZH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_WHErr,signalCB_WH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_WH),ROOT.RooArgSet(pt)))
            sigCB2d_ZH_Pt = ROOT.RooProdPdf("sigCB2d_ZH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_WH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ZHErr,signalCB_ZH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ZH),ROOT.RooArgSet(pt)))
            sigCB2d_ttH_Pt = ROOT.RooProdPdf("sigCB2d_ttH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"sigCB2d_ttH_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(sig_ttHErr,signalCB_ttH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ttH),ROOT.RooArgSet(pt)))


        ## --------------------------- MELA 2D PDFS ------------------------- ##
        discVarName = ""
        if self.useMEKDTemplates:
            discVarName = "mekdLD"
        else:
            discVarName = "melaLD"
    
        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        
        sigTempFile = ROOT.TFile(templateSigName)
        sigTemplate = sigTempFile.Get("h_mzzD")
        sigTemplate_Up = sigTempFile.Get("h_mzzD_up")
        sigTemplate_Down = sigTempFile.Get("h_mzzD_dn")

        dBins = sigTemplate.GetYaxis().GetNbins()
        dLow = sigTemplate.GetYaxis().GetXmin()
        dHigh = sigTemplate.GetYaxis().GetXmax()
        D = ROOT.RooRealVar(discVarName,discVarName,dLow,dHigh)
        D.setBins(dBins)
        print "discVarName ", discVarName, " dLow ", dLow, " dHigh ", dHigh, " dBins ", dBins
        
        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate)
        TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Up)
        TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Down)

        
        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ggH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ggH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_VBF = ROOT.RooHistPdf(TemplateName,TemplateName,RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_VBF_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_VBF_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_VBF_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_VBF_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_WH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_WH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_WH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_WH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ZH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ZH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ttH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ttH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)

        funcList_ggH = ROOT.RooArgList()  
        funcList_VBF = ROOT.RooArgList()
        funcList_WH  = ROOT.RooArgList()
        funcList_ZH  = ROOT.RooArgList()
        funcList_ttH = ROOT.RooArgList()


        if(self.isAltSig):
            #only ggH because if we do hypothesis testing we sum up over the channels in any case
              templateSigName = "{0}/Dsignal{2}_{1}.root".format(self.templateDir,self.appendName, self.appendHypType)
              print 'Taking 2D template for ALT signal from ',templateSigName
              sigTempFile = ROOT.TFile(templateSigName)
              sigTemplate = sigTempFile.Get("h_mzzD")
              sigTemplate_Up = sigTempFile.Get("h_mzzD_up")
              sigTemplate_Down = sigTempFile.Get("h_mzzD_dn")
              
              TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTempDataHist_ALT = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate)
              TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTempDataHist_ALT_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Up)
              TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTempDataHist_ALT_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Down)
              TemplateName = "sigTemplatePdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
              sigTemplatePdf_ggH_ALT = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_ALT)
              TemplateName = "sigTemplatePdf_ggH{2}_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
              sigTemplatePdf_ggH_ALT_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_ALT_Up)
              TemplateName = "sigTemplatePdf_ggH{2}_Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTemplatePdf_ggH_ALT_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_ALT_Down)

        funcList_ggH_ALT = ROOT.RooArgList() 

        if(self.sigMorph):
            
            funcList_ggH.add(sigTemplatePdf_ggH)
            funcList_ggH.add(sigTemplatePdf_ggH_Up)
            funcList_ggH.add(sigTemplatePdf_ggH_Down)  
            
            funcList_VBF.add(sigTemplatePdf_VBF)
            funcList_VBF.add(sigTemplatePdf_VBF_Up)
            funcList_VBF.add(sigTemplatePdf_VBF_Down)  
            
            funcList_WH.add(sigTemplatePdf_WH)
            funcList_WH.add(sigTemplatePdf_WH_Up)
            funcList_WH.add(sigTemplatePdf_WH_Down)  
            
            funcList_ZH.add(sigTemplatePdf_ZH)
            funcList_ZH.add(sigTemplatePdf_ZH_Up)
            funcList_ZH.add(sigTemplatePdf_ZH_Down)  
            
            funcList_ttH.add(sigTemplatePdf_ttH)
            funcList_ttH.add(sigTemplatePdf_ttH_Up)
            funcList_ttH.add(sigTemplatePdf_ttH_Down)  
            if(self.isAltSig):
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT)
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT_Up)
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT_Down)
        else:
            
            funcList_ggH.add(sigTemplatePdf_ggH)
            funcList_VBF.add(sigTemplatePdf_VBF)
            funcList_WH.add(sigTemplatePdf_WH)
            funcList_ZH.add(sigTemplatePdf_ZH)
            funcList_ttH.add(sigTemplatePdf_ttH)
            if(self.isAltSig):
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT)
    
        morphSigVarName = "CMS_zz4l_sigMELA_{0:.0f}".format(self.channel)
        alphaMorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if(self.sigMorph): alphaMorphSig.setConstant(False)
        else: alphaMorphSig.setConstant(True)
        
        morphVarListSig = ROOT.RooArgList()
    
        if(self.sigMorph): morphVarListSig.add(alphaMorphSig)  ## just one morphing for all signal processes (fully correlated systematics)
        
        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ggH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_VBF = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_VBF,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_WH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ZH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ttH,morphVarListSig,1.0,1)
        if(self.isAltSig):
            TemplateName = "sigTemplateMorphPdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.appendHypType)
            sigTemplateMorphPdf_ggH_ALT = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ggH_ALT,morphVarListSig,1.0,1)
    

	####  ----------------------- mekd  parametrized double gaussian stuffs  -------------------------
	discVarName = "mekd"
	MEKD = ROOT.RooRealVar(discVarName, discVarName, -5, 15);
	if theMEKD: 
		name = "mekd_sig_a0_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_sig_a0 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_sig_a0_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_sig_a1_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_sig_a1 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_sig_a1_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_sig_a2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_sig_a2 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_sig_a2_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_sig_a3_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_sig_a3 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_sig_a3_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_sig_a4_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_sig_a4 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_sig_a4_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		sigTemplateMorphPdf_ggH = ROOT.RooGenericPdf("mekd_sig_ggH", "mekd_sig_ggH", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_sig_a0, mekd_sig_a1, mekd_sig_a2, mekd_sig_a3, mekd_sig_a4))
		sigTemplateMorphPdf_VBF = sigTemplateMorphPdf_ggH 
		sigTemplateMorphPdf_WH = sigTemplateMorphPdf_ggH 
		sigTemplateMorphPdf_ZH = sigTemplateMorphPdf_ggH 
		sigTemplateMorphPdf_ttH = sigTemplateMorphPdf_ggH 
		print "\n \n mekd_sig_a2 channel ",self.channel
		m = 100
		while m >= 100 and m < 150:
			CMS_zz4l_mass.setVal(m)
			m = m + 0.1
			if mekd_sig_a2.getVal() < 0 : print m, mekd_sig_a2.getVal() 
			if mekd_sig_a2.getVal() > 1 : print m, mekd_sig_a2.getVal() 
		print "\n \n mekd_sig_a1 channel ",self.channel
		m = 100
		while m >= 100 and m < 150:
			CMS_zz4l_mass.setVal(m)
			m = m + 0.1
			if mekd_sig_a1.getVal() <= 0 : print m, mekd_sig_a1.getVal() 
		print "\n \n mekd_sig_a4 channel ",self.channel
		m = 100
		while m >= 100 and m < 150:
			CMS_zz4l_mass.setVal(m)
			m = m + 0.1
			if mekd_sig_a4.getVal() <= 0 : print m, mekd_sig_a4.getVal() 
	####  ----------------------- end mekd -----------------------------------------------------------
        sig2d_ggH = ROOT.RooProdPdf("sig2d_ggH","sig2d_ggH",ROOT.RooArgSet(self.getVariable(sig_ggH_HM,sig_ggH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sig2d_VBF = ROOT.RooProdPdf("sig2d_VBF","sig2d_VBF",ROOT.RooArgSet(self.getVariable(sig_VBF_HM,sig_VBF,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sig2d_WH = ROOT.RooProdPdf("sig2d_WH","sig2d_WH",ROOT.RooArgSet(self.getVariable(sig_WH_HM,sig_WH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sig2d_ZH = ROOT.RooProdPdf("sig2d_ZH","sig2d_ZH",ROOT.RooArgSet(self.getVariable(sig_ZH_HM,sig_ZH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sig2d_ttH = ROOT.RooProdPdf("sig2d_ttH","sig2d_ttH",ROOT.RooArgSet(self.getVariable(sig_ttH_HM,sig_ttH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                
        sigCB2d_ggH = ROOT.RooProdPdf("sigCB2d_ggH","sigCB2d_ggH",ROOT.RooArgSet(self.getVariable(sig_ggHErr,signalCB_ggH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sigCB2d_VBF = ROOT.RooProdPdf("sigCB2d_VBF","sigCB2d_VBF",ROOT.RooArgSet(self.getVariable(sig_VBFErr,signalCB_VBF,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sigCB2d_WH = ROOT.RooProdPdf("sigCB2d_WH","sigCB2d_WH",ROOT.RooArgSet(self.getVariable(sig_WHErr,signalCB_WH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sigCB2d_ZH = ROOT.RooProdPdf("sigCB2d_ZH","sigCB2d_ZH",ROOT.RooArgSet(self.getVariable(sig_ZHErr,signalCB_ZH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        sigCB2d_ttH = ROOT.RooProdPdf("sigCB2d_ttH","sigCB2d_ttH",ROOT.RooArgSet(self.getVariable(sig_ttHErr,signalCB_ttH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        
        if(self.isAltSig):
            sig2d_ggH_ALT = ROOT.RooProdPdf("sig2d_ggH{0}".format(self.appendHypType),"sig2d_ggH{0}".format(self.appendHypType),ROOT.RooArgSet(sig_ggH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH_ALT),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
            sigCB2d_ggH_ALT = ROOT.RooProdPdf("sigCB2d_ggH{0}".format(self.appendHypType),"sigCB2d_ggH{0}".format(self.appendHypType),ROOT.RooArgSet(signalCB_ggH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH_ALT),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))

        if(self.bVBF):
            if (self.Use3D):
                sig2d_ggH_VBF_KD = ROOT.RooProdPdf("sig2d_ggH_VBF_KD","sig2d_ggH_VBF_KD",ROOT.RooArgSet(self.getVariable(sig2d_ggH_Fisher,sig2d_ggH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sig2d_qqH_VBF_KD = ROOT.RooProdPdf("sig2d_qqH_VBF_KD","sig2d_qqH_VBF_KD",ROOT.RooArgSet(self.getVariable(sig2d_qqH_Fisher,sig2d_qqH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sig2d_WH_VBF_KD = ROOT.RooProdPdf("sig2d_WH_VBF_KD","sig2d_WH_VBF_KD",ROOT.RooArgSet(self.getVariable(sig2d_WH_Fisher,sig2d_WH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sig2d_ZH_VBF_KD = ROOT.RooProdPdf("sig2d_ZH_VBF_KD","sig2d_ZH_VBF_KD",ROOT.RooArgSet(self.getVariable(sig2d_ZH_Fisher,sig2d_ZH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sig2d_ttH_VBF_KD = ROOT.RooProdPdf("sig2d_ttH_VBF_KD","sig2d_ttH_VBF_KD",ROOT.RooArgSet(self.getVariable(sig2d_ttH_Fisher,sig2d_ttH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
            
                sigCB2d_ggH_VBF_KD = ROOT.RooProdPdf("sigCB2d_ggH_VBF_KD","sigCB2d_ggH_VBF_KD",ROOT.RooArgSet(self.getVariable(sigCB2d_ggH_Fisher,sigCB2d_ggH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sigCB2d_qqH_VBF_KD = ROOT.RooProdPdf("sigCB2d_qqH_VBF_KD","sigCB2d_qqH_VBF_KD",ROOT.RooArgSet(self.getVariable(sigCB2d_qqH_Fisher,sigCB2d_qqH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sigCB2d_WH_VBF_KD = ROOT.RooProdPdf("sigCB2d_WH_VBF_KD","sigCB2d_WH_VBF_KD",ROOT.RooArgSet(self.getVariable(sigCB2d_WH_Fisher,sigCB2d_WH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sigCB2d_ZH_VBF_KD = ROOT.RooProdPdf("sigCB2d_ZH_VBF_KD","sigCB2d_ZH_VBF_KD",ROOT.RooArgSet(self.getVariable(sigCB2d_ZH_Fisher,sigCB2d_ZH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
                sigCB2d_ttH_VBF_KD = ROOT.RooProdPdf("sigCB2d_ttH_VBF_KD","sigCB2d_ttH_VBF_KD",ROOT.RooArgSet(self.getVariable(sigCB2d_ttH_Fisher,sigCB2d_ttH_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))

        ## --------------------------- superMELA 1D PDFS ------------------------- ##

        superDiscVarName = "supermelaLD"
        SD = ROOT.RooRealVar(superDiscVarName,superDiscVarName,0,1)
    
        templateSDSigName = "{0}/Dsignal_superMELA_{1}.root".format(self.templateDir ,self.appendName)
        sigTempSDFile = ROOT.TFile(templateSDSigName)
        sigTemplateSD = sigTempSDFile.Get("hSuperD_sig") 
        
        TemplateSDName = "sigTempSDDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempSDDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD),sigTemplateSD)
        
        TemplateSDName = "sigTemplateSDPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ggH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_VBF = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_WH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ZH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ttH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        print sigTemplateSDPdf_ttH

        ##--------------##

        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution
        if not self.bVBF:
            name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
            name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
            name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
            name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
            name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
            name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
            name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
            name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
            name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
            name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
            CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
            name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
            CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
            name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
            CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
            name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
            CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
            name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
            CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        else:
            name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
            name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
            name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
            name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
            name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
            name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
            name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
            name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
            name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
            name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}_{2}".format( self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
            name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
            name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
            name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
            name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat )
            CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        

        if (DEBUG) :
            print "qqZZshape_a0 = ",theInputs['qqZZshape_a0']
            print "qqZZshape_a1 = ",theInputs['qqZZshape_a1']
            print "qqZZshape_a2 = ",theInputs['qqZZshape_a2']
            print "qqZZshape_a3 = ",theInputs['qqZZshape_a3']
            print "qqZZshape_a4 = ",theInputs['qqZZshape_a4']
            print "qqZZshape_a5 = ",theInputs['qqZZshape_a5']
            print "qqZZshape_a6 = ",theInputs['qqZZshape_a6']
            print "qqZZshape_a7 = ",theInputs['qqZZshape_a7']
            print "qqZZshape_a8 = ",theInputs['qqZZshape_a8']
            print "qqZZshape_a9 = ",theInputs['qqZZshape_a9']
            print "qqZZshape_a10 = ",theInputs['qqZZshape_a10']
            print "qqZZshape_a11 = ",theInputs['qqZZshape_a11']
            print "qqZZshape_a12 = ",theInputs['qqZZshape_a12']
            print "qqZZshape_a13 = ",theInputs['qqZZshape_a13']

        
        CMS_qqzzbkg_a0.setVal(theInputs['qqZZshape_a0'])
        CMS_qqzzbkg_a1.setVal(theInputs['qqZZshape_a1'])
        CMS_qqzzbkg_a2.setVal(theInputs['qqZZshape_a2'])
        CMS_qqzzbkg_a3.setVal(theInputs['qqZZshape_a3'])
        CMS_qqzzbkg_a4.setVal(theInputs['qqZZshape_a4'])
        CMS_qqzzbkg_a5.setVal(theInputs['qqZZshape_a5'])
        CMS_qqzzbkg_a6.setVal(theInputs['qqZZshape_a6'])
        CMS_qqzzbkg_a7.setVal(theInputs['qqZZshape_a7'])
        CMS_qqzzbkg_a8.setVal(theInputs['qqZZshape_a8'])
        CMS_qqzzbkg_a9.setVal(theInputs['qqZZshape_a9'])
        CMS_qqzzbkg_a10.setVal(theInputs['qqZZshape_a10'])
        CMS_qqzzbkg_a11.setVal(theInputs['qqZZshape_a11'])
        CMS_qqzzbkg_a12.setVal(theInputs['qqZZshape_a12'])
        CMS_qqzzbkg_a13.setVal(theInputs['qqZZshape_a13'])
        
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
        
        bkg_qqzz = ROOT.RooqqZZPdf_v2("bkg_qqzzTmp","bkg_qqzzTmp",CMS_zz4l_mass,CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13)
        
        ## ggZZ contribution
        if not self.bVBF:
            name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
            name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
            name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
            name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
            name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
            CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
            name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
            name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
            name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
            name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
            CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
            name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
            CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        else:
            name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
            name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
            name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
            name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
            name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat )
            CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
            name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
            name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
            name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
            name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat ) 
            CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
            name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}_{2}".format( self.channel, self.sqrts,self.VBFcat )
            CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        
        CMS_ggzzbkg_a0.setVal(theInputs['ggZZshape_a0'])
        CMS_ggzzbkg_a1.setVal(theInputs['ggZZshape_a1'])
        CMS_ggzzbkg_a2.setVal(theInputs['ggZZshape_a2'])
        CMS_ggzzbkg_a3.setVal(theInputs['ggZZshape_a3'])
        CMS_ggzzbkg_a4.setVal(theInputs['ggZZshape_a4'])
        CMS_ggzzbkg_a5.setVal(theInputs['ggZZshape_a5'])
        CMS_ggzzbkg_a6.setVal(theInputs['ggZZshape_a6'])
        CMS_ggzzbkg_a7.setVal(theInputs['ggZZshape_a7'])
        CMS_ggzzbkg_a8.setVal(theInputs['ggZZshape_a8'])
        CMS_ggzzbkg_a9.setVal(theInputs['ggZZshape_a9'])
        
        CMS_ggzzbkg_a0.setConstant(True)
        CMS_ggzzbkg_a1.setConstant(True)
        CMS_ggzzbkg_a2.setConstant(True)
        CMS_ggzzbkg_a3.setConstant(True)
        CMS_ggzzbkg_a4.setConstant(True)
        CMS_ggzzbkg_a5.setConstant(True)
        CMS_ggzzbkg_a6.setConstant(True)
        CMS_ggzzbkg_a7.setConstant(True)
        CMS_ggzzbkg_a8.setConstant(True)
        CMS_ggzzbkg_a9.setConstant(True)

        if (DEBUG) :
            print "ggZZshape_a0 = ",theInputs['ggZZshape_a0']
            print "ggZZshape_a1 = ",theInputs['ggZZshape_a1']
            print "ggZZshape_a2 = ",theInputs['ggZZshape_a2']
            print "ggZZshape_a3 = ",theInputs['ggZZshape_a3']
            print "ggZZshape_a4 = ",theInputs['ggZZshape_a4']
            print "ggZZshape_a5 = ",theInputs['ggZZshape_a5']
            print "ggZZshape_a6 = ",theInputs['ggZZshape_a6']
            print "ggZZshape_a7 = ",theInputs['ggZZshape_a7']
            print "ggZZshape_a8 = ",theInputs['ggZZshape_a8']
            print "ggZZshape_a9 = ",theInputs['ggZZshape_a9']
                   
        
        bkg_ggzz = ROOT.RooggZZPdf_v2("bkg_ggzzTmp","bkg_ggzzTmp",CMS_zz4l_mass,CMS_ggzzbkg_a0,CMS_ggzzbkg_a1,CMS_ggzzbkg_a2,CMS_ggzzbkg_a3,CMS_ggzzbkg_a4,CMS_ggzzbkg_a5,CMS_ggzzbkg_a6,CMS_ggzzbkg_a7,CMS_ggzzbkg_a8,CMS_ggzzbkg_a9)
    
        ## Reducible backgrounds
        val_meanL_3P1F = float(theInputs['zjetsShape_mean_3P1F'])
        val_sigmaL_3P1F = float(theInputs['zjetsShape_sigma_3P1F'])
        val_normL_3P1F = float(theInputs['zjetsShape_norm_3P1F'])
        
        val_meanL_2P2F = float(theInputs['zjetsShape_mean_2P2F'])
        val_sigmaL_2P2F = float(theInputs['zjetsShape_sigma_2P2F'])
        val_normL_2P2F = float(theInputs['zjetsShape_norm_2P2F'])
        val_pol0_2P2F = float(theInputs['zjetsShape_pol0_2P2F'])
        val_pol1_2P2F = float(theInputs['zjetsShape_pol1_2P2F'])
        
        val_meanL_2P2F_2 = float(theInputs['zjetsShape_mean_2P2F_2e2mu'])
        val_sigmaL_2P2F_2 = float(theInputs['zjetsShape_sigma_2P2F_2e2mu'])
        val_normL_2P2F_2 = float(theInputs['zjetsShape_norm_2P2F_2e2mu'])

        if not self.bVBF:
            if (self.channel == self.ID_4mu):
                name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL_2P2F)
                name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL_2P2F)
                print "mean 4mu: ",mlZjet.getVal()
                print "sigma 4mu: ",slZjet.getVal()
                bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",CMS_zz4l_mass,mlZjet,slZjet)
            elif (self.channel == self.ID_4e):
                name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
                name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
                name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F)
                name = "p0Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",val_pol0_2P2F)
                name = "p1Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",val_pol1_2P2F)
                print "mean 2p2f 4e: ",mlZjet_2p2f.getVal()
                print "sigma 2p2f 4e: ",slZjet_2p2f.getVal()
                print "norm 2p2f 4e: ",nlZjet_2p2f.getVal()
                print "pol0 2p2f 4e: ",p0Zjet_2p2f.getVal()
                print "pol1 2p2f 4e: ",p1Zjet_2p2f.getVal()
                bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*@3*(1.+ TMath::Exp(@4+@5*@0))",RooArgList(CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f,nlZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))
                
                name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
                name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
                name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F)
                print "mean 3p1f 4e: ",mlZjet_3p1f.getVal()
                print "sigma 3p1f 4e: ",slZjet_3p1f.getVal()
                print "norm 3p1f 4e: ",nlZjet_3p1f.getVal()
                bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)

                bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f))

            elif (self.channel == self.ID_2e2mu):
                name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
                name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
                name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F)
                print "mean 2p2f 2mu2e: ",mlZjet_2p2f.getVal()
                print "sigma 2p2f 2mu2e: ",slZjet_2p2f.getVal()
                print "norm 2p2f 2mu2e: ",nlZjet_2p2f.getVal()
                bkg_zjets_2p2f = ROOT.RooLandau("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f",CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f)

                name = "mlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                mlZjet_2p2f_2 = ROOT.RooRealVar(name,"mean landau Zjet 2p2f 2e2mu",val_meanL_2P2F_2)
                name = "slZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                slZjet_2p2f_2 = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f 2e2mu",val_sigmaL_2P2F_2)
                name = "nlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                nlZjet_2p2f_2 = ROOT.RooRealVar(name,"norm landau Zjet 2p2f 2e2mu",val_normL_2P2F_2)
                print "mean 2p2f 2e2mu: ",mlZjet_2p2f_2.getVal()
                print "sigma 2p2f 2e2mu: ",slZjet_2p2f_2.getVal()
                print "norm 2p2f 2e2mu: ",nlZjet_2p2f_2.getVal()
                bkg_zjets_2p2f_2 = ROOT.RooLandau("bkg_zjetsTmp_2p2f_2","bkg_zjetsTmp_2p2f_2",CMS_zz4l_mass,mlZjet_2p2f_2,slZjet_2p2f_2)
                
                name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
                name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
                name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
                nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F)
                print "mean 3p1f 2mu2e: ",mlZjet_3p1f.getVal()
                print "sigma 3p1f 2mu2e: ",slZjet_3p1f.getVal()
                print "norm 3p1f 2mu2e: ",nlZjet_3p1f.getVal()
                bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)

                bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f,bkg_zjets_2p2f_2),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f,nlZjet_2p2f_2))

        else:
            if (self.channel == self.ID_4mu):
                name = "mlZjet_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL_2P2F)
                name = "slZjet_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts, self.VBFcat)
                slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL_2P2F)
                print "mean 4mu: ",mlZjet.getVal()
                print "sigma 4mu: ",slZjet.getVal()
                bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",CMS_zz4l_mass,mlZjet,slZjet)
            elif (self.channel == self.ID_4e):
                name = "mlZjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts, self.VBFcat)
                mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
                name = "slZjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts, self.VBFcat)
                slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
                name = "nlZjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts, self.VBFcat)
                nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F)
                name = "p0Zjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts, self.VBFcat)
                p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",val_pol0_2P2F)
                name = "p1Zjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts, self.VBFcat)
                p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",val_pol1_2P2F)
                print "mean 2p2f 4e: ",mlZjet_2p2f.getVal()
                print "sigma 2p2f 4e: ",slZjet_2p2f.getVal()
                print "norm 2p2f 4e: ",nlZjet_2p2f.getVal()
                print "pol0 2p2f 4e: ",p0Zjet_2p2f.getVal()
                print "pol1 2p2f 4e: ",p1Zjet_2p2f.getVal()
                bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*@3*(1.+ TMath::Exp(@4+@5*@0))",RooArgList(CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f,nlZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))
                
                name = "mlZjet_3p1f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
                name = "slZjet_3p1f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
                name = "nlZjet_3p1f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F)
                print "mean 3p1f 4e: ",mlZjet_3p1f.getVal()
                print "sigma 3p1f 4e: ",slZjet_3p1f.getVal()
                print "norm 3p1f 4e: ",nlZjet_3p1f.getVal()
                bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)

                bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f))

            elif (self.channel == self.ID_2e2mu):
                name = "mlZjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
                name = "slZjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
                name = "nlZjet_2p2f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F)
                print "mean 2p2f 2mu2e: ",mlZjet_2p2f.getVal()
                print "sigma 2p2f 2mu2e: ",slZjet_2p2f.getVal()
                print "norm 2p2f 2mu2e: ",nlZjet_2p2f.getVal()
                bkg_zjets_2p2f = ROOT.RooLandau("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f",CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f)

                name = "mlZjet_2p2f_2_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                mlZjet_2p2f_2 = ROOT.RooRealVar(name,"mean landau Zjet 2p2f 2e2mu",val_meanL_2P2F_2)
                name = "slZjet_2p2f_2_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                slZjet_2p2f_2 = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f 2e2mu",val_sigmaL_2P2F_2)
                name = "nlZjet_2p2f_2_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                nlZjet_2p2f_2 = ROOT.RooRealVar(name,"norm landau Zjet 2p2f 2e2mu",val_normL_2P2F_2)
                print "mean 2p2f 2e2mu: ",mlZjet_2p2f_2.getVal()
                print "sigma 2p2f 2e2mu: ",slZjet_2p2f_2.getVal()
                print "norm 2p2f 2e2mu: ",nlZjet_2p2f_2.getVal()
                bkg_zjets_2p2f_2 = ROOT.RooLandau("bkg_zjetsTmp_2p2f_2","bkg_zjetsTmp_2p2f_2",CMS_zz4l_mass,mlZjet_2p2f_2,slZjet_2p2f_2)
                
                name = "mlZjet_3p1f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
                name = "slZjet_3p1f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
                name = "nlZjet_3p1f_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F)
                print "mean 3p1f 2mu2e: ",mlZjet_3p1f.getVal()
                print "sigma 3p1f 2mu2e: ",slZjet_3p1f.getVal()
                print "norm 3p1f 2mu2e: ",nlZjet_3p1f.getVal()
                bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)

                bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f,bkg_zjets_2p2f_2),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f,nlZjet_2p2f_2))

 
	bkg_qqzzErr = ROOT.RooProdPdf("bkg_qqzzErr","bkg_qqzzErr", ROOT.RooArgSet(bkg_qqzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrZZ), ROOT.RooArgSet(RelErr)));
	bkg_ggzzErr = ROOT.RooProdPdf("bkg_ggzzErr","bkg_ggzzErr", ROOT.RooArgSet(bkg_ggzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrZZ), ROOT.RooArgSet(RelErr)));
	bkg_zjetsErr = ROOT.RooProdPdf("bkg_zjetsErr","bkg_zjetsErr", ROOT.RooArgSet(bkg_zjets), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrZX), ROOT.RooArgSet(RelErr)));


      ##-------------------bVBF background ----------------------
        if(self.bVBF):
            bkg2d_qqZZ_Fisher = ROOT.RooProdPdf()
            bkg2d_ggZZ_Fisher = ROOT.RooProdPdf()
            bkg2d_ZX_Fisher = ROOT.RooProdPdf()
            bkg2d_qqZZ_Pt = ROOT.RooProdPdf()
            bkg2d_ggZZ_Pt = ROOT.RooProdPdf()
            bkg2d_ZX_Pt = ROOT.RooProdPdf()
        if(self.bVBF and self.VBFcat):
            
            qqZZtempFileName = "{0}/qqZZ_fisher.root".format(self.templateDir)
            qqZZtempFile = ROOT.TFile(qqZZtempFileName)
            qqZZtemplate = qqZZtempFile.Get("h_Fisher")
            qqZZtemplate_Up = qqZZtempFile.Get("h_Fisher_up")
            qqZZtemplate_Dn = qqZZtempFile.Get("h_Fisher_dn")

            Fisher_qqZZ_dataHist = ROOT.RooDataHist("temp_qqZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_qqZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),qqZZtemplate)
            Fisher_qqZZ_dataHist_Up = ROOT.RooDataHist("temp_qqZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_qqZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),qqZZtemplate_Up)
            Fisher_qqZZ_dataHist_Dn = ROOT.RooDataHist("temp_qqZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_qqZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),qqZZtemplate_Dn)
            
            TemplateName = "FisherTempDataHist_qqzz_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_qqZZ = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_qqZZ_dataHist)
            TemplateName = "FisherTempDataHist_qqzz_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_qqZZ_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_qqZZ_dataHist_Up)
            TemplateName = "FisherTempDataHist_qqzz_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_qqZZ_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_qqZZ_dataHist_Dn)

            ggZZtempFileName = "{0}/ggZZ_fisher.root".format(self.templateDir)
            ggZZtempFile = ROOT.TFile(ggZZtempFileName)
            ggZZtemplate = ggZZtempFile.Get("h_Fisher")
            
            Fisher_ggZZ_dataHist = ROOT.RooDataHist("temp_ggZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ggZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ggZZtemplate)
            
            TemplateName = "FisherTempDataHist_ggzz_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ggZZ = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ggZZ_dataHist)

            ZXtempFileName = "{0}/Z+X_fisher.root".format(self.templateDir)
            ZXtempFile = ROOT.TFile(ZXtempFileName)
            ZXtemplate = ZXtempFile.Get("h_Fisher")
            
            Fisher_ZX_dataHist = ROOT.RooDataHist("temp_ZX_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"temp_ZX_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,VD),ZXtemplate)
            
            TemplateName = "FisherTempDataHist_zx_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplatePdf_ZX = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,VD),Fisher_ZX_dataHist)

            FisherList_qqZZ = ROOT.RooArgList()  
            FisherList_ggZZ = ROOT.RooArgList()
            FisherList_ZX  = ROOT.RooArgList()
            FisherList_Zjets  = ROOT.RooArgList()

            if(self.FisherMorph):
                FisherList_qqZZ.add(FisherTemplatePdf_qqZZ)
                FisherList_qqZZ.add(FisherTemplatePdf_qqZZ_Up)
                FisherList_qqZZ.add(FisherTemplatePdf_qqZZ_Dn) 
                
                FisherList_ggZZ.add(FisherTemplatePdf_ggZZ)
                
                FisherList_ZX.add(FisherTemplatePdf_ZX)  
            
            
            else:
            
                FisherList_qqZZ.add(FisherTemplatePdf_qqZZ)
                FisherList_ggZZ.add(FisherTemplatePdf_ggZZ)
                FisherList_ZX.add(FisherTemplatePdf_ZX)
                

            morphFisherVarName = "CMS_zz4l_qqZZ_Fisher_sys"
            alphaMorphFisher_qqZZ = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            morphFisherVarName = "CMS_zz4l_ggZZ_Fisher_sys"
            alphaMorphFisher_ggZZ = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            morphFisherVarName = "CMS_zz4l_ZX_Fisher_sys"
            alphaMorphFisher_ZX = ROOT.RooRealVar(morphFisherVarName,morphFisherVarName,0,-3,3)
            if(self.FisherMorph):
                alphaMorphFisher_qqZZ.setConstant(False)
            else:
                alphaMorphFisher_qqZZ.setConstant(True)
            alphaMorphFisher_ggZZ.setConstant(True)
            alphaMorphFisher_ZX.setConstant(True)
        
            morphVarListFisher_qqZZ = ROOT.RooArgList()
            morphVarListFisher_ggZZ = ROOT.RooArgList()
            morphVarListFisher_ZX = ROOT.RooArgList()
                
            if(self.FisherMorph):
                morphVarListFisher_qqZZ.add(alphaMorphFisher_qqZZ)
        
            TemplateName = "FisherTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_qqZZ = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_qqZZ,morphVarListFisher_qqZZ,1.0,1)
        
            TemplateName = "FisherTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_ggZZ = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_ggZZ,morphVarListFisher_ggZZ,1.0,1)
        
            TemplateName = "FisherTemplateMorphPdf_zx_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            FisherTemplateMorphPdf_ZX = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,VD,true,FisherList_ZX,morphVarListFisher_ZX,1.0,1)

            bkg2d_qqZZ_Fisher = ROOT.RooProdPdf("bkg2d_qqZZ_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"bkg2d_qqZZ_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(bkg_qqzzErr,bkg_qqzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_qqZZ),ROOT.RooArgSet(VD)))
            bkg2d_ggZZ_Fisher = ROOT.RooProdPdf("bkg2d_ggZZ_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"bkg2d_ggZZ_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(bkg_ggzzErr,bkg_ggzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ggZZ),ROOT.RooArgSet(VD)))
            bkg2d_ZX_Fisher = ROOT.RooProdPdf("bkg2d_ZX_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"bkg2d_ZX_Fisher_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(bkg_zjetsErr,bkg_zjets,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(FisherTemplateMorphPdf_ZX),ROOT.RooArgSet(VD)))
        
            
        if(self.bVBF and not self.VBFcat):
            qqZZtempPtFileName = "{0}/Pt_mZZ_zz_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            qqZZtempPtFile = ROOT.TFile(qqZZtempPtFileName)
            qqZZtemplatePt = qqZZtempPtFile.Get("h_Ptmzz_mzz")            
            qqZZtemplatePt_Up = ggHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            qqZZtemplatePt_Dn = ggHtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_qqZZ_dataHist = ROOT.RooDataHist("tempPt_qqZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_qqZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),qqZZtemplatePt)
            Pt_qqZZ_dataHist_Up = ROOT.RooDataHist("tempPt_qqZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_qqZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),qqZZtemplatePt_Up)
            Pt_qqZZ_dataHist_Dn = ROOT.RooDataHist("tempPt_qqZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_qqZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),qqZZtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_qqZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_qqZZ = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_qqZZ_dataHist)
            TemplateName = "PtTempDataHist_qqZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_qqZZ_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_qqZZ_dataHist_Up)
            TemplateName = "PtTempDataHist_qqZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_qqZZ_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_qqZZ_dataHist_Dn)
            
            ggZZtempPtFileName = "{0}/Pt_mZZ_ggzz_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            ggZZtempPtFile = ROOT.TFile(ggZZtempPtFileName)
            ggZZtemplatePt = ggZZtempPtFile.Get("h_Ptmzz_mzz")
            ggZZtemplatePt_Up = ggZZtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            ggZZtemplatePt_Dn = ggZZtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_ggZZ_dataHist = ROOT.RooDataHist("tempPt_ggZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ggZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ggZZtemplatePt)
            Pt_ggZZ_dataHist_Up = ROOT.RooDataHist("tempPt_ggZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ggZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ggZZtemplatePt_Up)
            Pt_ggZZ_dataHist_Dn = ROOT.RooDataHist("tempPt_ggZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ggZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ggZZtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_ggZZ_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ggZZ = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ggZZ_dataHist)
            TemplateName = "PtTempDataHist_ggZZ_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ggZZ_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ggZZ_dataHist_Up)
            TemplateName = "PtTempDataHist_ggZZ_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ggZZ_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ggZZ_dataHist_Dn)

            ZXtempPtFileName = "{0}/Pt_mZZ_zx_{1}_{2:d}TeV.root".format(self.templateDir,self.appendName,int(self.sqrts))
            ZXtempPtFile = ROOT.TFile(ZXtempPtFileName)
            ZXtemplatePt = ZXtempPtFile.Get("h_Ptmzz_mzz")
            ZXtemplatePt_Up = ZXtempPtFile.Get("h_Ptmzz_mzz_OneSyst_up")
            ZXtemplatePt_Dn = ZXtempPtFile.Get("h_Ptmzz_mzz_OneSyst_down")            
            
            Pt_ZX_dataHist = ROOT.RooDataHist("tempPt_ZX_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ZX_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ZXtemplatePt)
            Pt_ZX_dataHist_Up = ROOT.RooDataHist("tempPt_ZX_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ZX_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ZXtemplatePt_Up)
            Pt_ZX_dataHist_Dn = ROOT.RooDataHist("tempPt_ZX_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"tempPt_ZX_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgList(CMS_zz4l_mass,pt),ZXtemplatePt_Dn)
            
            TemplateName = "PtTempDataHist_ZX_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ZX = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ZX_dataHist)
            TemplateName = "PtTempDataHist_ZX_Up_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ZX_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ZX_dataHist_Up)
            TemplateName = "PtTempDataHist_ZX_Dn_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplatePdf_ZX_Dn = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,pt),Pt_ZX_dataHist_Dn)

            PtList_qqZZ = ROOT.RooArgList()  
            PtList_ggZZ = ROOT.RooArgList()
            PtList_ZX  = ROOT.RooArgList()

            if(self.PtMorph):
                PtList_qqZZ.add(PtTemplatePdf_qqZZ)
                PtList_qqZZ.add(PtTemplatePdf_qqZZ_Up)
                PtList_qqZZ.add(PtTemplatePdf_qqZZ_Dn)  
                
                PtList_ggZZ.add(PtTemplatePdf_ggZZ)
                PtList_ggZZ.add(PtTemplatePdf_ggZZ_Up)
                PtList_ggZZ.add(PtTemplatePdf_ggZZ_Dn)  
                
                PtList_ZX.add(PtTemplatePdf_ZX)
                PtList_ZX.add(PtTemplatePdf_ZX_Up)
                PtList_ZX.add(PtTemplatePdf_ZX_Dn)  
            
            
            else:
            
                PtList_qqZZ.add(PtTemplatePdf_qqZZ)
                PtList_ggZZ.add(PtTemplatePdf_ggZZ)
                PtList_ZX.add(PtTemplatePdf_ZX)
                

            morphPtVarName = "CMS_zz4l_qqZZ_Pt_sys"
            alphaMorphPt_qqZZ = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            morphPtVarName = "CMS_zz4l_ggZZ_Pt_sys"
            alphaMorphPt_ggZZ = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            morphPtVarName = "CMS_zz4l_ZX_Pt_sys"
            alphaMorphPt_ZX = ROOT.RooRealVar(morphPtVarName,morphPtVarName,0,-3,3)
            if(self.PtMorph):
                alphaMorphPt_qqZZ.setConstant(False)
                alphaMorphPt_ggZZ.setConstant(False)
                alphaMorphPt_ZX.setConstant(False)
            else:
                alphaMorphPt_qqZZ.setConstant(True)
                alphaMorphPt_ggZZ.setConstant(True)
                alphaMorphPt_ZX.setConstant(True)
        
            morphVarListPt_qqZZ = ROOT.RooArgList()
            morphVarListPt_ggZZ = ROOT.RooArgList()
            morphVarListPt_ZX = ROOT.RooArgList()
                
            if(self.PtMorph):
                morphVarListPt_qqZZ.add(alphaMorphPt_qqZZ)
                morphVarListPt_ggZZ.add(alphaMorphPt_ggZZ)
                morphVarListPt_ZX.add(alphaMorphPt_ZX)
        
            PtTemplateName = "PtTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_qqZZ = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_qqZZ,morphVarListPt_qqZZ,1.0,1)
        
            PtTemplateName = "PtTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_ggZZ = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_ggZZ,morphVarListPt_ggZZ,1.0,1)
        
            PtTemplateName = "PtTemplateMorphPdf_zx_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            PtTemplateMorphPdf_ZX = ROOT.FastVerticalInterpHistPdf2D(PtTemplateName,PtTemplateName,CMS_zz4l_mass,pt,true,PtList_ZX,morphVarListPt_ZX,1.0,1)

            bkg2d_qqZZ_Pt = ROOT.RooProdPdf("bkg2d_qqZZ_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"bkg2d_qqZZ_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(bkg_qqzzErr,bkg_qqzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_qqZZ),ROOT.RooArgSet(pt)))
            bkg2d_ggZZ_Pt = ROOT.RooProdPdf("bkg2d_ggZZ_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"bkg2d_ggZZ_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(bkg_ggzzErr,bkg_ggzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ggZZ),ROOT.RooArgSet(pt)))
            bkg2d_ZX_Pt = ROOT.RooProdPdf("bkg2d_ZX_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),"bkg2d_ZX_Pt_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat),ROOT.RooArgSet(self.getVariable(bkg_zjetsErr,bkg_zjets,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(PtTemplateMorphPdf_ZX),ROOT.RooArgSet(pt)))
            
      ## ----------------- 2D BACKGROUND SHAPES --------------- ##
        if self.useMEKDTemplates:
            templateBkgName = "{0}/Dbackground_ZX_{1}.root".format(self.templateDir, self.appendName)
        else:
            templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
            
        print templateBkgName, "file used for ZZ"
        
        bkgTempFile = ROOT.TFile(templateBkgName)
        bkgTemplate = bkgTempFile.Get("h_mzzD")
        
        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate)
        
        
        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        ggbkgTempFile = ROOT.TFile(templateggBkgName)
        ggbkgTemplate = ggbkgTempFile.Get("h_mzzD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),ggbkgTemplate)
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist)
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),ggbkgTempDataHist)

        if self.useMEKDTemplates:
            templatezxBkgName = "{0}/Dbackground_ZX_{1}.root".format(self.templateDir, self.appendName)
        else:
            templatezxBkgName = "{0}/Dbackground_ZX_{1}.root".format(self.templateDir ,self.appendName)
            
        print templatezxBkgName, "file used for ZX"

        zxbkgTempFile = ROOT.TFile(templatezxBkgName)
        zxbkgTemplate = zxbkgTempFile.Get("h_mzzD")
        zxbkgTemplate_Up = zxbkgTempFile.Get("h_mzzD_up")
        zxbkgTemplate_Down = zxbkgTempFile.Get("h_mzzD_dn")

        TemplateName = "zjetsTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        zxbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),zxbkgTemplate)
        TemplateName = "zxbkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        zxbkgTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),zxbkgTemplate_Up)
        TemplateName = "zxbkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        zxbkgTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),zxbkgTemplate_Down)
        
        TemplateName = "zxbkgTemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),zxbkgTempDataHist)
        TemplateName = "zxbkgTemplatePdf_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),zxbkgTempDataHist_Up)
        TemplateName = "zxbkgTemplatePdf_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),zxbkgTempDataHist_Down)
        
        funcList_zjets = ROOT.RooArgList()  
        morphBkgVarName = "CMS_zz4l_bkgMELA"    
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()
        
        if(self.bkgMorph):
            funcList_zjets.add(bkgTemplatePdf_zjets)
            funcList_zjets.add(bkgTemplatePdf_zjets_Up)
            funcList_zjets.add(bkgTemplatePdf_zjets_Down)  
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)  
        else:
            funcList_zjets.add(bkgTemplatePdf_zjets)
            alphaMorphBkg.setConstant(True)


        TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_zjets,morphVarListBkg,1.0,1)

	####  ----------------------- mekd  parametrized double gaussian stuffs  -------------------------
	if theMEKD: 
		name = "mekd_qqZZ_a0_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		print "mekd_qqZZ_a0_shape=",theInputs['mekd_qqZZ_a0_shape'] 
		print "mekd_sig_a0_shape=",theInputs['mekd_sig_a0_shape'] 
		mekd_qqZZ_a0 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_qqZZ_a0_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_qqZZ_a1_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_qqZZ_a1 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_qqZZ_a1_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_qqZZ_a2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_qqZZ_a2 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_qqZZ_a2_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_qqZZ_a3_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_qqZZ_a3 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_qqZZ_a3_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		name = "mekd_qqZZ_a4_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
		mekd_qqZZ_a4 = ROOT.RooFormulaVar(name,"("+theInputs['mekd_qqZZ_a4_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
		bkgTemplateMorphPdf_qqzz = ROOT.RooGenericPdf("mekd_qqZZ", "mekd_qqZZ", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_qqZZ_a0, mekd_qqZZ_a1, mekd_qqZZ_a2, mekd_qqZZ_a3, mekd_qqZZ_a4))
		bkgTemplateMorphPdf_ggzz = ROOT.RooGenericPdf("mekd_ggZZ", "mekd_ggZZ", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_qqZZ_a0, mekd_qqZZ_a1, mekd_qqZZ_a2, mekd_qqZZ_a3, mekd_qqZZ_a4))
		bkgTemplateMorphPdf_zjets= ROOT.RooGenericPdf("mekd_zjets", "mekd_zjets", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_qqZZ_a0, mekd_qqZZ_a1, mekd_qqZZ_a2, mekd_qqZZ_a3, mekd_qqZZ_a4))
		m = 100
		while m >= 100 and m < 150:
			CMS_zz4l_mass.setVal(m)
			m = m + 0.1
			if mekd_qqZZ_a2.getVal() < 0 : print m, mekd_qqZZ_a2.getVal() 
			if mekd_qqZZ_a2.getVal() > 1 : print m, mekd_qqZZ_a2.getVal() 
		print "\n \n mekd_qqZZ_a1 channel ",self.channel
		m = 100
		while m >= 100 and m < 150:
			CMS_zz4l_mass.setVal(m)
			m = m + 0.1
			if mekd_qqZZ_a1.getVal() <= 0 : print m, mekd_qqZZ_a1.getVal() 
		print "\n \n mekd_qqZZ_a4 channel ",self.channel
		m = 100
		while m >= 100 and m < 150:
			CMS_zz4l_mass.setVal(m)
			m = m + 0.1
			if mekd_qqZZ_a4.getVal() <= 0 : print m, mekd_qqZZ_a4.getVal() 
	####  ----------------------- end mekd -----------------------------------------------------------
        bkg2d_qqzz = ROOT.RooProdPdf("bkg2d_qqzz","bkg2d_qqzz",ROOT.RooArgSet(self.getVariable(bkg_qqzzErr,bkg_qqzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        bkg2d_ggzz = ROOT.RooProdPdf("bkg2d_ggzz","bkg2d_ggzz",ROOT.RooArgSet(self.getVariable(bkg_ggzzErr,bkg_ggzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))
        bkg2d_zjets = ROOT.RooProdPdf("bkg2d_zjets","bkg2d_zjets",ROOT.RooArgSet(self.getVariable(bkg_zjetsErr,bkg_zjets,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(self.getVariable(MEKD,D,self.bMEKD))))

        if(self.bVBF):
            if(self.Use3D):
                bkg2d_qqZZ_VBF_KD = ROOT.RooProdPdf("bkg2d_qqZZ_VBF_KD","bkg2d_qqZZ_VBF_KD",ROOT.RooArgSet(self.getVariable(bkg2d_qqZZ_Fisher,bkg2d_qqZZ_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(D)))
                bkg2d_ggZZ_VBF_KD = ROOT.RooProdPdf("bkg2d_ggZZ_VBF_KD","bkg2d_ggZZ_VBF_KD",ROOT.RooArgSet(self.getVariable(bkg2d_ggZZ_Fisher,bkg2d_ggZZ_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(D)))
                bkg2d_ZX_VBF_KD = ROOT.RooProdPdf("bkg2d_ZX_VBF_KD","bkg2d_ZX_VBF_KD",ROOT.RooArgSet(self.getVariable(bkg2d_ZX_Fisher,bkg2d_ZX_Pt,self.VBFcat)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(D)))

        ## ----------------- SUPERMELA BACKGROUND SHAPES --------------- ##
        
        #templateSDBkgName = "{0}/Dbackground_qqZZ_superMELA_{1}.root".format(self.templateDir ,self.appendName) 
        #bkgTempSDFile = ROOT.TFile(templateSDBkgName)
        #bkgTemplateSD = bkgTempSDFile.Get("hSuperD_bkg") 
        
        #TemplateSDName = "bkgTempSDDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        #bkgTempSDDataHist = ROOT.RooDataHist(TemplateSDName,TemplateSDName,ROOT.RooArgList(SD),bkgTemplateSD)
        
        #TemplateSDName = "bkgTemplateSDPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        #bkgTemplateSDPdf_qqzz = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)

        #TemplateSDName = "bkgTemplateSDPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        #bkgTemplateSDPdf_ggzz = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)
        #TemplateSDName = "bkgTemplateSDPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        #bkgTemplateSDPdf_zjets = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)
        

        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##
        if not self.bVBF:
            canv_name = "czz_{0}_{1}".format(self.mH,self.appendName)
        else:
            canv_name = "czz_{0}_{1}_{2}".format(self.mH,self.appendName,self.VBFcat)
        czz = ROOT.TCanvas( canv_name, canv_name, 750, 700 )
        czz.cd()
        zzframe_s = CMS_zz4l_mass.frame(45)
        if self.bUseCBnoConvolution: super(RooDoubleCB,signalCB_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        elif self.isHighMass : super(ROOT.RooFFTConvPdf,sig_ggH_HM).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        else : super(ROOT.RooFFTConvPdf,sig_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        super(ROOT.RooqqZZPdf_v2,bkg_qqzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )
        super(ROOT.RooggZZPdf_v2,bkg_ggzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(6) )
        super(ROOT.RooAbsPdf,bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
        zzframe_s.Draw()
        if not self.bVBF:
            figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
        else:
            figName = "{0}/figs/mzz_{1}_{2}_{3}.png".format(self.outputDir, self.mH, self.appendName,self.VBFcat)
        czz.SaveAs(figName)
        del czz
        
        ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        
        CMS_zz4l_mass.setRange("shape",self.low_M,self.high_M)
        
        fr_low_M = self.low_M
        fr_high_M = self.high_M        
        if (self.mH >= 450): 
            fr_low_M = 100
            fr_high_M = 1000
        if (self.mH >= 750):
            fr_low_M = 100
            fr_high_M = 1400
            

        CMS_zz4l_mass.setRange("fullrangesignal",fr_low_M,fr_high_M)
        CMS_zz4l_mass.setRange("fullrange",100,1400)
        
       
        if not self.bVBF:
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
            rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
            rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
            rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
            rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
            rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
            rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
            rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
            rrvg1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g1'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
            rrvg2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g2'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
            rrvg3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g3'])

            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
            rrva1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa1'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
            rrva2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa2'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
            rrva3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa3'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
            rrva4_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa4'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
            rrvb1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb1'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
            rrvb2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb2'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
            rrvb3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb3'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
            rrvg1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg1'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
            rrvg2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg2'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
            rrvg3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg3'])

            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
            rrva1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa1'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
            rrva2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa2'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
            rrva3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa3'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
            rrva4_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa4'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
            rrvb1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb1'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
            rrvb2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb2'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
            rrvb3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb3'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
            rrvg1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg1'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
            rrvg2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg2'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
            rrvg3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg3'])

            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
            rrva1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa1'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
            rrva2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa2'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
            rrva3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa3'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
            rrva4_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa4'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
            rrvb1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb1'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
            rrvb2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb2'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
            rrvb3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb3'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
            rrvg1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg1'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
            rrvg2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg2'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
            rrvg3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg3'])

            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
            rrva1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa1'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
            rrva2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa2'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
            rrva3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa3'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
            rrva4_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa4'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
            rrvb1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb1'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
            rrvb2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb2'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
            rrvb3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb3'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
            rrvg1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg1'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
            rrvg2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg2'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
            rrvg3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg3'])
        else:
            
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts, self.VBFcat)
            rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts, self.VBFcat)
            rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts, self.VBFcat)
            rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts, self.VBFcat)
            rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts, self.VBFcat)
            rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts, self.VBFcat)
            rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts, self.VBFcat)
            rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts, self.VBFcat)
            rrvg1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g1'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts, self.VBFcat)
            rrvg2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g2'])
            sigEffName = "hzz4lggHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts, self.VBFcat)
            rrvg3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g3'])

            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts, self.VBFcat)
            rrva1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa1'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts, self.VBFcat)
            rrva2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa2'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts, self.VBFcat)
            rrva3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa3'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts, self.VBFcat)
            rrva4_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa4'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts, self.VBFcat)
            rrvb1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb1'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts, self.VBFcat)
            rrvb2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb2'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts, self.VBFcat)
            rrvb3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb3'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts, self.VBFcat)
            rrvg1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg1'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts, self.VBFcat)
            rrvg2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg2'])
            sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts, self.VBFcat)
            rrvg3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg3'])

            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts, self.VBFcat)
            rrva1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa1'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts, self.VBFcat)
            rrva2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa2'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts, self.VBFcat)
            rrva3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa3'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts, self.VBFcat)
            rrva4_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa4'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts, self.VBFcat)
            rrvb1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb1'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts, self.VBFcat)
            rrvb2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb2'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts, self.VBFcat)
            rrvb3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb3'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts, self.VBFcat)
            rrvg1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg1'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts, self.VBFcat)
            rrvg2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg2'])
            sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts, self.VBFcat)
            rrvg3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg3'])

            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts, self.VBFcat)
            rrva1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa1'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts, self.VBFcat)
            rrva2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa2'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts, self.VBFcat)
            rrva3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa3'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts, self.VBFcat)
            rrva4_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa4'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts, self.VBFcat)
            rrvb1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb1'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts, self.VBFcat)
            rrvb2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb2'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts, self.VBFcat)
            rrvb3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb3'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts, self.VBFcat)
            rrvg1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg1'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts, self.VBFcat)
            rrvg2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg2'])
            sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts, self.VBFcat)
            rrvg3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg3'])

            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts, self.VBFcat)
            rrva1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa1'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts, self.VBFcat)
            rrva2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa2'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts, self.VBFcat)
            rrva3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa3'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts, self.VBFcat)
            rrva4_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa4'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts, self.VBFcat)
            rrvb1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb1'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts, self.VBFcat)
            rrvb2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb2'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts, self.VBFcat)
            rrvb3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb3'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts, self.VBFcat)
            rrvg1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg1'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts, self.VBFcat)
            rrvg2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg2'])
            sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts, self.VBFcat)
            rrvg3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg3'])

        if(DEBUG):
            print "sigEff_a1 = ",theInputs['sigEff_a1']
            print "sigEff_a2 = ",theInputs['sigEff_a2']
            print "sigEff_a3 = ",theInputs['sigEff_a3']
            print "sigEff_a4 = ",theInputs['sigEff_a4']
            print "sigEff_b1 = ",theInputs['sigEff_b1']
            print "sigEff_b2 = ",theInputs['sigEff_b2']
            print "sigEff_b3 = ",theInputs['sigEff_b3']
            print "sigEff_g1 = ",theInputs['sigEff_g1']
            print "sigEff_g2 = ",theInputs['sigEff_g2']
            print "sigEff_g3 = ",theInputs['sigEff_g3']

            print "sigEff_qqHa1 = ",theInputs['sigEff_qqHa1']
            print "sigEff_qqHa2 = ",theInputs['sigEff_qqHa2']
            print "sigEff_qqHa3 = ",theInputs['sigEff_qqHa3']
            print "sigEff_qqHa4 = ",theInputs['sigEff_qqHa4']
            print "sigEff_qqHb1 = ",theInputs['sigEff_qqHb1']
            print "sigEff_qqHb2 = ",theInputs['sigEff_qqHb2']
            print "sigEff_qqHb3 = ",theInputs['sigEff_qqHb3']
            print "sigEff_qqHg1 = ",theInputs['sigEff_qqHg1']
            print "sigEff_qqHg2 = ",theInputs['sigEff_qqHg2']
            print "sigEff_qqHg3 = ",theInputs['sigEff_qqHg3']

            print "sigEff_ZHa1 = ",theInputs['sigEff_ZHa1']
            print "sigEff_ZHa2 = ",theInputs['sigEff_ZHa2']
            print "sigEff_ZHa3 = ",theInputs['sigEff_ZHa3']
            print "sigEff_ZHa4 = ",theInputs['sigEff_ZHa4']
            print "sigEff_ZHb1 = ",theInputs['sigEff_ZHb1']
            print "sigEff_ZHb2 = ",theInputs['sigEff_ZHb2']
            print "sigEff_ZHb3 = ",theInputs['sigEff_ZHb3']
            print "sigEff_ZHg1 = ",theInputs['sigEff_ZHg1']
            print "sigEff_ZHg2 = ",theInputs['sigEff_ZHg2']
            print "sigEff_ZHg3 = ",theInputs['sigEff_ZHg3']

            print "sigEff_WHa1 = ",theInputs['sigEff_WHa1']
            print "sigEff_WHa2 = ",theInputs['sigEff_WHa2']
            print "sigEff_WHa3 = ",theInputs['sigEff_WHa3']
            print "sigEff_WHa4 = ",theInputs['sigEff_WHa4']
            print "sigEff_WHb1 = ",theInputs['sigEff_WHb1']
            print "sigEff_WHb2 = ",theInputs['sigEff_WHb2']
            print "sigEff_WHb3 = ",theInputs['sigEff_WHb3']
            print "sigEff_WHg1 = ",theInputs['sigEff_WHg1']
            print "sigEff_WHg2 = ",theInputs['sigEff_WHg2']
            print "sigEff_WHg3 = ",theInputs['sigEff_WHg3']

            print "sigEff_ttHa1 = ",theInputs['sigEff_ttHa1']
            print "sigEff_ttHa2 = ",theInputs['sigEff_ttHa2']
            print "sigEff_ttHa3 = ",theInputs['sigEff_ttHa3']
            print "sigEff_ttHa4 = ",theInputs['sigEff_ttHa4']
            print "sigEff_ttHb1 = ",theInputs['sigEff_ttHb1']
            print "sigEff_ttHb2 = ",theInputs['sigEff_ttHb2']
            print "sigEff_ttHb3 = ",theInputs['sigEff_ttHb3']
            print "sigEff_ttHg1 = ",theInputs['sigEff_ttHg1']
            print "sigEff_ttHg2 = ",theInputs['sigEff_ttHg2']
            print "sigEff_ttHg3 = ",theInputs['sigEff_ttHg3']
           

        if not self.bVBF:
            sigEffName_ggH = "hzz4lggHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            sigEffName_qqH = "hzz4lqqHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            sigEffName_WH = "hzz4lWHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            sigEffName_ZH = "hzz4lZHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            sigEffName_ttH = "hzz4lttHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        else:
            sigEffName_ggH = "hzz4lggHeff_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
            sigEffName_qqH = "hzz4lqqHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.VBFcat)
            sigEffName_WH = "hzz4lWHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.VBFcat)
            sigEffName_ZH = "hzz4lZHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.VBFcat)
            sigEffName_ttH = "hzz4lttHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.VBFcat)

        listEff = ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH)
        listEff.add(rrvg1)
        listEff.add(rrvg2)
        listEff.add(rrvg3)

        listEff_qqh = ROOT.RooArgList(rrva1_qqh,rrva2_qqh,rrva3_qqh,rrva4_qqh,rrvb1_qqh,rrvb2_qqh,rrvb3_qqh,self.MH)
        listEff_qqh.add(rrvg1_qqh)
        listEff_qqh.add(rrvg2_qqh)
        listEff_qqh.add(rrvg3_qqh)

        listEff_wh = ROOT.RooArgList(rrva1_wh,rrva2_wh,rrva3_wh,rrva4_wh,rrvb1_wh,rrvb2_wh,rrvb3_wh,self.MH)
        listEff_wh.add(rrvg1_wh)
        listEff_wh.add(rrvg2_wh)
        listEff_wh.add(rrvg3_wh)

        listEff_zh = ROOT.RooArgList(rrva1_zh,rrva2_zh,rrva3_zh,rrva4_zh,rrvb1_zh,rrvb2_zh,rrvb3_zh,self.MH)
        listEff_zh.add(rrvg1_zh)
        listEff_zh.add(rrvg2_zh)
        listEff_zh.add(rrvg3_zh)

        listEff_tth = ROOT.RooArgList(rrva1_tth,rrva2_tth,rrva3_tth,rrva4_tth,rrvb1_tth,rrvb2_tth,rrvb3_tth,self.MH)
        listEff_tth.add(rrvg1_tth)
        listEff_tth.add(rrvg2_tth)
        listEff_tth.add(rrvg3_tth)
        
        rfvSigEff_ggH = ROOT.RooFormulaVar(sigEffName_ggH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff) #ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH,rrvg1,rrvg2,rrvg3))

        rfvSigEff_qqH = ROOT.RooFormulaVar(sigEffName_qqH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_qqh)
        rfvSigEff_ZH = ROOT.RooFormulaVar(sigEffName_ZH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_zh)
        rfvSigEff_WH = ROOT.RooFormulaVar(sigEffName_WH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_wh)
        rfvSigEff_ttH = ROOT.RooFormulaVar(sigEffName_ttH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",listEff_tth)
        #from TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., xMax);
        
        ## following printout is needed ,  dont remove it
        print " @@@@@@@@ ggHeff ",rfvSigEff_ggH.getVal()
        print " @@@@@@@@ qqHeff ",rfvSigEff_qqH.getVal()
        print " @@@@@@@@ ZHeff ",rfvSigEff_ZH.getVal()
        print " @@@@@@@@ WHeff ",rfvSigEff_WH.getVal()
        print " @@@@@@@@ ttHeff ",rfvSigEff_ttH.getVal()
    
        CS_ggH = myCSW.HiggsCS(1,self.mH,self.sqrts)
        CS_VBF = myCSW.HiggsCS(2,self.mH,self.sqrts)
        CS_WH = myCSW.HiggsCS(3,self.mH,self.sqrts)
        CS_ZH = myCSW.HiggsCS(4,self.mH,self.sqrts)
        CS_ttH = myCSW.HiggsCS(5,self.mH,self.sqrts)
    
        BRH2e2mu = myCSW.HiggsBR(13,self.mH)
        BRH4mu = myCSW.HiggsBR(12,self.mH)
        BRH4e = myCSW.HiggsBR(12,self.mH)
        BR = 0.0
        if( self.channel == self.ID_4mu ): BR = BRH4mu
        if( self.channel == self.ID_4e ): BR = BRH4e
        if( self.channel == self.ID_2e2mu ): BR = BRH2e2mu

        #HZZ Branching ratio for ZH,WH,ttH samples
        BRZZ = myCSW.HiggsBR(11,self.mH)
    
        sigEfficiency_ggH = rfvSigEff_ggH.getVal()
        sigEfficiency_qqH = rfvSigEff_qqH.getVal()
        sigEfficiency_ZH = rfvSigEff_ZH.getVal()
        sigEfficiency_WH = rfvSigEff_WH.getVal()
        sigEfficiency_ttH = rfvSigEff_ttH.getVal()

        if(DEBUG):
            print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
            print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e,", BRZZ: ",BRZZ

        
        ## SIG YIELDS
        sigRate_ggH = CS_ggH*BR*sigEfficiency_ggH*1000.*self.lumi
        sigRate_VBF = CS_VBF*BR*sigEfficiency_qqH*1000.*self.lumi
        sigRate_WH = CS_WH*BRZZ*sigEfficiency_WH*1000.*self.lumi
        sigRate_ZH = CS_ZH*BRZZ*sigEfficiency_ZH*1000.*self.lumi
        sigRate_ttH = CS_ttH*BRZZ*sigEfficiency_ttH*1000.*self.lumi

        rfvTag_Ratio_ggH = ROOT.RooFormulaVar()
        rfvTag_Ratio_qqH = ROOT.RooFormulaVar()
        rfvTag_Ratio_WH = ROOT.RooFormulaVar()
        rfvTag_Ratio_ZH = ROOT.RooFormulaVar()
        rfvTag_Ratio_ttH = ROOT.RooFormulaVar()

        if self.bVBF:
            if self.VBFcat:
                tag_Ratio_Name = "hzz4l_tag_ratio_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_ggH = ROOT.RooFormulaVar(tag_Ratio_Name,theInputs['tagged_ggH_ratio'],ROOT.RooArgList(self.MH))
                tag_Ratio_Name = "hzz4l_tag_ratio_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_qqH = ROOT.RooFormulaVar(tag_Ratio_Name,theInputs['tagged_qqH_ratio'],ROOT.RooArgList(self.MH))
                tag_Ratio_Name = "hzz4l_tag_ratio_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_WH = ROOT.RooFormulaVar(tag_Ratio_Name,theInputs['tagged_WH_ratio'],ROOT.RooArgList(self.MH))
                tag_Ratio_Name = "hzz4l_tag_ratio_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_ZH = ROOT.RooFormulaVar(tag_Ratio_Name,theInputs['tagged_ZH_ratio'],ROOT.RooArgList(self.MH))
                tag_Ratio_Name = "hzz4l_tag_ratio_ttH_{0:.0f}_{1:.0f}_{2}_a1".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_ttH = ROOT.RooFormulaVar(tag_Ratio_Name,theInputs['tagged_ttH_ratio'],ROOT.RooArgList(self.MH))

            else:
                tag_Ratio_Name = "hzz4l_tag_ratio_ggH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_ggH = ROOT.RooFormulaVar(tag_Ratio_Name,"(@1-("+theInputs['tagged_ggH_ratio']+"))",ROOT.RooArgList(self.MH,one))
                tag_Ratio_Name = "hzz4l_tag_ratio_qqH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_qqH = ROOT.RooFormulaVar(tag_Ratio_Name,"(@1-("+theInputs['tagged_qqH_ratio']+"))",ROOT.RooArgList(self.MH,one))
                tag_Ratio_Name = "hzz4l_tag_ratio_WH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_WH = ROOT.RooFormulaVar(tag_Ratio_Name,"(@1-("+theInputs['tagged_WH_ratio']+"))",ROOT.RooArgList(self.MH,one))
                tag_Ratio_Name = "hzz4l_tag_ratio_ZH_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_ZH = ROOT.RooFormulaVar(tag_Ratio_Name,"(@1-("+theInputs['tagged_ZH_ratio']+"))",ROOT.RooArgList(self.MH,one))
                tag_Ratio_Name = "hzz4l_tag_ratio_ttH_{0:.0f}_{1:.0f}_{2}_a1".format(self.channel,self.sqrts,self.VBFcat)
                rfvTag_Ratio_ttH = ROOT.RooFormulaVar(tag_Ratio_Name,"(@1-("+theInputs['tagged_ttH_ratio']+"))",ROOT.RooArgList(self.MH,one))

            print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvTag_Ratio_ggH.getVal()
            print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvTag_Ratio_qqH.getVal()
            print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvTag_Ratio_WH.getVal()
            print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvTag_Ratio_ZH.getVal()
            print "@@@@@@@@@@@@@@@@@@@@@@ ", rfvTag_Ratio_ttH.getVal()
            sigRate_ggH *= rfvTag_Ratio_ggH.getVal()
            sigRate_VBF *= rfvTag_Ratio_qqH.getVal()
            sigRate_WH *= rfvTag_Ratio_WH.getVal()
            sigRate_ZH *= rfvTag_Ratio_ZH.getVal()
            sigRate_ttH *= rfvTag_Ratio_ttH.getVal()
       
        tmpNormSigNoConv = signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigConv = sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigHM   = sig_ggH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
      
        normalizationSignal = 0.0
        if self.isHighMass : normalizationSignal = tmpNormSigHM
        else : normalizationSignal = self.getVariable(tmpNormSigNoConv,tmpNormSigConv,self.bUseCBnoConvolution)
            
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        print "#################### ",sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### norm Signal",normalizationSignal
        
        sclFactorSig_ggH = sigRate_ggH/normalizationSignal
        sclFactorSig_VBF = sigRate_VBF/normalizationSignal
        sclFactorSig_WH = sigRate_WH/normalizationSignal
        sclFactorSig_ZH = sigRate_ZH/normalizationSignal
        sclFactorSig_ttH = sigRate_ttH/normalizationSignal

        integral_ggH = 0.0
        integral_VBF = 0.0
        integral_WH  = 0.0
        integral_ZH  = 0.0
        integral_ttH = 0.0

        if self.isHighMass : integral_ggH = sig_ggH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ggH = self.getVariable(signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_VBF = sig_VBF_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_VBF = self.getVariable(signalCB_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_WH = sig_WH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_WH = self.getVariable(signalCB_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        

        if self.isHighMass : integral_ZH = sig_ZH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ZH = self.getVariable(signalCB_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_ttH = sig_ttH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ttH = self.getVariable(signalCB_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        
        sigRate_ggH_Shape = sclFactorSig_ggH*integral_ggH
        sigRate_VBF_Shape = sclFactorSig_VBF*integral_VBF
        sigRate_WH_Shape = sclFactorSig_WH*integral_WH
        sigRate_ZH_Shape = sclFactorSig_ZH*integral_ZH
        sigRate_ttH_Shape = sclFactorSig_ttH*integral_ttH
        
        if not self.bVBF:
            normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        else:
            normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}_{2}".format(self.channel,self.sqrts,self.VBFcat)
        rrvNormSig = ROOT.RooRealVar()

        

        if self.isHighMass :
            rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, sig_ggH_HM.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal())
        else :
            rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.getVariable(signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),self.bUseCBnoConvolution))
        rrvNormSig.setConstant(True)
        print "!!!%%%*** ",rrvNormSig.getVal()
        print "!!!%%%*** ",integral_ggH
        

        #rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvSigEff_ggH, rhfXsBrFuncV_1))

        rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*1000*{0}*{2}/{1}*{3}".format(self.lumi,rrvNormSig.getVal(),integral_ggH,self.getVariable(rfvTag_Ratio_ggH.getVal(),one.getVal(),self.bVBF)),ROOT.RooArgList(rfvSigEff_ggH, rhfXsBrFuncV_1))
        
        print "Compare integrals: integral_ggH=",integral_ggH,"  ; calculated=",self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)
        
        rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_norm","@0*@1*1000*{0}*{2}/{1}*{3}".format(self.lumi,rrvNormSig.getVal(),integral_VBF,self.getVariable(rfvTag_Ratio_qqH.getVal(),one.getVal(),self.bVBF)),ROOT.RooArgList(rfvSigEff_qqH, rhfXsBrFuncV_2))
        
        
        rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm","@0*@1*1000*{0}*{2}/{1}*{3}".format(self.lumi,rrvNormSig.getVal(),integral_WH,self.getVariable(rfvTag_Ratio_WH.getVal(),one.getVal(),self.bVBF)),ROOT.RooArgList(rfvSigEff_WH, rhfXsBrFuncV_3))
        
        
        rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm","@0*@1*1000*{0}*{2}/{1}*{3}".format(self.lumi,rrvNormSig.getVal(),integral_ZH,self.getVariable(rfvTag_Ratio_ZH.getVal(),one.getVal(),self.bVBF)),ROOT.RooArgList(rfvSigEff_ZH, rhfXsBrFuncV_4))
        
        
        rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm","@0*@1*1000*{0}*{2}/{1}*{3}".format(self.lumi,rrvNormSig.getVal(),integral_ttH,self.getVariable(rfvTag_Ratio_ttH.getVal(),one.getVal(),self.bVBF)),ROOT.RooArgList(rfvSigEff_ttH, rhfXsBrFuncV_5))
        

        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal()
        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal()
        if (self.all_chan):
            print "Requested to sum up over the 5 chans: the norm in rfvSigRate_ggH should be the sum of the values of sigRate_XYZ_Shape variables:"
        print " @@@@@@@ norm sig = ",rrvNormSig.getVal()
        print " @@@@@@@ rfvSigRate_ggH = ",rfvSigRate_ggH.getVal()
        print " sigRate_ggH_Shape=",sigRate_ggH_Shape
        print " @@@@@@@ rfvSigRate_VBF = ",rfvSigRate_VBF.getVal()
        print " sigRate_VBF_Shape=",sigRate_VBF_Shape
        print " @@@@@@@ rfvSigRate_WH = ",rfvSigRate_WH.getVal()
        print " sigRate_WH_Shape=",sigRate_WH_Shape
        print " @@@@@@@ rfvSigRate_ZH = ",rfvSigRate_ZH.getVal()
        print " sigRate_ZH_Shape=",sigRate_ZH_Shape
        print " @@@@@@@ rfvSigRate_ttH = ",rfvSigRate_ttH.getVal()
        print " sigRate_ttH_Shape=",sigRate_ttH_Shape
        print "Sum of sigRate_XYZ_Shape=",sigRate_ggH_Shape+sigRate_VBF_Shape+sigRate_WH_Shape+sigRate_ZH_Shape+sigRate_ttH_Shape
        ## SET RATES TO 1 
        ## DC RATES WILL BE MULTIPLIED
        ## BY RATES IMPORTED TO WS
        sigRate_ggH_Shape = 1
        sigRate_VBF_Shape = 1
        sigRate_WH_Shape = 1
        sigRate_ZH_Shape = 1
        sigRate_ttH_Shape = 1

             
        ## ----------------------- BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling
        bkgRate_qqzz = theInputs['qqZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_ggzz = theInputs['ggZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_zjets = theInputs['zjets_rate']/theInputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()

        print "channel: "+self.appendName
        print "fullrange zjets: ",normalizationBackground_zjets
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorBkg_ggzz = self.lumi*bkgRate_ggzz/normalizationBackground_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets

        CMS_zz4l_mass.setRange("tempregion2",100.,200.)
        print "100-200: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion2")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion3",100.,300.)
        print "100-300: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion3")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion4",100.,400.)
        print "100-400: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion4")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion5",100.,500.)
        print "100-500: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion5")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion6",100.,600.)
        print "100-600: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion6")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion7",100.,700.)
        print "100-700: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion7")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion8",100.,800.)
        print "100-800: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion8")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion9",100.,900.)
        print "100-900: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion9")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion10",100.,1000.)
        print "100-1000: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion10")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion11",100.,1100.)
        print "100-1100: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion11")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion12",100.,1200.)
        print "100-1200: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion12")).getVal()/normalizationBackground_zjets
        CMS_zz4l_mass.setRange("tempregion13",100.,1300.)
        print "100-1300: ", bkg_zjets.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("tempregion13")).getVal()/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_ggzz_Shape = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_zjets_Shape = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        
        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",low_M," - ",high_M
            CMS_zz4l_mass.setRange("lowmassregion",100.,160.)
            bkgRate_qqzz_lowmassregion = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_ggzz_lowmassregion = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_zjets_lowmassregion = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            lowmassyield = bkgRate_qqzz_lowmassregion + bkgRate_ggzz_lowmassregion + bkgRate_zjets_lowmassregion
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        if not self.bVBF:
            dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
        else:
            dataFileName = "{0}/hzz{1}_{2}_{3}.root".format(dataFileDir,self.appendName,self.lumi,self.VBFcat)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        
        if (self.is2D == 0):
            if(self.bIncludingError): data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass, RelErr))
            else: data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass))
		
        if (self.is2D == 1):
            if(self.bIncludingError and not self.Use3D):
                data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,self.getVariable(MEKD,D,self.bMEKD), RelErr))
            elif(self.Use3D and not self.bIncludingError):
                data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,self.getVariable(MEKD,D,self.bMEKD),self.getVariable(VD,pt,self.VBFcat)))
            elif(self.Use3D and self.bIncludingError):
                data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,self.getVariable(MEKD,D,self.bMEKD),self.getVariable(VD,pt,self.VBFcat)), RelErr)
            elif(not self.bIncludingError and not self.Use3D):
                data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,self.getVariable(MEKD,D,self.bMEKD)))

        if (self.is2D == 2):
            data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(SD))
            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""
        name_ShapeWSXSBR = ""

        if not self.bVBF:
            if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)

            if (endsInP5): name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

            if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

            name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)
        else:
            if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)

            if (endsInP5): name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)
            else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)

            if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)
            else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)

            name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV_{2}.input.root".format(self.appendName,self.sqrts,self.VBFcat)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
        w.importClassCode(RooggZZPdf_v2.Class(),True)
        w.importClassCode(RooRelBWUFParam.Class(),True)
        w.importClassCode(RooDoubleCB.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)
        if self.isHighMass :
            w.importClassCode(RooRelBWHighMass.Class(),True)

        if( FactorizedShapes ):
            if( self.channel == self.ID_4mu ):
                w.importClassCode(RooFourMuMassShapePdf2.Class(),True)
                w.importClassCode(RooFourMuMassRes.Class(),True)
            elif( self.channel == self.ID_4e ):
                w.importClassCode(RooFourEMassShapePdf2.Class(),True)
                w.importClassCode(RooFourEMassRes.Class(),True)
            elif( self.channel == self.ID_2e2mu ):
                w.importClassCode(RooTwoETwoMuMassShapePdf2.Class(),True)
                w.importClassCode(RooTwoETwoMuMassRes.Class(),True)
            
                
                
        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
    
        if(self.bUseCBnoConvolution) :
            if (self.is2D == 0):
		if not self.bIncludingError:
                	signalCB_ggH.SetNameTitle("ggH","ggH")
                	signalCB_VBF.SetNameTitle("qqH","qqH")
                	signalCB_WH.SetNameTitle("WH","WH")
                	signalCB_ZH.SetNameTitle("ZH","ZH")
                	signalCB_ttH.SetNameTitle("ttH","ttH")
                
                	getattr(w,'import')(signalCB_ggH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_VBF, ROOT.RooFit.RecycleConflictNodes())
               		getattr(w,'import')(signalCB_WH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_ZH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_ttH, ROOT.RooFit.RecycleConflictNodes())
		else:
                	sig_ggHErr.SetNameTitle("ggH","ggH")
                	sig_VBFErr.SetNameTitle("qqH","qqH")
                	sig_WHErr.SetNameTitle("WH","WH")
                	sig_ZHErr.SetNameTitle("ZH","ZH")
                	sig_ttHErr.SetNameTitle("ttH","ttH")
                
                	getattr(w,'import')(sig_ggHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_VBFErr, ROOT.RooFit.RecycleConflictNodes())
               		getattr(w,'import')(sig_WHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_ZHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_ttHErr, ROOT.RooFit.RecycleConflictNodes())
                
            if (self.is2D == 1):
                if not self.Use3D:
                    sigCB2d_ggH.SetNameTitle("ggH","ggH")
                    sigCB2d_VBF.SetNameTitle("qqH","qqH")
                    sigCB2d_WH.SetNameTitle("WH","WH")
                    sigCB2d_ZH.SetNameTitle("ZH","ZH")
                    sigCB2d_ttH.SetNameTitle("ttH","ttH")
                
                    getattr(w,'import')(sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_VBF, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_WH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_ZH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_ttH, ROOT.RooFit.RecycleConflictNodes())
                    if(self.isAltSig):
                        sigCB2d_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                        getattr(w,'import')(sigCB2d_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())

                else:
                    sigCB2d_ggH_VBF_KD.SetNameTitle("ggH","ggH")
                    sigCB2d_qqH_VBF_KD.SetNameTitle("qqH","qqH")
                    sigCB2d_WH_VBF_KD.SetNameTitle("WH","WH")
                    sigCB2d_ZH_VBF_KD.SetNameTitle("ZH","ZH")
                    sigCB2d_ttH_VBF_KD.SetNameTitle("ttH","ttH")
                    getattr(w,'import')(sigCB2d_ggH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_qqH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_WH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_ZH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_ttH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())

            if (self.is2D == 2):
                sigTemplateSDPdf_ggH.SetNameTitle("ggH","ggH")
                sigTemplateSDPdf_VBF.SetNameTitle("qqH","qqH")
                sigTemplateSDPdf_WH.SetNameTitle("WH","WH")
                sigTemplateSDPdf_ZH.SetNameTitle("ZH","ZH")
                sigTemplateSDPdf_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sigTemplateSDPdf_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ttH, ROOT.RooFit.RecycleConflictNodes())

        else:
                
            if (self.is2D == 0):

                if self.isHighMass:
                    sig_ggH_HM.SetNameTitle("ggH","ggH")
                    sig_VBF_HM.SetNameTitle("qqH","qqH")
                    sig_WH_HM.SetNameTitle("WH","WH")
                    sig_ZH_HM.SetNameTitle("ZH","ZH")
                    sig_ttH_HM.SetNameTitle("ttH","ttH")
                    
                    getattr(w,'import')(sig_ggH_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_VBF_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_WH_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ZH_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ttH_HM, ROOT.RooFit.RecycleConflictNodes())

                else :
                    sig_ggH.SetNameTitle("ggH","ggH")
                    sig_VBF.SetNameTitle("qqH","qqH")
                    sig_WH.SetNameTitle("WH","WH")
                    sig_ZH.SetNameTitle("ZH","ZH")
                    sig_ttH.SetNameTitle("ttH","ttH")
                    
                    getattr(w,'import')(sig_ggH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_VBF, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_WH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ZH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ttH, ROOT.RooFit.RecycleConflictNodes())
                    
            if (self.is2D == 1):
                if not self.Use3D:
                    sig2d_ggH.SetNameTitle("ggH","ggH")
                    sig2d_VBF.SetNameTitle("qqH","qqH")
                    sig2d_WH.SetNameTitle("WH","WH")
                    sig2d_ZH.SetNameTitle("ZH","ZH")
                    sig2d_ttH.SetNameTitle("ttH","ttH")
                
                    getattr(w,'import')(sig2d_ggH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_VBF, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_WH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_ZH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_ttH, ROOT.RooFit.RecycleConflictNodes())
                    if(self.isAltSig):
                        sigCB2d_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                        getattr(w,'import')(sigCB2d_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())

                else:
                    sig2d_ggH_VBF_KD.SetNameTitle("ggH","ggH")
                    sig2d_qqH_VBF_KD.SetNameTitle("qqH","qqH")
                    sig2d_WH_VBF_KD.SetNameTitle("WH","WH")
                    sig2d_ZH_VBF_KD.SetNameTitle("ZH","ZH")
                    sig2d_ttH_VBF_KD.SetNameTitle("ttH","ttH")
                    getattr(w,'import')(sig2d_ggH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_qqH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_WH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_ZH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig2d_ttH_VBF_KD,ROOT.RooFit.RecycleConflictNodes())

            if (self.is2D == 2): 
                sigTemplateSDPdf_ggH.SetNameTitle("ggH","ggH")
                sigTemplateSDPdf_VBF.SetNameTitle("qqH","qqH")
                sigTemplateSDPdf_WH.SetNameTitle("WH","WH")
                sigTemplateSDPdf_ZH.SetNameTitle("ZH","ZH")
                sigTemplateSDPdf_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sigTemplateSDPdf_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ttH, ROOT.RooFit.RecycleConflictNodes())

 
        if (self.is2D == 0):
		if not self.bIncludingError:
			bkg_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
			bkg_ggzz.SetNameTitle("bkg_ggzz","bkg_ggzz")
			bkg_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
            		getattr(w,'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_ggzz, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
		else:
			bkg_qqzzErr.SetNameTitle("bkg_qqzz","bkg_qqzz")
			bkg_ggzzErr.SetNameTitle("bkg_ggzz","bkg_ggzz")
			bkg_zjetsErr.SetNameTitle("bkg_zjets","bkg_zjets")
            		getattr(w,'import')(bkg_qqzzErr, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_ggzzErr, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_zjetsErr, ROOT.RooFit.RecycleConflictNodes())
            
        if (self.is2D == 1):
            if not self.Use3D:
                getattr(w,'import')(bkg2d_qqzz,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(bkg2d_ggzz,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(bkg2d_zjets,ROOT.RooFit.RecycleConflictNodes())
            else:
                bkg2d_qqZZ_VBF_KD.SetNameTitle("bkg2d_qqzz","bkg2d_qqzz")
                bkg2d_ggZZ_VBF_KD.SetNameTitle("bkg2d_ggzz","bkg2d_ggzz")
                bkg2d_ZX_VBF_KD.SetNameTitle("bkg2d_zjets","bkg2d_zjets")
                getattr(w,'import')(bkg2d_qqZZ_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(bkg2d_ggZZ_VBF_KD,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(bkg2d_ZX_VBF_KD,ROOT.RooFit.RecycleConflictNodes())

        if (self.is2D == 2): 
            bkgTemplateSDPdf_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
            bkgTemplateSDPdf_ggzz.SetNameTitle("bkg_ggzz","bkg_ggzz")
            bkgTemplateSDPdf_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
            getattr(w,'import')(bkgTemplateSDPdf_ggzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkgTemplateSDPdf_qqzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkgTemplateSDPdf_zjets, ROOT.RooFit.RecycleConflictNodes())

        
        getattr(w,'import')(rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_VBF, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_WH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ZH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ttH, ROOT.RooFit.RecycleConflictNodes())
        if(self.isAltSig):
            rfvSigRate_ggH_ALT=ROOT.RooFormulaVar(rfvSigRate_ggH,"ggH{0}_norm".format(self.appendHypType))
            print 'Compare signal rates: STD=',rfvSigRate_ggH.getVal(),"   ALT=",rfvSigRate_ggH_ALT.getVal()
            getattr(w,'import')(rfvSigRate_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
            
        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)
        
        ## --------------------------- DATACARDS -------------------------- ##

        systematics.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape, bkgRate_zjets_Shape)
        systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape,bkgRate_zjets_Shape)

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan and not self.all_chan :  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = bkgRate_qqzz_Shape
        rates['ggZZ']  = bkgRate_ggzz_Shape
        rates['zjets'] = bkgRate_zjets_Shape
        rates['ttbar'] = 0
        rates['zbb']   = 0
        

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        if not self.bVBF:
            systematics.WriteSystematics(fo, theInputs)
            systematics.WriteShapeSystematics(fo,theInputs)
        else:
            systematics.WriteSystematics(fo, theInputs,self.VBFcat,self.Use3D)
            systematics.WriteShapeSystematics(fo,theInputs)
        fo.close()

        if(self.isAltSig):
            if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
            fo = open( name_Shape, "wb")
            self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D,True,self.appendHypType )
            systematics.WriteSystematics(fo, theInputs)
            systematics.WriteShapeSystematics(fo,theInputs)
            fo.close()

        ## forXSxBR
        if not self.bVBF:
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
        else:
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.bVBF)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.bVBF)
            
        fo = open( name_Shape, "wb" )
        if not self.bVBF:
            self.WriteDatacard(fo, theInputs,name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        else:
            self.WriteDatacard(fo,theInputs,name_ShapeWS2,rates,data_obs.numEntries(),self.is2D,False,"",True,self.VBFcat)
        systematics_forXSxBR.WriteSystematics(fo, theInputs)
        systematics_forXSxBR.WriteShapeSystematics(fo,theInputs)
        fo.close()

        if(self.isAltSig):
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.appendHypType)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.appendHypType)
            fo = open( name_Shape, "wb")
            self.WriteDatacard(fo,theInputs,name_ShapeWS2,rates,data_obs.numEntries(),self.is2D,True,self.appendHypType )
            systematics.WriteSystematics(fo, theInputs)
            systematics.WriteShapeSystematics(fo,theInputs)
            fo.close()
        


    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents,is2D,isAltCard=False,AltLabel="",bVBF=False,VBFcat=""):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)
        
        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        
        if not self.bVBF:
            file.write("bin a{0} \n".format(self.channel))
        else:
            file.write("bin a{0}_{1} \n".format(self.channel,self.VBFcat))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("bin ")        

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']

        channelName1D=['ggH','qqH','WH','ZH','ttH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
        channelName2D=['ggH','qqH','WH','ZH','ttH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

#            channelList=['ggH{0}'.format(AltLabel),'qqZZ','ggZZ','zjets','ttbar','zbb']

        if theInputs["all"]:
            channelList=['ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
            if isAltCard :
                channelName2D=['ggH{0}'.format(AltLabel),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
            else:
                channelName2D=['ggH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
         
        for chan in channelList:
            if theInputs[chan]:
                if not self.bVBF:
                    file.write("a{0} ".format(self.channel))
                else:
                    file.write("a{0}_{1} ".format(self.channel,self.VBFcat))
            else:
                if chan.startswith("ggH") and theInputs["all"] :
                    file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0
        if not (self.is2D == 1):
            for chan in channelList:
                if theInputs[chan]:
                    file.write("{0} ".format(channelName1D[i]))
                i+=1
        else:
            for chan in channelList:
#                print 'checking if ',chan,' is in the list of to-do'
                if theInputs[chan]:
                    file.write("{0} ".format(channelName2D[i]))
 #                   print 'writing in card index=',i,'  chan=',chan
                    i+=1
                else:
                    if chan.startswith("ggH") and theInputs["all"] :
                        file.write("{0} ".format(channelName2D[i]))
  #                      print 'writing in card TOTAL SUM, index=',i,'  chan=',chan,'  ',channelName2D[i]
                        i+=1
        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        for chan in channelList:
            if theInputs[chan] or (chan.startswith("ggH") and theInputs["all"]):
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


        
    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggH']: counter+=1
        if inputs['qqH']: counter+=1
        if inputs['WH']:  counter+=1
        if inputs['ZH']:  counter+=1
        if inputs['ttH']: counter+=1
        if inputs['all']: counter+=1
        
        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1
        
        return counter

