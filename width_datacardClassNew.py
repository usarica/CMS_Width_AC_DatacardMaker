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

class width_datacardClass:

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
 
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
    
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theLowSide, theOutputDir, theInputs, theTemplateDir="templates2D"):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = 125.6   ## FIXED
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.outputDir = theOutputDir
        self.templateDir = theTemplateDir

        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.ggZZ_signal_chan = theInputs['ggZZ_signal']
        self.ggZZ_bkg_chan = theInputs['ggZZ_bkg']
        self.ggZZ_interf_chan = theInputs['ggZZ_interf']
        self.zjets_chan = theInputs['zjets']
        
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
        
        self.low_M = theLowSide
        self.high_M = 1400
        
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
            
            
        
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)
        systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        bins = 260
        bins2 = (self.high_M-self.low_M)/5
        # if(self.bUseCBnoConvolution): bins = 200

        CMS_zz4l_mass_name = "CMS_zz4l_mass"
            
        CMS_zz4l_mass = ROOT.RooRealVar(CMS_zz4l_mass_name,CMS_zz4l_mass_name,self.low_M,self.high_M)
        CMS_zz4l_mass.setBins(bins2)

        ## use this variable only For Integration (FI)
        CMS_zz4l_mass_name = "CMS_zz4l_mass_FI"
            
        CMS_zz4l_mass_FI = ROOT.RooRealVar(CMS_zz4l_mass_name,CMS_zz4l_mass_name,100,1400)    
        CMS_zz4l_mass_FI.setBins(bins)

        x_name = "CMS_zz4l_GGsm"

        x = ROOT.RooRealVar(x_name,x_name,0.0001,100)
        x.setVal(1)
        x.setBins(100)

        mu_name = "CMS_zz4l_mu"

        mu = ROOT.RooRealVar(mu_name,mu_name,0.1,10)
        mu.setVal(1)
        mu.setBins(100)

        D2Name = "CMS_zz4l_widthKD"

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",300.)   
        self.MH.setConstant(True)
        
        print '1D signal shapes for Width'
      ##  templateSigName = "{0}/prova_{1}_{2:.0f}TeV_m4l.root".format(self.templateDir, self.appendName, self.sqrts)
        templateSigName = "{0}/templ1D_{1}_{2:.0f}TeV_m4l.root".format(self.templateDir, self.appendName, self.sqrts)
        sigTempFile = ROOT.TFile(templateSigName)
        print templateSigName
        Sig_T_1 = sigTempFile.Get("mZZ_bkg")
        Sig_T_2 = sigTempFile.Get("mZZ_sig")
        Sig_T_4 = sigTempFile.Get("mZZ_inter")

        dBinsX = Sig_T_1.GetXaxis().GetNbins()
        print "X bins: ",dBinsX
        dLowX = Sig_T_1.GetXaxis().GetXmin()
        dHighX = Sig_T_1.GetXaxis().GetXmax()
        wBinsX = Sig_T_1.GetXaxis().GetBinWidth(1)
        
        dBinsY = Sig_T_1.GetYaxis().GetNbins()
        print "Y bins: ",dBinsY
        dLowY = Sig_T_1.GetYaxis().GetXmin()
        dHighY = Sig_T_1.GetYaxis().GetXmax()
        
        D2 = ROOT.RooRealVar(D2Name,D2Name,dLowY,dHighY)
        CMS_zz4l_mass.setBins(dBinsX)
        CMS_zz4l_mass_FI.setBins(dBinsX)
        D2.setBins(dBinsY)

        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)

        ## -------------------------- gg2ZZ SHAPES ---------------------------------- ##
        
        sigRateName = "signal_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigRates = ROOT.RooRealVar(sigRateName,"sigRates",0.0,10000.0)
        bkgRateName = "bkg_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgRates = ROOT.RooRealVar(bkgRateName,"bkgRates",0.0,10000.0)
        interfRateName = "interf_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        interfRates = ROOT.RooRealVar(interfRateName,"interfRates",0.0,10000.0)
        sigRateNameNorm = "signalNorm_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigRatesNorm = ROOT.RooFormulaVar(sigRateNameNorm,"@0*@1/(@0*@1-sqrt(@0*@1)+1)",ROOT.RooArgList(x,mu))
        interfRateNameNorm = "interfNorm_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        interfRatesNorm = ROOT.RooFormulaVar(interfRateNameNorm,"-sqrt(@0*@1)/(@0*@1-sqrt(@0*@1)+1)",ROOT.RooArgList(x,mu))

        #ggZZpdfName = "ggZZ_RooWidth_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        #ggZZpdf = ROOT.HZZ4lWidth(ggZZpdfName,ggZZpdfName,CMS_zz4l_mass,one,x,bkgRates,sigRates,interfRates,Sig_T_1,Sig_T_2,Sig_T_4)

        TemplateName = "ggZZsignal_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZsignal_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass),Sig_T_2)
        PdfName = "ggZZsignal_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZsignal_TemplatePdf = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_mass),ggZZsignal_TempDataHist)

        TemplateName = "ggZZbkg_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZbkg_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass),Sig_T_1)
        PdfName = "ggZZbkg_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZbkg_TemplatePdf = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_mass),ggZZbkg_TempDataHist)
        
        TemplateName = "ggZZinterf_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZinterf_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass),Sig_T_4)
        PdfName = "ggZZinterf_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZinterf_TemplatePdf = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_mass),ggZZinterf_TempDataHist)

        ggZZpdfName = "ggZZ_RooWidth_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZpdf = ROOT.RooRealSumPdf(ggZZpdfName,ggZZpdfName,ROOT.RooArgList(ggZZsignal_TemplatePdf,ggZZinterf_TemplatePdf,ggZZbkg_TemplatePdf),ROOT.RooArgList(sigRatesNorm,interfRatesNorm))

        ## for integration
        TemplateName = "ggZZsignal_TempDataHist_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
        ggZZsignal_TempDataHist_FI = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass_FI),Sig_T_2)
        PdfName = "ggZZsignal_TemplatePdf_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
        ggZZsignal_TemplatePdf_FI = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_mass_FI),ggZZsignal_TempDataHist_FI)
      ##   PDFName = "ggZZsignal_TemplatePDF_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
      ##   ggZZsignal_TemplatePDF_FI = ROOT.RooHistPdf(PDFName,PDFName,ROOT.RooArgSet(CMS_zz4l_mass_FI),ggZZsignal_TempDataHist_FI)

        TemplateName = "ggZZbkg_TempDataHist_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
        ggZZbkg_TempDataHist_FI = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass_FI),Sig_T_1)
        PdfName = "ggZZbkg_TemplatePdf_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
        ggZZbkg_TemplatePdf_FI = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_mass_FI),ggZZbkg_TempDataHist_FI)
        
        TemplateName = "ggZZinterf_TempDataHist_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
        ggZZinterf_TempDataHist_FI = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass_FI),Sig_T_4)
        PdfName = "ggZZinterf_TemplatePdf_{0:.0f}_{1:.0f}_FI".format(self.channel,self.sqrts)
        ggZZinterf_TemplatePdf_FI = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_mass_FI),ggZZinterf_TempDataHist_FI)

        ggZZpdfName = "ggZZ_RooWidth_FI_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZpdf_FI = ROOT.RooRealSumPdf(ggZZpdfName,ggZZpdfName,ROOT.RooArgList(ggZZsignal_TemplatePdf_FI,ggZZinterf_TemplatePdf_FI,ggZZbkg_TemplatePdf_FI),ROOT.RooArgList(sigRatesNorm,interfRatesNorm))

        ## -------------------------- OTHER BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution
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

        ## for integration
        bkg_qqzz_FI = ROOT.RooqqZZPdf_v2("bkg_qqzzTmp_FI","bkg_qqzzTmp_FI",CMS_zz4l_mass_FI,CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13)
            
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


        if (self.channel == self.ID_4mu):
            name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL_2P2F)
            name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL_2P2F)
            print "mean 4mu: ",mlZjet.getVal()
            print "sigma 4mu: ",slZjet.getVal()
            bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",CMS_zz4l_mass,mlZjet,slZjet)
            bkg_zjets_FI = ROOT.RooLandau("bkg_zjetsTmp_FI","bkg_zjetsTmp_FI",CMS_zz4l_mass_FI,mlZjet,slZjet)
        elif (self.channel == self.ID_4e):

            summa = val_normL_2P2F + val_normL_3P1F
            
            name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
            name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
            name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F/summa)
            name = "p0Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",val_pol0_2P2F)
            name = "p1Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",val_pol1_2P2F)
            print "mean 2p2f 4e: ",mlZjet_2p2f.getVal()
            print "sigma 2p2f 4e: ",slZjet_2p2f.getVal()
            print "norm 2p2f 4e: ",nlZjet_2p2f.getVal()
            print "pol0 2p2f 4e: ",p0Zjet_2p2f.getVal()
            print "pol1 2p2f 4e: ",p1Zjet_2p2f.getVal()
            bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",RooArgList(CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))

            bkg_zjets_2p2f_FI = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f_FI","bkg_zjetsTmp_2p2f_FI","(TMath::Landau(@0,@1,@2))*(1.+ TMath::Exp(@3+@4*@0))",RooArgList(CMS_zz4l_mass_FI,mlZjet_2p2f,slZjet_2p2f,p0Zjet_2p2f,p1Zjet_2p2f))
            
            name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
            name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
            name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F/summa)
            print "mean 3p1f 4e: ",mlZjet_3p1f.getVal()
            print "sigma 3p1f 4e: ",slZjet_3p1f.getVal()
            print "norm 3p1f 4e: ",nlZjet_3p1f.getVal()
            bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)

            bkg_zjets_3p1f_FI = ROOT.RooLandau("bkg_zjetsTmp_3p1f_FI","bkg_zjetsTmp_3p1f_FI",CMS_zz4l_mass_FI,mlZjet_3p1f,slZjet_3p1f)
            
            bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f))
            bkg_zjets_FI = ROOT.RooAddPdf("bkg_zjetsTmp_FI","bkg_zjetsTmp_FI",ROOT.RooArgList(bkg_zjets_2p2f_FI,bkg_zjets_3p1f_FI),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f))
            
        elif (self.channel == self.ID_2e2mu):

            summa = val_normL_2P2F + val_normL_2P2F_2 + val_normL_3P1F
            
            name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",val_meanL_2P2F)
            name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",val_sigmaL_2P2F)
            name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",val_normL_2P2F/summa)
            print "mean 2p2f 2mu2e: ",mlZjet_2p2f.getVal()
            print "sigma 2p2f 2mu2e: ",slZjet_2p2f.getVal()
            print "norm 2p2f 2mu2e: ",nlZjet_2p2f.getVal()
            bkg_zjets_2p2f = ROOT.RooLandau("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f",CMS_zz4l_mass,mlZjet_2p2f,slZjet_2p2f)
            
            bkg_zjets_2p2f_FI = ROOT.RooLandau("bkg_zjetsTmp_2p2f_FI","bkg_zjetsTmp_2p2f_FI",CMS_zz4l_mass_FI,mlZjet_2p2f,slZjet_2p2f)
            
            name = "mlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_2p2f_2 = ROOT.RooRealVar(name,"mean landau Zjet 2p2f 2e2mu",val_meanL_2P2F_2)
            name = "slZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_2p2f_2 = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f 2e2mu",val_sigmaL_2P2F_2)
            name = "nlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_2p2f_2 = ROOT.RooRealVar(name,"norm landau Zjet 2p2f 2e2mu",val_normL_2P2F_2/summa)
            print "mean 2p2f 2e2mu: ",mlZjet_2p2f_2.getVal()
            print "sigma 2p2f 2e2mu: ",slZjet_2p2f_2.getVal()
            print "norm 2p2f 2e2mu: ",nlZjet_2p2f_2.getVal()
            bkg_zjets_2p2f_2 = ROOT.RooLandau("bkg_zjetsTmp_2p2f_2","bkg_zjetsTmp_2p2f_2",CMS_zz4l_mass,mlZjet_2p2f_2,slZjet_2p2f_2)
            
            bkg_zjets_2p2f_2_FI = ROOT.RooLandau("bkg_zjetsTmp_2p2f_2_FI","bkg_zjetsTmp_2p2f_2_FI",CMS_zz4l_mass_FI,mlZjet_2p2f_2,slZjet_2p2f_2)
            
            name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",val_meanL_3P1F)
            name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",val_sigmaL_3P1F)
            name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",val_normL_3P1F/summa)
            print "mean 3p1f 2mu2e: ",mlZjet_3p1f.getVal()
            print "sigma 3p1f 2mu2e: ",slZjet_3p1f.getVal()
            print "norm 3p1f 2mu2e: ",nlZjet_3p1f.getVal()
            bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",CMS_zz4l_mass,mlZjet_3p1f,slZjet_3p1f)
       
            bkg_zjets_3p1f_FI = ROOT.RooLandau("bkg_zjetsTmp_3p1f_FI","bkg_zjetsTmp_3p1f_FI",CMS_zz4l_mass_FI,mlZjet_3p1f,slZjet_3p1f)
            
            bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(bkg_zjets_2p2f,bkg_zjets_3p1f,bkg_zjets_2p2f_2),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f,nlZjet_2p2f_2))

            bkg_zjets_FI = ROOT.RooAddPdf("bkg_zjetsTmp_FI","bkg_zjetsTmp_FI",ROOT.RooArgList(bkg_zjets_2p2f_FI,bkg_zjets_3p1f_FI,bkg_zjets_2p2f_2_FI),ROOT.RooArgList(nlZjet_2p2f,nlZjet_3p1f,nlZjet_2p2f_2))

        
        ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- RANGES ----------------------- ##
        
        CMS_zz4l_mass_FI.setRange("shape",self.low_M,self.high_M)            

        CMS_zz4l_mass_FI.setRange("fullrangesignal",100,1400)
        CMS_zz4l_mass_FI.setRange("fullrange",100,1400)

             
        ## ----------------------- SIGNAL AND BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling
        bkgRate_qqzz = theInputs['qqZZ_rate']/theInputs['qqZZ_lumi']
        totalRate_ggzz = theInputs['ggZZ_rate']/theInputs['qqZZ_lumi']
        ## rate_signal_ggzz = theInputs['ggZZ_signal_rate']/theInputs['qqZZ_lumi']
        ## print " @@@@@@ initial rate: ",theInputs['ggZZ_signal_rate']
        ## print " @@@@@@ corrected rate: ",rate_signal_ggzz
        ## rate_bkg_ggzz = theInputs['ggZZ_bkg_rate']/theInputs['qqZZ_lumi']
        ## rate_interf_ggzz = theInputs['ggZZ_interf_rate']/theInputs['qqZZ_lumi']
        bkgRate_zjets = theInputs['zjets_rate']/theInputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = bkg_qqzz_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("fullrange") ).getVal()
        print 
        normalizationBackground_ggzz = ggZZpdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("fullrange") ).getVal()
        #normalizationsignal_ggzz = ggZZsignal_TemplatePdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("fullrange") ).getVal()
        print " @@@@@@ total normalization: ",normalizationBackground_ggzz
        #normalizationbkg_ggzz = ggZZbkg_TemplatePdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("fullrange") ).getVal()
        #normalizationinterf_ggzz = ggZZinterf_TemplatePdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = bkg_zjets_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("fullrange") ).getVal()

        print " @@@@@@ channel: "+self.appendName
        print " @@@@@@ fullrange zz: ",normalizationBackground_qqzz
        print " @@@@@@ fullrange zjets: ",normalizationBackground_zjets
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorTotal_ggzz = self.lumi*totalRate_ggzz/normalizationBackground_ggzz
        #sclFactor_signal_ggzz = self.lumi*rate_signal_ggzz/normalizationsignal_ggzz
        print " @@@@@@ scale factor: ",sclFactorTotal_ggzz
        #sclFactor_bkg_ggzz = self.lumi*rate_bkg_ggzz/normalizationbkg_ggzz
        #sclFactor_interf_ggzz = self.lumi*rate_interf_ggzz/normalizationinterf_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * bkg_qqzz_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("shape") ).getVal()
        #totalRate_ggzz_Shape = sclFactorTotal_ggzz * ggZZpdfFI.createIntegral( ROOT.RooArgSet(CMS_zz4l_massFI), ROOT.RooFit.Range("shape") ).getVal()
        rate_signal_ggzz_Shape = sclFactorTotal_ggzz * ggZZsignal_TemplatePdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("shape") ).getVal() 
        rate_bkg_ggzz_Shape = sclFactorTotal_ggzz * ggZZbkg_TemplatePdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("shape") ).getVal() 
        rate_interf_ggzz_Shape = sclFactorTotal_ggzz * ggZZinterf_TemplatePdf_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("shape") ).getVal() 
        totalRate_ggzz_Shape = rate_signal_ggzz_Shape + rate_bkg_ggzz_Shape - rate_interf_ggzz_Shape;
        bkgRate_zjets_Shape = sclFactorBkg_zjets * bkg_zjets_FI.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass_FI), ROOT.RooFit.Range("shape") ).getVal()
        print " @@@@@@ signal normalization in signal region: ",rate_signal_ggzz_Shape
       # print " @@@@@@ signal normalization in signal region TEST: ",rate_signal_ggzz_Test
        print " @@@@@@ interf normalization in signal region: ",rate_interf_ggzz_Shape
        print " @@@@@@ bkg normalization in signal region: ",rate_bkg_ggzz_Shape

        sigRates.setVal(rate_signal_ggzz_Shape)
        sigRates.setConstant(True)
        bkgRates.setVal(rate_bkg_ggzz_Shape)
        bkgRates.setConstant(True)
        interfRates.setVal(rate_interf_ggzz_Shape)
        interfRates.setConstant(True)
        
        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",low_M," - ",high_M
            CMS_zz4l_mass.setRange("shapiro",100.,160.)
            bkgRate_qqzz_shapiro = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shapiro") ).getVal()
            totalRate_ggzz_shapiro = sclFactorTotal_ggzz * ggZZpdf.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shapiro") ).getVal()
            bkgRate_zjets_shapiro = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shapiro") ).getVal()
            lowmassyield = bkgRate_qqzz_shapiro + totalRate_ggzz_shapiro + bkgRate_zjets_shapiro
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
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
        

        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass))
        data_obs_red = data_obs.reduce("CMS_zz4l_mass > {0}".format(self.low_M))
  #      data_obs_red.append(data_obs_red)
  #      data_obs_red.append(data_obs_red)  ## 4 times the data!

            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.low_M
        if ( math.fabs(math.floor(tmpMH)-self.low_M) > 0.001): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""
        name_ShapeWSXSBR = ""

        if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        
        if (endsInP5): name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        
        if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        
        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
##        w.importClassCode(RooggZZPdf_v2.Class(),True)
##         w.importClassCode(RooRelBWUFParam.Class(),True)
##         w.importClassCode(RooDoubleCB.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)
##         if self.isHighMass :
##             w.importClassCode(RooRelBWHighMass.Class(),True)
            
                
                
        getattr(w,'import')(data_obs_red,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?

        ggZZpdf.SetNameTitle("ggzz","ggzz")
        getattr(w,'import')(ggZZpdf, ROOT.RooFit.RecycleConflictNodes())
        ggZZpdfNormName = "ggZZ_RooWidth_{0:.0f}_{1:.0f}_norm".format(self.channel,self.sqrts)
        ggZZpdf_norm = ROOT.RooFormulaVar(ggZZpdfNormName,"@0*@3*@4-@1*sqrt(@3*@4)+@2",ROOT.RooArgList(sigRates,interfRates,bkgRates,x,mu))
        ggZZpdf_norm.SetNameTitle("ggzz_norm","ggzz_norm")
        getattr(w,'import')(ggZZpdf_norm, ROOT.RooFit.RecycleConflictNodes())
       
        bkg_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
        ##ggZZsignal_TemplatePdf.SetNameTitle("ggsignalzz","ggsignalzz")
        ##ggZZbkg_TemplatePdf.SetNameTitle("ggbkgzz","ggbkgzz")
        ##ggZZinterf_TemplatePdf.SetNameTitle("gginterfzz","gginterfzz")
        bkg_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
        getattr(w,'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())
        ##getattr(w,'import')(ggZZsignal_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        ##getattr(w,'import')(ggZZbkg_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        ##getattr(w,'import')(ggZZinterf_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
  
        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)
        
        ## --------------------------- DATACARDS -------------------------- ##

        systematics.setSystematics(bkgRate_qqzz_Shape, totalRate_ggzz_Shape, bkgRate_zjets_Shape)
        systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape, totalRate_ggzz_Shape,bkgRate_zjets_Shape)

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan:  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  totalRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = bkgRate_qqzz_Shape
        rates['ggZZ']  = 1
        rates['ggZZ_signal']  = rate_signal_ggzz_Shape
        rates['ggZZ_bkg']  = rate_bkg_ggzz_Shape
        rates['ggZZ_interf']  = rate_interf_ggzz_Shape 
        rates['zjets'] = bkgRate_zjets_Shape
        rates['ttbar'] = totalRate_ggzz_Shape
        rates['zbb']   = 0
        

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs_red.numEntries())
        
        systematics.WriteSystematics(fo, theInputs)
        systematics.WriteShapeSystematics(fo,theInputs)

        print "GO THERE"
        fo.close()
        print "GOT HERE"
        

        ## forXSxBR

        if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.low_M,self.sqrts)	
        else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.low_M,self.sqrts)
            
        fo = open( name_Shape, "wb" )

        self.WriteDatacard(fo, theInputs,name_ShapeWS2, rates, data_obs.numEntries())
        
        systematics_forXSxBR.WriteSystematics(fo, theInputs)
        systematics_forXSxBR.WriteShapeSystematics(fo,theInputs)
        fo.close()
        


    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)
        
        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        

        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("## signal,bkg,interf,tot rates [{0:.4f}, {1:.4f}, -{2:.4f}, {3:.4f}] \n".format(theRates["ggZZ_signal"],theRates["ggZZ_bkg"],theRates["ggZZ_interf"],theRates["ttbar"]))
        file.write("bin ")        

        #channelList=['ggZZ_signal','ggZZ_interf','ggZZ_bkg','qqZZ','zjets']
        channelList=['ggZZ','qqZZ','zjets'] 

        #channelName=['ggsignalzz','gginterfzz','ggbkgzz','bkg_qqzz','bkg_zjets']
        channelName=['ggzz','bkg_qqzz','bkg_zjets'] 
         
        for chan in channelList:
            if theInputs[chan]:
                file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0

        for chan in channelList:
            #print 'checking if ',chan,' is in the list of to-do'
            #print "{0} ".format(channelName[i])
            if theInputs[chan]:
                file.write("{0} ".format(channelName[i]))
                #print 'writing in card index=',i,'  chan=',chan
                #print "{0} ".format(channelName[i])
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
            if theInputs[chan]:
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


        
    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggZZ']:  counter+=1
        if inputs['ggZZ_signal']: counter+=1
        if inputs['ggZZ_interf']: counter+=1
  ##       if inputs['qqH']: counter+=1
##         if inputs['WH']:  counter+=1
##         if inputs['ZH']:  counter+=1
##         if inputs['ttH']: counter+=1
        
        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        ##if inputs['ggZZ']:  counter+=1
        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ_bkg']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1
        
        return counter

