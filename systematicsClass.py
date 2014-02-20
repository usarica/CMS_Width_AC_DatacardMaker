#! /usr/bin/env python
import os
import re
import math
from ROOT import *
import ROOT
from array import array


## ------------------------------------
##  systematics class
## ------------------------------------

class systematicsClass:

    def __init__(self,theMass,theForXSxBR,theisFSR,theInputs):

        self.ID_4mu = 1
        self.ID_4e = 2
        self.ID_2e2mu = 3

        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.mH = theMass
        self.isForXSxBR = theForXSxBR
        self.isFSR = theisFSR
        self.model = theInputs['model']


        self.muSelError = 0.0
        self.eSelError = 0.0
        self.muSelError2e2mu = 0.0
        self.eSelError2e2mu = 0.0
        self.muSelErrorZZ2e2mu = 0.0
        self.eSelErrorZZ2e2mu = 0.0

        self.qqVV_scaleSys = 0.0
        self.qqVV_pdfSys = 0.0
        self.ggVV_scaleSys = 0.0
        self.ggVV_pdfSys = 0.0

        self.lumiUncertainty = theInputs['lumiUnc']
        self.sel_muontrig = theInputs['muonTrigUnc']
        self.sel_muonfull = theInputs['muonFullUnc']

        self.sel_eletrig = theInputs['elecTrigUnc']
        self.sel_elefull = theInputs['elecFullUnc']

        self.zjetKappaLow = theInputs['zjetsKappaLow']
        self.zjetKappaHigh = theInputs['zjetsKappaHigh']

        self.QCD_scale_ggH_2j_sys = theInputs['QCD_scale_ggH_2j_sys']
        self.QCD_scale_qqH_2j_sys = theInputs['QCD_scale_qqH_2j_sys']
        #self.QCD_scale_qqZZ_2j_sys = theInputs['QCD_scale_qqZZ_2j_sys']

        self.theoryHighMass = 1
        
        if theInputs['muonTrigCutoff'] > 100 and theInputs['muonTrigUnc_HM'] > 0:
            if self.mH > theInputs['muonTrigCutoff']:
                self.sel_muontrig = theInputs['muonTrigUnc_HM']
                
        if theInputs['muonFullCutoff'] > 100 and theInputs['muonFullUnc_HM'] > 0:
            if self.mH > theInputs['muonFullCutoff']:
                self.sel_muonfull = theInputs['muonFullUnc_HM']
 
        if theInputs['elecTrigCutoff'] > 100 and theInputs['elecTrigUnc_HM'] > 0:
            if self.mH > theInputs['elecTrigCutoff']:
                self.sel_eletrig = theInputs['elecTrigUnc_HM']
                
        if theInputs['elecFullCutoff'] > 100 and theInputs['elecFullUnc_HM'] > 0:
            if self.mH > theInputs['elecFullCutoff']:
                self.sel_elefull = theInputs['elecFullUnc_HM']

        self.qqVV_scaleSys = 1. + 0.01*math.sqrt((self.mH - 20.)/13.)
        self.qqVV_pdfSys = 1. + 0.0035*math.sqrt(self.mH - 30.)
        self.ggVV_scaleSys = 1.04 + 0.10*math.sqrt((self.mH + 40.)/40.)
        self.ggVV_pdfSys = 1. + 0.0066*math.sqrt(self.mH - 10.)
            
          


    def setSystematics(self,theRateBkg_qqZZ,theRateBkg_ggZZ,theRateBkg_zjets ):

        self.rateBkg_qqZZ = theRateBkg_qqZZ
        self.rateBkg_ggZZ = theRateBkg_ggZZ
        self.rateBkg_zjets = theRateBkg_zjets

        if not self.model == "SM4" and not self.model == "SM" and not self.model == "FF" :
            
            print "In setSystematics in Systematics -------> Unknown model ",self.model
            print "Choices are SM, SM4, or FF"
            
        
        #ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        #ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")
        
        self.myCSW = HiggsCSandWidth()
        self.myCSWSM4 = HiggsCSandWidthSM4()

        if not self.isForXSxBR:
        
            self.CSpdfErrPlus_gg = self.myCSW.HiggsCSpdfErrPlus(1,self.mH,self.sqrts)
            self.CSpdfErrMinus_gg = self.myCSW.HiggsCSpdfErrMinus(1,self.mH,self.sqrts)
            self.CSpdfErrPlus_vbf = self.myCSW.HiggsCSpdfErrPlus(2,self.mH,self.sqrts)
            self.CSpdfErrMinus_vbf = self.myCSW.HiggsCSpdfErrMinus(2,self.mH,self.sqrts)
            self.CSpdfErrPlus_wh = self.myCSW.HiggsCSpdfErrPlus(3,self.mH,self.sqrts)
            self.CSpdfErrMinus_wh = self.myCSW.HiggsCSpdfErrMinus(3,self.mH,self.sqrts)
            self.CSpdfErrPlus_zh = self.myCSW.HiggsCSpdfErrPlus(4,self.mH,self.sqrts)
            self.CSpdfErrMinus_zh = self.myCSW.HiggsCSpdfErrMinus(4,self.mH,self.sqrts)
            self.CSpdfErrPlus_tth = self.myCSW.HiggsCSpdfErrPlus(5,self.mH,self.sqrts)
            self.CSpdfErrMinus_tth = self.myCSW.HiggsCSpdfErrMinus(5,self.mH,self.sqrts)
            
            self.CSscaleErrPlus_gg = self.myCSW.HiggsCSscaleErrPlus(1,self.mH,self.sqrts)
            self.CSscaleErrMinus_gg = self.myCSW.HiggsCSscaleErrMinus(1,self.mH,self.sqrts)
            self.CSscaleErrPlus_vbf = self.myCSW.HiggsCSscaleErrPlus(2,self.mH,self.sqrts)
            self.CSscaleErrMinus_vbf = self.myCSW.HiggsCSscaleErrMinus(2,self.mH,self.sqrts)
            self.CSscaleErrPlus_wh = self.myCSW.HiggsCSscaleErrPlus(3,self.mH,self.sqrts)
            self.CSscaleErrMinus_wh = self.myCSW.HiggsCSscaleErrMinus(3,self.mH,self.sqrts)
            self.CSscaleErrPlus_zh = self.myCSW.HiggsCSscaleErrPlus(4,self.mH,self.sqrts)
            self.CSscaleErrMinus_zh = self.myCSW.HiggsCSscaleErrMinus(4,self.mH,self.sqrts)
            self.CSscaleErrPlus_tth = self.myCSW.HiggsCSscaleErrPlus(5,self.mH,self.sqrts)
            self.CSscaleErrMinus_tth = self.myCSW.HiggsCSscaleErrMinus(5,self.mH,self.sqrts)
      

            if( self.mH >= 200): self.theoryHighMass = 1 + 1.5*(self.mH/1000)*(self.mH/1000)*(self.mH/1000)
      
            self.BRErr_Hff = self.myCSWSM4.HiggsBRErr_Hff(11,self.mH,self.sqrts)
            self.BRErr_HVV = self.myCSWSM4.HiggsBRErr_HVV(11,self.mH,self.sqrts)
            self.BRErr_Hgg = self.myCSWSM4.HiggsBRErr_Hgluglu(11,self.mH,self.sqrts)
      

        
        self.muSelError = 1 + math.sqrt( self.sel_muonfull*self.sel_muonfull + self.sel_muontrig*self.sel_muontrig )
        self.eSelError = 1 + math.sqrt( self.sel_elefull*self.sel_elefull + self.sel_eletrig*self.sel_eletrig )
        self.muSelError2e2mu = 1 + math.sqrt( self.sel_muonfull*self.sel_muonfull + self.sel_muontrig*self.sel_muontrig )
        self.eSelError2e2mu = 1 + math.sqrt( self.sel_elefull*self.sel_elefull + self.sel_eletrig*self.sel_eletrig )

    def Write_Systematics_Line(self,systLine,theFile,theInputs):
        print "~~~~~~~~~~~~~~~~~"
        channelList=['ggH','qqH','WH','ZH','ttH','ggZZ','ggZZ_signal','ggZZ_bkg','ggZZ_interf','qqZZ','zjets','ttbar','zbb']
        if theInputs["all"]:
            channelList=['ggH','ggZZ','ggZZ_signal','ggZZ_bkg','ggZZ_interf','qqZZ','zjets','ttbar','zbb']
        
        for chan in channelList:
            if theInputs[chan] or (chan.startswith("gg") and theInputs["all"]):
                print chan, systLine[chan]
                theFile.write(systLine[chan])
        theFile.write("\n")
        
    def Build_lumi(self,theFile,theInputs):
        if(self.sqrts == 7):
            theFile.write("lumi_7TeV lnN ")
        elif (self.sqrts == 8):
            theFile.write("lumi_8TeV lnN ")
        else:
            raise RuntimeError, "Unknown sqrts in systematics!"

        systLine={'ggH':"{0} ".format(self.lumiUncertainty)}
        systLine['qqH']  = "{0} ".format(self.lumiUncertainty)
        systLine['WH']   = "{0} ".format(self.lumiUncertainty)
        systLine['ZH']   = "{0} ".format(self.lumiUncertainty)
        systLine['ttH']  = "{0} ".format(self.lumiUncertainty)
        systLine['ggZZ'] = "{0} ".format(self.lumiUncertainty)
        systLine['ggZZ_signal'] = "{0} ".format(self.lumiUncertainty)
        systLine['ggZZ_bkg'] = "{0} ".format(self.lumiUncertainty)
        systLine['ggZZ_interf'] = "{0} ".format(self.lumiUncertainty)
        systLine['qqZZ'] = "{0} ".format(self.lumiUncertainty)
        systLine['zjets']= "- "
        systLine['ttbar']= "{0} ".format(self.lumiUncertainty)
        systLine['zbb']  = "{0} ".format(self.lumiUncertainty)

        self.Write_Systematics_Line(systLine,theFile,theInputs)
        
    def Write_pdf_qqbar(self,theFile,theInputs):

        theFile.write("pdf_qqbar lnN ")
            
        if not self.isForXSxBR:
            systLine={'ggH':"- "}
            systLine['qqH']  = "{0:.4f} ".format(1. + (self.CSpdfErrPlus_vbf-self.CSpdfErrMinus_vbf)/2.)
            systLine['WH']   = "{0:.4f} ".format(1. + (self.CSpdfErrPlus_wh-self.CSpdfErrMinus_wh)/2.)
            systLine['ZH']   = "{0:.4f} ".format(1. + (self.CSpdfErrPlus_zh-self.CSpdfErrMinus_zh)/2.)
            systLine['ttH']  = "- "
            systLine['ggZZ'] = "- "
            systLine['ggZZ_signal'] = "- "
            systLine['ggZZ_bkg'] = "- "
            systLine['ggZZ_interf'] = "- "
            systLine['qqZZ'] = "{0:.4f} ".format(self.qqVV_pdfSys)
            systLine['zjets']= "- "
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)
            
        elif self.isForXSxBR:
            systLine={'ggH':"- "}
            systLine['qqH']  = "- "
            systLine['WH']   = "- "
            systLine['ZH']   = "- "
            systLine['ttH']  = "- " 
            systLine['ggZZ'] = "- "
            systLine['ggZZ_signal'] = "- "
            systLine['ggZZ_bkg'] = "- "
            systLine['ggZZ_interf'] = "- "
            systLine['qqZZ'] = "{0:.4f} ".format(self.qqVV_pdfSys)
            systLine['zjets']= "- "
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)
                            
    def Write_pdf_gg(self,theFile,theInputs):

        theFile.write("pdf_gg lnN ")

        if not self.isForXSxBR and not self.model=="FF":
            systLine={'ggH':"{0:.4f} ".format(1. + (self.CSpdfErrPlus_gg-self.CSpdfErrMinus_gg)/2.)}
            systLine['qqH']  = "- "
            systLine['WH']   = "- " 
            systLine['ZH']   = "- "
            systLine['ttH']  = "{0:.4f} ".format(1. + (self.CSpdfErrPlus_tth-self.CSpdfErrMinus_tth)/2)
            systLine['qqZZ'] = "- " 
            systLine['ggZZ'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['ggZZ_signal'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['ggZZ_bkg'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['ggZZ_interf'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['zjets']= "- " 
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)
            
        elif self.isForXSxBR:
            systLine={'ggH':"- "}
            systLine['qqH']  = "- "
            systLine['WH']   = "- " 
            systLine['ZH']   = "- "
            systLine['ttH']  = "- " 
            systLine['qqZZ'] = "- " 
            systLine['ggZZ'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['ggZZ_signal'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['ggZZ_bkg'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['ggZZ_interf'] = "{0:.4f} ".format(self.ggVV_pdfSys)
            systLine['zjets']= "- " 
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)            
            
    def Write_pdf_hzz4l_accept(self,theFile,theInputs):
        theFile.write("pdf_hzz4l_accept lnN ")
        systLine={'ggH':"1.02 "}
        systLine['qqH']  = "1.02 "
        systLine['WH']   = "1.02 " 
        systLine['ZH']   = "1.02 "
        systLine['ttH']  = "1.02 "
        systLine['ggZZ_signal'] = "1.02 "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)            

    def Write_QCDscale_ggH(self,theFile,theInputs):
        
        theFile.write("QCDscale_ggH lnN ")

        systLine={'ggH':"{0:.4f} ".format(1. + (self.CSscaleErrPlus_gg-self.CSscaleErrMinus_gg)/2.)}
        systLine['qqH']  = "- "
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "- " 
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)            

    def Write_QCDscale_qqH(self,theFile,theInputs):

        theFile.write("QCDscale_qqH lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "{0:.4f} ".format(1. + (self.CSscaleErrPlus_vbf-self.CSscaleErrMinus_vbf)/2.)
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "- " 
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)            
        
    def Write_QCDscale_VH(self,theFile,theInputs):
        theFile.write("QCDscale_VH lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- "
        systLine['WH']   = "{0:.4f} ".format(1. + (self.CSscaleErrPlus_wh-self.CSscaleErrMinus_wh)/2.)
        systLine['ZH']   = "{0:.4f} ".format(1. + (self.CSscaleErrPlus_zh-self.CSscaleErrMinus_zh)/2.)
        systLine['ttH']  = "- " 
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)                    
        
    def Write_QCDscale_ttH(self,theFile,theInputs):

        theFile.write("QCDscale_ttH lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- "
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "{0:.4f} ".format(1. + (self.CSscaleErrPlus_tth-self.CSscaleErrMinus_tth)/2.)
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)                    
        
    def Write_QCDscale_ggVV(self,theFile,theInputs):

        theFile.write("QCDscale_ggVV lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- "
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "- "
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "{0:.4f} ".format(self.ggVV_scaleSys)
        systLine['ggZZ_signal'] = "{0:.4f} ".format(self.ggVV_scaleSys)
        systLine['ggZZ_bkg'] = "{0:.4f} ".format(self.ggVV_scaleSys)
        systLine['ggZZ_interf'] = "{0:.4f} ".format(self.ggVV_scaleSys)
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)                         

    def Write_QCDscale_VV(self,theFile,theInputs):

        theFile.write("QCDscale_VV lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- "
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "- "
        systLine['qqZZ'] = "{0:.4f} ".format(self.qqVV_scaleSys) 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)                         
        
    def Write_theoryUncXS_HighMH(self,theFile,theInputs):

        theFile.write("theoryUncXS_HighMH lnN ")

        systLine={'ggH':"{0:.3f} ".format(self.theoryHighMass)}
        systLine['qqH']  = "{0:.3f} ".format(self.theoryHighMass)
        systLine['WH']   = "{0:.3f} ".format(self.theoryHighMass)
        systLine['ZH']   = "{0:.3f} ".format(self.theoryHighMass)
        systLine['ttH']  = "{0:.3f} ".format(self.theoryHighMass)
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- "
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)
                
    def Write_BRhiggs_hzz4l(self,theFile,theInputs):

        theFile.write("BRhiggs_hzz4l lnN ")

        systLine={'ggH':"1.02 "}
        systLine['qqH']  = "1.02 " 
        systLine['WH']   = "1.02 " 
        systLine['ZH']   = "1.02 " 
        systLine['ttH']  = "1.02 " 
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- "
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)
        
    def Write_gamma_Hff(self,theFile,theInputs):
        theFile.write("gamma_Hff lnN ")
        
        systLine={'ggH':"{0:.4f} ".format(self.BRErr_Hff)}
        systLine['qqH']  = "{0:.4f} ".format(self.BRErr_Hff)
        systLine['WH']   = "{0:.4f} ".format(self.BRErr_Hff)
        systLine['ZH']   = "{0:.4f} ".format(self.BRErr_Hff)
        systLine['ttH']  = "{0:.4f} ".format(self.BRErr_Hff)
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- "
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)
        
    def Write_gamma_HVV(self,theFile,theInputs):
      
        theFile.write("gamma_HVV lnN ")

        systLine={'ggH':"{0:.4f} ".format(self.BRErr_HVV)}
        systLine['qqH']  = "{0:.4f} ".format(self.BRErr_HVV)
        systLine['WH']   = "{0:.4f} ".format(self.BRErr_HVV)
        systLine['ZH']   = "{0:.4f} ".format(self.BRErr_HVV)
        systLine['ttH']  = "{0:.4f} ".format(self.BRErr_HVV)
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- "
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)
        
    def Write_gamma_Hgluglu(self,theFile,theInputs):                    
        theFile.write("gamma_Hgluglu lnN ")

        systLine={'ggH':"{0:.4f} ".format(self.BRErr_Hgg)}
        systLine['qqH']  = "{0:.4f} ".format(self.BRErr_Hgg)
        systLine['WH']   = "{0:.4f} ".format(self.BRErr_Hgg)
        systLine['ZH']   = "{0:.4f} ".format(self.BRErr_Hgg)
        systLine['ttH']  = "{0:.4f} ".format(self.BRErr_Hgg)
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "- "
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_eff_e(self,theFile,theInputs):
        if self.channel == self.ID_4e:
            theFile.write("CMS_eff_e lnN ")

            systLine={'ggH':"{0:.3f} ".format(self.eSelError)}
            systLine['qqH']  = "{0:.3f} ".format(self.eSelError)
            systLine['WH']   = "{0:.3f} ".format(self.eSelError)
            systLine['ZH']   = "{0:.3f} ".format(self.eSelError)
            systLine['ttH']  = "{0:.3f} ".format(self.eSelError)
            systLine['qqZZ'] = "{0:.3f} ".format(self.eSelError)
            systLine['ggZZ'] = "{0:.3f} ".format(self.eSelError)
            systLine['ggZZ_signal'] = "{0:.3f} ".format(self.eSelError)
            systLine['ggZZ_bkg'] = "{0:.3f} ".format(self.eSelError)
            systLine['ggZZ_interf'] = "{0:.3f} ".format(self.eSelError)
            systLine['zjets']= "- "
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)
            
        elif self.channel == self.ID_2e2mu:
            theFile.write("CMS_eff_e lnN ")

            systLine={'ggH':"{0:.3f} ".format(self.eSelError2e2mu)}
            systLine['qqH']  = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['WH']   = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['ZH']   = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['ttH']  = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['qqZZ'] = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['ggZZ'] = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['ggZZ_signal'] = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['ggZZ_bkg'] = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['ggZZ_interf'] = "{0:.3f} ".format(self.eSelError2e2mu)
            systLine['zjets']= "- "
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_eff_m(self,theFile,theInputs):
        if self.channel == self.ID_4mu:
            theFile.write("CMS_eff_m lnN ")

            systLine={'ggH':"{0:.3f} ".format(self.muSelError)}
            systLine['qqH']  = "{0:.3f} ".format(self.muSelError)
            systLine['WH']   = "{0:.3f} ".format(self.muSelError)
            systLine['ZH']   = "{0:.3f} ".format(self.muSelError)
            systLine['ttH']  = "{0:.3f} ".format(self.muSelError)
            systLine['qqZZ'] = "{0:.3f} ".format(self.muSelError)
            systLine['ggZZ'] = "{0:.3f} ".format(self.muSelError)
            systLine['ggZZ_signal'] = "{0:.3f} ".format(self.muSelError)
            systLine['ggZZ_bkg'] = "{0:.3f} ".format(self.muSelError)
            systLine['ggZZ_interf'] = "{0:.3f} ".format(self.muSelError)
            systLine['zjets']= "- "
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)

        elif self.channel == self.ID_2e2mu:
            theFile.write("CMS_eff_m lnN ")

            systLine={'ggH':"{0:.3f} ".format(self.muSelError2e2mu)}
            systLine['qqH']  = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['WH']   = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['ZH']   = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['ttH']  = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['qqZZ'] = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['ggZZ'] = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['ggZZ_signal'] = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['ggZZ_bkg'] = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['ggZZ_interf'] = "{0:.3f} ".format(self.muSelError2e2mu)
            systLine['zjets']= "- "
            systLine['ttbar']= "- "
            systLine['zbb']  = "- "
            
            self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_CMS_hzz2e2mu_Zjets(self,theFile,theInputs):

        theFile.write("CMS_hzz2e2mu_Zjets lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- " 
        systLine['WH']   = "- " 
        systLine['ZH']   = "- " 
        systLine['ttH']  = "- " 
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "{0}/{1} ".format(self.zjetKappaLow,self.zjetKappaHigh)
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_CMS_hzz4mu_Zjets(self,theFile,theInputs):

        theFile.write("CMS_hzz4mu_Zjets lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- " 
        systLine['WH']   = "- " 
        systLine['ZH']   = "- " 
        systLine['ttH']  = "- " 
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "{0}/{1} ".format(self.zjetKappaLow,self.zjetKappaHigh)
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_CMS_hzz4e_Zjets(self,theFile,theInputs):
        
        theFile.write("CMS_hzz4e_Zjets lnN ")

        systLine={'ggH':"- "}
        systLine['qqH']  = "- " 
        systLine['WH']   = "- " 
        systLine['ZH']   = "- " 
        systLine['ttH']  = "- " 
        systLine['qqZZ'] = "- " 
        systLine['ggZZ'] = "- "
        systLine['ggZZ_signal'] = "- "
        systLine['ggZZ_bkg'] = "- "
        systLine['ggZZ_interf'] = "- "
        systLine['zjets']= "{0}/{1} ".format(self.zjetKappaLow,self.zjetKappaHigh)
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_CMS_zz4l_bkgMELA(self,theFile,theInputs):
        theFile.write("CMS_zz4l_bkgMELA param 0  1  [-3,3]\n")

    def Write_CMS_zz4l_sigMELA(self,theFile,theInputs):
        theFile.write("CMS_zz4l_sigMELA param 0  1  [-3,3]\n")

    def Write_CMS_zz4l_Jet_Split(self,theFile,theInputs,theVBFcat=False):
        theFile.write("QCDscale_ggH2in lnN ")
        if theVBFcat:
            systLine={'ggH':"{0:.3f} ".format(1+self.QCD_scale_ggH_2j_sys)}
        else:
            self.MH = ROOT.RooRealVar("MH","MH",self.mH)
            self.MH.setConstant(True)
            print theInputs['dijetRatio']
            rfv_dijetRatio = ROOT.RooFormulaVar("dijetRatio",theInputs['dijetRatio'],ROOT.RooArgList(self.MH))            
            systLine={'ggH':"{0:.3f} ".format(1-self.QCD_scale_ggH_2j_sys*rfv_dijetRatio.getVal())}
            #systLine={'ggH':"- "}
            
        systLine['qqH']  = "- "
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "- "
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)

        theFile.write("QCDscale_qqH2in lnN ")

        systLine={'ggH':"- "}
        if theVBFcat:
            systLine['qqH']  = "{0:.3f} ".format(1+self.QCD_scale_qqH_2j_sys)
        else:
            systLine['qqH']  = "- "
        systLine['WH']   = "- " 
        systLine['ZH']   = "- "
        systLine['ttH']  = "- "
        systLine['qqZZ'] = "- "
        systLine['ggZZ'] = "- "
        systLine['zjets']= "- " 
        systLine['ttbar']= "- "
        systLine['zbb']  = "- "
        
        self.Write_Systematics_Line(systLine,theFile,theInputs)

        #theFile.write("QCDscale_qqZZ_2j lnN ")

        #systLine={'ggH':"- "}
        #systLine['qqH']  = "- "
        #systLine['WH']   = "- " 
        #systLine['ZH']   = "- "
        #systLine['ttH']  = "- "
        #if theVBFcat:
        #    systLine['qqZZ'] = "{0.3f} ".format(1+self.QCD_scale_qqZZ_2j_sys)
        #else:
        #    systLine['qqZZ'] = "- "
        #systLine['ggZZ'] = "- "
        #systLine['zjets']= "- " 
        #systLine['ttbar']= "- "
        #systLine['zbb']  = "- "
        #
        #self.Write_Systematics_Line(systLine,theFile,theInputs)

    def Write_CMS_zz4l_Fisher_sys(self,theFile,theInputs):
        theFile.write("CMS_zz4l_ggH_Fisher_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_qqH_Fisher_sys param 0  1  [-3,3]\n")
        #theFile.write("CMS_zz4l_ttH_Fisher_sys param 0  1  [-3,3]\n")
        #theFile.write("CMS_zz4l_WH_Fisher_sys param 0  1  [-3,3]\n")
        #theFile.write("CMS_zz4l_ZH_Fisher_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_qqZZ_Fisher_sys param 0  1  [-3,3]\n")
        #theFile.write("CMS_zz4l_ggZZ_Fisher_sys param 0  1  [-3,3]\n")
        #theFile.write("CMS_zz4l_ZX_Fisher_sys param 0  1  [-3,3]\n")

    def Write_CMS_zz4l_Pt_sys(self,theFile,theInputs):
        theFile.write("CMS_zz4l_ggH_Pt_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_qqH_Pt_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_ttH_Pt_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_VH_Pt_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_qqZZ_Pt_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_ggZZ_Pt_sys param 0  1  [-3,3]\n")
        theFile.write("CMS_zz4l_ZX_Pt_sys param 0  1  [-3,3]\n")
    
    def WriteSystematics(self,theFile,theInputs, theVBFcat=False, theUse3D=False):

        #basic systematics

        if theInputs['useLumiUnc']:
            self.Build_lumi(theFile,theInputs)

        if theInputs['useCMS_eff']:
            self.Write_eff_m(theFile,theInputs)
            self.Write_eff_e(theFile,theInputs)

	if (self.channel == self.ID_4mu) and theInputs['useCMS_hzz4l_Zjets']:
            self.Write_CMS_hzz4mu_Zjets(theFile,theInputs)
	  
	if (self.channel == self.ID_4e) and theInputs['useCMS_hzz4l_Zjets']:
            self.Write_CMS_hzz4e_Zjets(theFile,theInputs)
        
	if (self.channel == self.ID_2e2mu) and theInputs['useCMS_hzz4l_Zjets']:
            self.Write_CMS_hzz2e2mu_Zjets(theFile,theInputs)

        #non-basic
        
        #if theInputs['usePdf_gg']:
        #    self.Write_pdf_gg(theFile,theInputs)

        if theInputs['usePdf_qqbar']:
            self.Write_pdf_qqbar(theFile,theInputs)
		
        if theInputs['usePdf_hzz4l_accept']:
            self.Write_pdf_hzz4l_accept(theFile,theInputs)
	    
        #non-non-basic

        #if not self.model == "FF" and theInputs['useQCDscale_ggH'] and not theVBFcat :
        #    self.Write_QCDscale_ggH(theFile,theInputs)

        #if theInputs['useQCDscale_qqH'] and not theVBFcat :
        #    self.Write_QCDscale_qqH(theFile,theInputs)

        #if theInputs['useQCDscale_VH']:
        #    self.Write_QCDscale_VH(theFile,theInputs)
	    
        #if not self.model == "FF" and theInputs['useQCDscale_ttH']:
        #    self.Write_QCDscale_ttH(theFile,theInputs)

        #if theInputs['useTheoryUncXS_HighMH']:
        #    self.Write_theoryUncXS_HighMH(theFile,theInputs)

        #if theInputs['useQCDscale_ggVV']:
        #    self.Write_QCDscale_ggVV(theFile,theInputs)

        if theInputs['useQCDscale_VV']:
            self.Write_QCDscale_VV(theFile,theInputs)
	
	## Higgs BR
        if(self.model == "SM" or self.model == "FF") and theInputs['useBRhiggs_hzz4l']:
            self.Write_BRhiggs_hzz4l(theFile,theInputs)
	  
	#elif(self.model == "SM4"):
        #    self.Write_gamma_Hff(theFile,theInputs)
        #    self.Write_gamma_HVV(theFile,theInputs)
        #    self.Write_gamma_Hgluglu(theFile,theInputs)
	  
	##  ----------- SELECTION EFFICIENCIES ----------

        #if theInputs['useCMS_zz4l_bkgMELA']:
        #    self.Write_CMS_zz4l_bkgMELA(theFile,theInputs)
            
        #if theInputs['useCMS_zz4l_sigMELA']:
        #    self.Write_CMS_zz4l_sigMELA(theFile,theInputs)

        #if theInputs['useCMS_zz4l_doVBFtest']:
        #    self.Write_CMS_zz4l_Jet_Split(theFile,theInputs, theVBFcat)

        #if (theVBFcat and theUse3D and theInputs['useCMS_zz4l_Fisher_sys']):
        #    self.Write_CMS_zz4l_Fisher_sys(theFile,theInputs)
            
        #if (not theVBFcat and theUse3D and theInputs['useCMS_zz4l_Pt_sys']):
        #    self.Write_CMS_zz4l_Pt_sys(theFile,theInputs)
            

    def WriteShapeSystematics(self,theFile,theInputs):
  
        #meanCB_e_errPerCent = theInputs['CMS_zz4l_mean_e_sig']
        #sigmaCB_e_errPerCent = theInputs['CMS_zz4l_sigma_e_sig']
        #N_CB_errPerCent = theInputs['CMS_zz4l_n_sig']
        #meanCB_m_errPerCent = theInputs['CMS_zz4l_mean_m_sig']
        #sigmaCB_m_errPerCent = theInputs['CMS_zz4l_sigma_m_sig']
        #Gamma_BW_errPerCent = theInputs['CMS_zz4l_gamma_sig']

        #theFile.write("CMS_zz4l_mu param 0.935929  -0.23/+0.26 \n") #stat only
        #theFile.write("CMS_zz4l_mu param 1.0 -0.255977/+0.30849 \n") #expected uncertainty
        theFile.write("CMS_zz4l_mu param 1.00  -0.24/+0.27 \n")
        #theFile.write("CMS_zz4l_mu param 0.93  -0.24/+0.26 \n")
        #theFile.write("CMS_zz4l_mu param 0.935929  -0.244308/+0.296613 \n") #stat + syst
        theFile.write("CMS_zz4l_kbkg param 1.0  0.10 \n")
        theFile.write("CMS_zz4l_scale_syst shape1 1 - - \n")
        #theFile.write("QCDScale shape 1 - - \n")
        #theFile.write("CMS_zz4l_pdfUnc shape 1 - - \n")
        #theFile.write("CMS_zz4l_syst param 0.0 1 [-3,3] \n") 
        #if( self.channel == self.ID_4mu):

            #if theInputs['useCMS_zz4l_mean']:
            #    theFile.write("CMS_zz4l_mean_m_sig param 0.0 1.0 \n")
            #    theFile.write("## CMS_zz4l_mean_m_sig = {0} \n".format(meanCB_m_errPerCent))
            #if theInputs['useCMS_zz4l_sigma']:
            #    theFile.write("CMS_zz4l_sigma_m_sig param 0.0 {0} \n".format(sigmaCB_m_errPerCent))
            #if theInputs['useCMS_zz4l_n']:
            #    theFile.write("CMS_zz4l_n_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,N_CB_errPerCent))
            #if theInputs['useCMS_zz4l_gamma']:
            #    theFile.write("interf_ggH param 0 1 [-1,1] \n".format(self.channel,self.sqrts,Gamma_BW_errPerCent))
                #theFile.write("CMS_zz4l_gamma_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,Gamma_BW_errPerCent))
            
        #if( self.channel == self.ID_4e):

            #if theInputs['useCMS_zz4l_mean']:
            #    theFile.write("CMS_zz4l_mean_e_sig param 0.0 1.0 \n")
            #    theFile.write("## CMS_zz4l_mean_e_sig = {0} \n".format(meanCB_e_errPerCent))
            #if theInputs['useCMS_zz4l_sigma']:
            #    theFile.write("CMS_zz4l_sigma_e_sig param 0.0 {0} \n".format(sigmaCB_e_errPerCent))
            #if theInputs['useCMS_zz4l_n']:
            #    theFile.write("CMS_zz4l_n_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,N_CB_errPerCent))
            #if theInputs['useCMS_zz4l_gamma']:
            #    theFile.write("interf_ggH param 0 1 [-1,1] \n".format(self.channel,self.sqrts,Gamma_BW_errPerCent))
                #theFile.write("CMS_zz4l_gamma_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,Gamma_BW_errPerCent))
            
        #if( self.channel == self.ID_2e2mu):

            #if theInputs['useCMS_zz4l_mean']:
            #    theFile.write("CMS_zz4l_mean_m_sig param 0.0 1.0 \n")
            #    theFile.write("## CMS_zz4l_mean_m_sig = {0} \n".format(meanCB_m_errPerCent))
            #    theFile.write("CMS_zz4l_mean_e_sig param 0.0 1.0 \n".format(meanCB_e_errPerCent))
            #    theFile.write("## CMS_zz4l_mean_e_sig = {0} \n".format(meanCB_e_errPerCent))
            #if theInputs['useCMS_zz4l_sigma']:
            #    theFile.write("CMS_zz4l_sigma_m_sig param 0.0 {0} \n".format(sigmaCB_m_errPerCent))
            #    theFile.write("CMS_zz4l_sigma_e_sig param 0.0 {0} \n".format(sigmaCB_e_errPerCent))
            #if theInputs['useCMS_zz4l_n']:
            #    theFile.write("CMS_zz4l_n_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,N_CB_errPerCent))
            #if theInputs['useCMS_zz4l_gamma']:
            #    theFile.write("interf_ggH param 0 1 [-1,1] \n".format(self.channel,self.sqrts,Gamma_BW_errPerCent))
                #theFile.write("CMS_zz4l_gamma_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,Gamma_BW_errPerCent))

