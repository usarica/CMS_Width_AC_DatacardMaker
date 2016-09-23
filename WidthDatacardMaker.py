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
      self.channel = self.theInputCard.decayChan
      self.theChannelName = self.theInputCard.decayChan

      self.theCategorizer = theCategorizer
      self.iCat = iCat
      self.theSystematizer = theSystematizer
      self.theSystVarsDict = self.theSystematizer.getVariableDict()
      self.theOutputDir = theOutputDir

      self.workspace = ROOT.RooWorkspace("w", "w")
      self.workspace.importClassCode(AsymPow.Class(),True)
      self.workspace.importClassCode(AsymQuad.Class(),True)
      self.workspace.importClassCode(RooqqZZPdf_v2.Class(), True)
      self.workspace.importClassCode(RooFormulaVar.Class(), True)
      self.workspace.importClassCode(RooRealFlooredSumPdf.Class(),True)
      self.workspace.importClassCode(VerticalInterpPdf.Class(),True)

      self.theEqnsMaker = theEqnsMaker
      # RooRealVars from the equations maker class
      self.theLumi = self.theEqnsMaker.theLumi
      self.muF = self.theEqnsMaker.muF
      self.muV = self.theEqnsMaker.muV
      self.kbkg_gg = self.theEqnsMaker.kbkg_gg
      self.kbkg_VBF = self.theEqnsMaker.kbkg_VBF
      self.fai1 = self.theEqnsMaker.fai1
      self.phiai1 = self.theEqnsMaker.phiai1
      self.phia1_gg = self.theEqnsMaker.phia1_gg
      self.phia1_VBF = self.theEqnsMaker.phia1_VBF
      self.varm4l = self.theEqnsMaker.varm4l
      self.varKD = self.theEqnsMaker.varKD
      self.varKD2 = self.theEqnsMaker.varKD2

      self.theBunches = [] # To keep track of which files are open
      self.pdfList = []
      self.rateList = []
      self.normList = []

      self.dataFileDir = "CMSdata"
      if (self.dataAppendDir != ''):
         self.dataFileDir = "{0}_{1}".format(self.dataFileDir,self.dataAppendDir)
      self.dataTreeName = "data_obs"
      self.dataFileName = "{0}/hzz{1}_{2}_{3:.0f}TeV.root".format(
         self.dataFileDir, self.theChannelName, self.theCategorizer.catNameList[self.iCat], self.sqrts
      )
      self.datacardName = "{0}/HCG/hzz{1}_{2}_{3:.0f}TeV_{4}.txt".format(
         self.outputDir, self.theChannelName, self.theCategorizer.catNameList[self.iCat], self.sqrts, self.appendName
      )
      self.workspaceFileName = "{0}/HCG/hzz{1}_{2}_{3:.0f}TeV_{4}.root".format(
         self.outputDir, self.theChannelName, self.theCategorizer.catNameList[self.iCat], self.sqrts, self.appendName
      )

      for proc in self.theInputCard.channels:
         procname = proc[0]
         proctype = proc[3]

         bunchNominal = None
         bunchVariations = []
         procPdf = None # Using VerticalInterpPdf
         procRate = None # Using AsymQuad
         procNorm = None # procRate*theLumi

         # Template file name core piece
         templateFileNameMain = "HtoZZ{}_{}_ModifiedSmoothTemplates_".format(self.theChannelName, self.theCategorizer.catNameList[self.iCat])

         # Nominal pdf and rate
         systName = "Nominal"
         templateFileName = "{0}/{1}{2}_{3}{4}".format(self.templateDir,templateFileNameMain,procname,systName,".root")
         if proctype>0:
            bunchNominal = BkgTemplateHelper(self.options,self,self.theCategorizer,procname,templateFileName,self.iCat,systName)
            bunchNominal.getTemplates()
         else:
            bunchNominal = BSITemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,procname,templateFileName,self.iCat,systName)
            bunchNominal.getTemplates(processName=procname)
         self.theBunches.append(bunchNominal)

         # Systematic "template" variations
         for syst in self.theInpurCard.systematics:
            systType = syst[1]
            if systType.lower() == "template":
               systVariableName = syst[0] # The name of variable that is supposed to exist for the variation
               for systchan in syst[2]: # Loop over channels that are relevant for this systematic variation
                  if systchan[0] == procname: # Make sure to pick the correct channel
                     tmplist = []
                     for isyst in range(1,3): # Up/Dn variations
                        systName = systchan[isyst]
                        templateFileName = "{0}/{1}{2}_{3}{4}".format(self.templateDir,templateFileNameMain,procname,systName,".root")
                        bunchVar = None
                        if proctype>0:
                           bunchVar = BkgTemplateHelper(self.options,self,self.theCategorizer,procname,templateFileName,self.iCat,systName)
                           bunchVar.getTemplates()
                        else:
                           bunchVar = BSITemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,procname,templateFileName,self.iCat,systName)
                           bunchVar.getTemplates(processName=procname)
                        if bunchVar is not None:
                           tmplist.append(bunchVar)
                           self.theBunches.append(bunchVar)
                     if len(tmplist)==2: # Expect exactly two templates
                        if self.theSystVarsDict[systVariableName] is not None:
                           tmplist.append(self.theSystVarsDict[systVariableName])
                           bunchVariations.append(tmplist)
                        else: raise RuntimeError("{} does not exist in the systematic template variables dictionary!".format(systVariableName))
                     else: raise RuntimeError("{} does not have exactly 2 template variations!".format(systVariableName))

         # Construct the ultimate pdf and rate for the process
         morphVarList = ROOT.RooArgList()
         morphPdfList = ROOT.RooArgList()
         morphRateList = ROOT.RooArgList()
         morphPdfList.add(bunchNominal.getThePdf())
         morphRateList.add(bunchNominal.getTheRate())
         for systvar in bunchVariations:
            for isyst in range(0,2):
               morphPdfList.add(systvar[isyst].getThePdf())
               morphRateList.add(systvar[isyst].getTheRate())
            morphVarList.add(systvar[2])
         procPdf = ROOT.VerticalInterpPdf(procname, procname, morphPdfList, morphVarList,1.0)
         ratename = bunchNominal.getTheRate().GetName() + "_AsymQuad"
         procRate = ROOT.AsymQuad(procname, procname, morphRateList, morphVarList, 1.0)
         normname = procname + "_norm"
         procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1,1e-15)", ROOT.RooArgList(procRate, self.theLumi))

         getattr(self.workspace, 'import')(procPdf,ROOT.RooFit.RecycleConflictNodes())
         self.pdfList.append(procPdf)
         getattr(self.workspace, 'import')(procNorm,ROOT.RooFit.RecycleConflictNodes())
         self.rateList.append(procRate)
         self.normList.append(procNorm)

      # Get the data
      self.theDataFile = ROOT.TFile.Open(self.dataFileName,"read")
      self.theDataTree = self.theDataFile.Get(self.dataTreeName)
      if not self.theDataTree:
         print "File, \"", self.dataFileName, "\" or tree \"", self.dataTreeName, "\" is not found."
      else:
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_full"
        data_obs = ROOT.RooDataSet(datasetName, datasetName, self.theDataTree, ROOT.RooArgSet(self.varm4l, self.varKD, self.varKD2))
        self.theDataRDS = data_obs.reduce("{0}>={1:.2f} && {0}<{1:.2f}".format(self.varm4l.GetName(), self.low_M, self.high_M))
        self.theDataRDS.SetName("data_obs")
        getattr(self.workspace, 'import')(self.theDataRDS, ROOT.RooFit.Rename("data_obs"))

      # Write datacards
      self.WriteDatacard()

      # Write the workspace
      #self.workspace.writeToFile(self.workspaceFileName)
      self.theWorkspaceFile = ROOT.TFile.Open(self.workspaceFileName,"recreate")
      self.theWorkspaceFile.cd()
      self.theWorkspaceFile.WriteTObject(self.workspace)
      self.theWorkspaceFile.Close()

      # Garbage collection
      for bunch in self.theBunches:
         if bunch.templateFile is not None:
            if bunch.templateFile.IsOpen():
               bunch.templateFile.Close()


   def WriteDatacard(self, theFile, self.theInputs, nameWS, theRates, obsEvents):
      self.theDatacardFile = open(self.datacardName, "wb")
      tmplist = self.workspaceFileName.split('/')
      nameWS = tmplist[len(tmplist)-1]

      nSigProcs = self.theInputCard.getNSigProcs()
      nBkgProcs = self.theInputCard.getNBkgProcs()

      self.theDatacardFile.write("imax {0:.0f}\n".format(nSigProcs))
      self.theDatacardFile.write("jmax {0:.0f}\n".format(nBkgProcs))
      self.theDatacardFile.write("kmax *\n")

      self.theDatacardFile.write("------------\n")
      self.theDatacardFile.write("shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
      self.theDatacardFile.write("------------\n")

      binname = "a{0:.0f}".format(self.channel)

      self.theDatacardFile.write("bin {} \n".format(binname))
      self.theDatacardFile.write("observation {0:.0f} \n".format(int(self.theDataRDS.numEntries())))

      self.theDatacardFile.write("------------\n")
      self.theDatacardFile.write("##########################################################################################################################\n")
      self.theDatacardFile.write("## Combine manual:                       https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit     ##\n")
      self.theDatacardFile.write("## Non-standard combine use cases:       https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/SWGuideNonStandardCombineUses ##\n")
      self.theDatacardFile.write("## Latest Higgs combination conventions: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions   ##\n")
      self.theDatacardFile.write("## CMS StatCom twiki:                    https://twiki.cern.ch/twiki/bin/view/CMS/StatisticsCommittee                   ##\n")
      self.theDatacardFile.write("##########################################################################################################################\n")
      self.theDatacardFile.write("## Mass window [{0:.1f},{1:.1f}] \n".format(self.low_M, self.high_M))
      #self.theDatacardFile.write("## signal,bkg,interf,tot rates [{0:.4f}, {1:.4f}, {2:.4f}, {3:.4f}] \n".format(
      #   theRates["ggZZ_signal"], theRates["ggZZbkg"], theRates["ggZZ_interf"], theRates["ggZZ_tot"])
      #)
      #self.theDatacardFile.write("## vbfsig,vbfbkg,vbfinterf,vbftot rates [{0:.4f}, {1:.4f}, {2:.4f}, {3:.4f}] \n".format(
      #   theRates["VBF_offshell_signal"], theRates["VBF_offshell_bkg"], theRates["VBF_offshell_interf"], theRates["VBF_offshell_tot"])
      #)

      self.theDatacardFile.write("bin ")
      for proc in self.theInputCard.channels:
         self.theDatacardFile.write("{} ".format(binname))
      self.theDatacardFile.write("\n")

      self.theDatacardFile.write("process ")
      for proc in self.theInputCard.channels:
         self.theDatacardFile.write("{} ".format(proc[0]))
      self.theDatacardFile.write("\n")

      self.theDatacardFile.write("process ")
      ctr_bkg = 0
      ctr_sig = -nSigProcs
      for proc in self.theInputCard.channels:
         if proc[3]>0:
            self.theDatacardFile.write("{0:.0f} ".format(ctr_bkg))
            ctr_bkg += 1
         else:
            self.theDatacardFile.write("{0:.0f} ".format(ctr_sig))
            ctr_sig += 1
      self.theDatacardFile.write("\n")

      self.theDatacardFile.write("rate ")
      for proc in self.theInputCard.channels:
         deflumi = proc[2]
         defrate = proc[1]
         if deflumi>0.:
            defrate = defrate/deflumi*self.theInputCard.lumi
         self.theDatacardFile.write("{0:.5f} ".format(defrate))
      self.theDatacardFile.write("\n")

      self.theDatacardFile.write("------------\n")
      self.theSystematizer.writeSystematics(self.theDatacardFile)
      self.theDatacardFile.close()



