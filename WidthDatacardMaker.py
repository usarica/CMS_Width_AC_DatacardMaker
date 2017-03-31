#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array
from InputCardReader import InputCardReader
from CategoryHelper import CategoryHelper
from EquationsMaker import EquationsMaker
from SystematicsHelper import SystematicsHelper
from SystematicsHelper import FloatToString
from BSITemplateHelper import BSITemplateHelper
from BkgTemplateHelper import BkgTemplateHelper


class WidthDatacardMaker:
   def __init__(self, options, theInputCard, theEqnsMaker, theCategorizer, theSystematizer, iCat, theOutputDir):
      self.iCat = iCat
      self.theOutputDir = theOutputDir

      self.options = options
      self.theInputCard = theInputCard
      self.theEqnsMaker = theEqnsMaker
      self.theCategorizer = theCategorizer
      self.theSystematizer = theSystematizer

      self.templateDir = self.options.templateDir
      self.dataAppendDir = self.options.dataDirAppend
      self.iCatScheme = self.options.iCatScheme
      self.mH = self.options.mPOLE
      self.mLow = self.options.mLow
      self.mHigh = self.options.mHigh
      self.coordList = options.coordList

      self.sqrts = self.theInputCard.sqrts
      self.channel = self.theInputCard.decayChan
      self.theChannelName = self.theInputCard.decayChanName
      self.catName = self.theCategorizer.catNameList[self.iCat]

      self.theSystVarsDict = self.theSystematizer.getVariableDict()

      self.workspace = ROOT.RooWorkspace("w", "w")
      self.workspace.importClassCode(ROOT.AsymPow.Class(),True)
      self.workspace.importClassCode(ROOT.AsymQuad.Class(),True)
      self.workspace.importClassCode(ROOT.RooqqZZPdf_v2.Class(), True)
      self.workspace.importClassCode(ROOT.RooFormulaVar.Class(), True)
      self.workspace.importClassCode(ROOT.FastHistoFunc_f.Class(),True)
      self.workspace.importClassCode(ROOT.FastHisto2DFunc_f.Class(),True)
      self.workspace.importClassCode(ROOT.FastHisto3DFunc_f.Class(),True)
      self.workspace.importClassCode(ROOT.RooRealFlooredSumPdf.Class(),True)
      self.workspace.importClassCode(ROOT.VerticalInterpPdf.Class(),True)

      # Other input-independent RooRealVars are taken from the equations maker class
      eqnrrvars = dict(self.theEqnsMaker.rrvars)
      self.theLumi = eqnrrvars["lumi"]
      self.mass = eqnrrvars["mass"]

      self.onevar = ROOT.RooConstVar("VarOne","VarOne",1.0)

      self.KD1=None
      self.KD2=None
      self.KD3=None
      self.dataVars = ROOT.RooArgSet()
      self.dataVars.add(self.mass)
      for coord in self.coordList:
         if self.KD3 is not None:
            raise RuntimeError("There are >3 KDs in the list of coordinates, which is nt cupported.")
         for key, value in eqnrrvars.iteritems():
            if key==coord:
               if key!="mass":
                  self.dataVars.add(value)
               if self.KD1 is None:
                  self.KD1=value
               elif self.KD2 is None:
                  self.KD2=value
               elif self.KD3 is None:
                  self.KD3=value
               else:
                  sys.exit("Too many KDs!")

      self.theBunches = [] # To keep track of which files are open
      self.pdfList = []
      self.rateList = []
      self.extraRateList = []
      self.normList = []
      self.extraVars = [] # To keep track of which variables are created on the fly

      self.dataFileDir = "CMSdata"
      if (self.dataAppendDir != ''):
         self.dataFileDir = "{0}_{1}".format(self.dataFileDir,self.dataAppendDir)
      self.dataTreeName = "data_obs"
      self.dataFileName = "{0}/hzz{1}_{2}_{3:.0f}TeV.root".format(
         self.dataFileDir, self.theChannelName, self.catName, self.sqrts
      )
      self.datacardName = "{0}/HCG/{3:.0f}TeV/hzz{1}_{2}.txt".format(
         self.theOutputDir, self.theChannelName, self.catName, self.sqrts
      )
      self.workspaceFileName = "{0}/HCG/{3:.0f}TeV/hzz{1}_{2}.input.root".format(
         self.theOutputDir, self.theChannelName, self.catName, self.sqrts
      )

      for proc in self.theInputCard.channels:
         procname = proc[0]
         proctype = proc[3]

         bunchNominal = None
         bunchVariations = []
         procPdf = None # Using VerticalInterpPdf
         procRate = None # Using AsymQuad
         procNorm = None # procRate*theLumi

         quadNVars = []

         # Template file name core piece
         templateFileNameMain = "HtoZZ{}_{}_FinalTemplates_".format(self.theChannelName, self.catName)

         # Nominal pdf and rate
         systName = "Nominal"
         templateFileName = "{0}/{1:.0f}TeV/{2}{3}_{4}{5}".format(self.templateDir, self.sqrts, templateFileNameMain, procname, systName, ".root")
         if proctype==0:
            bunchNominal = BkgTemplateHelper(self.options,self,self.theCategorizer,procname,templateFileName,self.iCat,systName)
            bunchNominal.getTemplates()
         else:
            bunchNominal = BSITemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,procname,proctype,templateFileName,self.iCat,systName)
            bunchNominal.getTemplates(processName=procname)
         self.theBunches.append(bunchNominal)

         # Systematic "template" variations
         for syst in self.theInputCard.systematics:
            systType = syst[1]
            if systType.lower() == "template":
               systVariableName = syst[0] # The name of variable that is supposed to exist for the variation
               for systchan in syst[2]: # Loop over channels that are relevant for this systematic variation
                  if systchan[0] == procname: # Make sure to pick the correct channel
                     tmplist = []
                     for isyst in range(1,3): # Up/Dn variations
                        systName = systVariableName
                        if isyst==1:
                           systName = "{}Up".format(systName)
                        else:
                           systName = "{}Down".format(systName)
                        templateFileName = "{0}/{1:.0f}TeV/{2}{3}_{4}{5}".format(self.templateDir,self.sqrts,templateFileNameMain,procname,systName,".root")
                        bunchVar = None
                        if proctype==0:
                           bunchVar = BkgTemplateHelper(self.options,self,self.theCategorizer,procname,templateFileName,self.iCat,systName)
                           bunchVar.getTemplates()
                        else:
                           bunchVar = BSITemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,procname,proctype,templateFileName,self.iCat,systName)
                           bunchVar.getTemplates(processName=procname)
                        tmplist.append(bunchVar)
                        self.theBunches.append(bunchVar)
                     if len(tmplist)==2: # Expect exactly two templates
                        if self.theSystVarsDict[systVariableName] is not None:
                           tmplist.append(self.theSystVarsDict[systVariableName]) # Append the systematic control RooRealVar as the last element
                           bunchVariations.append(tmplist)
                        else: raise RuntimeError("The RooRealVar {} does not exist in the systematic template variables dictionary!".format(systVariableName))
                     else: raise RuntimeError("{} does not have exactly 2 template variations!".format(systVariableName))
                     break
            elif systType.lower() == "quadn":
               systVariableName = syst[0] # The name of variable that is supposed to exist for the variation
               for systchan in syst[2]: # Loop over channels that are relevant for this systematic variation
                  if systchan[0] == procname: # Make sure to pick the correct channel
                     if(len(systchan)!=3):
                        raise RuntimeError("{} variation for {} needs to have 2 variations!".format(systVariableName, procname))
                     tmplist = []
                     for isyst in range(1,3): # Up/Dn variations
                        systName = ""
                        if (isyst==1):
                           systName = "{0}UpConst_{1}_{2}_{3}_{4:.0f}TeV".format(systVariableName,procname,self.catName,self.theChannelName,self.sqrts)
                        else:
                           systName = "{0}DownConst_{1}_{2}_{3}_{4:.0f}TeV".format(systVariableName,procname,self.catName,self.theChannelName,self.sqrts)
                        tmpvar = ROOT.RooConstVar(systName,systName,float(systchan[isyst]))
                        tmplist.append(tmpvar)
                        self.extraVars.append(tmpvar)
                     if len(tmplist)==2: # Expect exactly two templates
                        if self.theSystVarsDict[systVariableName] is not None:
                           tmplist.append(self.theSystVarsDict[systVariableName]) # Append the systematic control RooRealVar as the last element
                           quadNVars.append(tmplist)
                        else: raise RuntimeError("The RooRealVar {} does not exist in the systematic quadN variables dictionary!".format(systVariableName))
                     else: raise RuntimeError("{} does not have exactly 2 quadN variations!".format(systVariableName))
                     break

         procRateExtra = None
         if len(quadNVars)>0:
            extraMorphVarList = ROOT.RooArgList()
            extraMorphRateList = ROOT.RooArgList()
            self.extraVars.append(self.onevar)
            extraMorphRateList.add(self.onevar)
            for systvar in quadNVars:
               for isyst in range(0,2):
                  extraMorphRateList.add(systvar[isyst])
               extraMorphVarList.add(systvar[2])
            ratename = bunchNominal.getTheRate().GetName() + "_ExtraAsymQuad"
            procRateExtra = ROOT.AsymQuad(ratename, ratename, extraMorphRateList, extraMorphVarList, 1.0)
            self.extraRateList.append(procRateExtra)

         # Construct the ultimate pdf and rate for the process
         procPdf = None
         procRate = None
         if len(bunchVariations)>0:
            print "Bunch variations are found. Constructing VerticalInterpPdf pdf and AsymQuad rate..."
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
            procRate = ROOT.AsymQuad(ratename, ratename, morphRateList, morphVarList, 1.0)
         else:
            print "Bunch variations do not exist. Constructing pdf and rate from nominal bunch..."
            procPdf = bunchNominal.getThePdf()
            procPdf.SetName(procname)
            procPdf.SetTitle(procname)
            procRate = bunchNominal.getTheRate()

         self.pdfList.append(procPdf)
         self.rateList.append(procRate)

         normname = procname + "_norm"
         procNorm = None
         if procRateExtra is None:
            procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1,1e-15)", ROOT.RooArgList(procRate, self.theLumi))
         else:
            procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1*@2,1e-15)", ROOT.RooArgList(procRate, procRateExtra, self.theLumi))
         self.normList.append(procNorm)

         print "Last check on pdf value:",procPdf.getVal(),"at"
         if self.KD1 is not None:
            print "\t- KD1 =",self.KD1.getVal()
         if self.KD2 is not None:
            print "\t- KD2 =",self.KD2.getVal()
         if self.KD3 is not None:
            print "\t- KD3 =",self.KD3.getVal()
         print "Last check on pdf norm:",procNorm.getVal()
         getattr(self.workspace, 'import')(procPdf,ROOT.RooFit.RecycleConflictNodes())
         getattr(self.workspace, 'import')(procNorm,ROOT.RooFit.RecycleConflictNodes())

      # Get the data
      data_obs = ROOT.RooDataSet("dummy", "dummy", self.dataVars)
      self.theDataFile = ROOT.TFile.Open(self.dataFileName,"read")
      if not self.theDataFile:
         print "Data tree is not found!"
         self.theDataRDS = data_obs
      else:
         self.theDataTree = self.theDataFile.Get(self.dataTreeName)
         if not self.theDataTree:
            print "File, \"", self.dataFileName, "\" or tree \"", self.dataTreeName, "\" is not found."
            self.theDataRDS = data_obs
         else:
            datasetName = "data_obs_full"
            data_obs = ROOT.RooDataSet(datasetName, datasetName, self.theDataTree, self.dataVars)
            self.theDataRDS = data_obs.reduce("{0}>={1} && {0}<{2}".format(self.mass.GetName(), FloatToString(self.mLow), FloatToString(self.mHigh)))
      self.theDataRDS.SetName("data_obs")
      getattr(self.workspace, 'import')(self.theDataRDS, ROOT.RooFit.Rename("data_obs"))

      self.workspace.Print("v")
      print "Now writing workspace"
      # Write the workspace
      self.theWorkspaceFile = ROOT.TFile.Open(self.workspaceFileName,"recreate")
      self.theWorkspaceFile.WriteTObject(self.workspace)
      self.theWorkspaceFile.Close()

      print "Testing the workspace"
      ftmp = ROOT.TFile.Open(self.workspaceFileName,"read")
      wtmp = ftmp.Get("w")
      wtmp.Print("v")
      ftmp.Close()

      print "Now writing datacards"

      # Write datacards
      self.WriteDatacard()

      # Garbage collection
      for bunch in self.theBunches:
         bunch.close()


   def WriteDatacard(self):
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
      self.theDatacardFile.write("## Mass window [{0:.1f},{1:.1f}] \n".format(self.mLow, self.mHigh))
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
         if proc[3]==0:
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



