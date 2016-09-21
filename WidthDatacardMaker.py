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
      for proc in self.theInputCard.channels:
         procname = proc[0]
         proctype = proc[3]

         bunchNominal = None
         bunchVariations = []
         procPdf = None # Using VerticalInterpPdf
         procRate = None # Using AsymQuad
         procNorm = None # procRate*theLumi

         # Template file name core piece
         templateFileNameMain = "HtoZZ4l_ModifiedSmoothTemplates_"

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

      # LEFT HERE: Now handle writing of systematics and importing of the data tree


      # Write Datacards
      fo = open(name_Shape, "wb")
      self.WriteDatacard(
      fo, self.theInputs, name_ShapeWS2, rates, data_obs_red.numEntries())

      systematics.WriteSystematics(fo, self.theInputs)
      systematics.WriteShapeSystematics(fo, self.theInputs)

      fo.close()

   def WriteDatacard(self, file, self.theInputs, nameWS, theRates, obsEvents):

      numberSig = self.numberOfSigChan(self.theInputs)
      numberBg = self.numberOfBgChan(self.theInputs)

      file.write("imax 1\n")
      file.write("jmax {0}\n".format(numberSig + numberBg - 1))
      file.write("kmax *\n")

      file.write("------------\n")
      file.write(
      "shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
      file.write("------------\n")

      file.write("bin a{0} \n".format(self.channel))
      file.write("observation {0} \n".format(obsEvents))

      file.write("------------\n")
      file.write(
      "## mass window [{0},{1}] \n".format(self.low_M, self.high_M))
      file.write("## signal,bkg,interf,tot rates [{0:.4f}, {1:.4f}, {2:.4f}, {3:.4f}] \n".format(
      theRates["ggZZ_signal"], theRates["ggZZbkg"], theRates["ggZZ_interf"], theRates["ggZZ_tot"]))
      file.write("## vbfsig,vbfbkg,vbfinterf,vbftot rates [{0:.4f}, {1:.4f}, {2:.4f}, {3:.4f}] \n".format(
      theRates["VBF_offshell_signal"], theRates["VBF_offshell_bkg"], theRates["VBF_offshell_interf"], theRates["VBF_offshell_tot"]))
      file.write("bin ")

      # channelList=['ggZZ_signal','ggZZ_interf','ggZZbkg','qqZZ','zjets']
      channelList = ['ggZZ', 'VBF_offshell', 'qqZZ', 'zjets']

      # channelName=['ggsignalzz','gginterfzz','ggZZbkg','bkg_qqzz','bkg_zjets']
      channelName = ['ggzz', 'vbf_offshell', 'bkg_qqzz', 'bkg_zjets']

      for chan in channelList:
      if self.theInputs[chan]:
      file.write("a{0} ".format(self.channel))
      file.write("\n")

      file.write("process ")

      i = 0

      for chan in channelList:
      if self.theInputs[chan]:
      file.write("{0} ".format(channelName[i]))
      i += 1

      file.write("\n")

      processLine = "process "

      for x in range(-numberSig + 1, 1):
      processLine += "{0} ".format(x)

      for y in range(1, numberBg + 1):
      processLine += "{0} ".format(y)

      file.write(processLine)
      file.write("\n")

      file.write("rate ")
      for chan in channelList:
      if self.theInputs[chan]:
      file.write("{0:.4f} ".format(theRates[chan]))
      file.write("\n")
      file.write("------------\n")

   def numberOfSigChan(self, inputs):

      counter = 0

      if inputs['ggZZ']:
      counter += 1
      if inputs['ggZZ_signal']:
      counter += 1
      if inputs['ggZZ_interf']:
      counter += 1
      if inputs['VBF_offshell']:
      counter += 1

      return counter

   def numberOfBgChan(self, inputs):

      counter = 0

      if inputs['qqZZ']:
      counter += 1
      if inputs['ggZZbkg']:
      counter += 1
      if inputs['zjets']:
      counter += 1
      if inputs['ttbar']:
      counter += 1
      if inputs['zbb']:
      counter += 1

      return counter
