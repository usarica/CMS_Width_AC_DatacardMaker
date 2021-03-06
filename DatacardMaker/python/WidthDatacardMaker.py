#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array
from CMS_Width_AC_DatacardMaker.DatacardMaker.InputCardReader import InputCardReader
from CMS_Width_AC_DatacardMaker.DatacardMaker.CategoryHelper import CategoryHelper
from CMS_Width_AC_DatacardMaker.DatacardMaker.EquationsMaker import EquationsMaker
from CMS_Width_AC_DatacardMaker.DatacardMaker.SystematicsHelper import SystematicsHelper
from CMS_Width_AC_DatacardMaker.DatacardMaker.WidthHelperFunctions import FloatToString
from CMS_Width_AC_DatacardMaker.DatacardMaker.BSITemplateHelper import BSITemplateHelper
from CMS_Width_AC_DatacardMaker.DatacardMaker.SimpleTemplateHelper import SimpleTemplateHelper
from CMS_Width_AC_DatacardMaker.DatacardMaker.ExternalShapeHelper import ExternalShapeHelper


def FindMinMax(rate,var):
   varmin = var.getBinning().lowBound()
   varmax = var.getBinning().highBound()
   varval = var.getVal()
   npoints = 101
   ratemin = 999999999.
   ratemax = -999999999.
   for ix in range(0,npoints):
      var.setVal(varmin + float(ix)*(varmax-varmin)/float(npoints-1))
      rateval = rate.getVal()
      if ratemin>rateval:
         ratemin=rateval
      if ratemax<rateval:
         ratemax=rateval
   var.setVal(varval)
   return (ratemin,ratemax)

def PlotRate(rate,args,path):
   for xvar in args:
      canvasname = "c_{}_{}".format(rate.GetName(),xvar.GetName())
      cproj = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      (rmin,rmax) = FindMinMax(rate,xvar)
      plot = xvar.frame()
      rate.plotOn(plot, ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(1))
      plot.GetXaxis().CenterTitle()
      plot.GetYaxis().SetTitleOffset(1.2)
      plot.GetYaxis().CenterTitle()
      plot.GetXaxis().SetTitle(xvar.GetName())
      plot.GetYaxis().SetTitle(rate.GetName())
      plot.GetXaxis().SetNdivisions(510)
      plot.SetTitle("{} along {}".format(rate.GetName(),xvar.GetName()))
      plot.GetYaxis().SetRangeUser(
         rmin,
         rmax
      )
      if rmax>=20.*rmin and rmin>0.:
         cproj.SetLogy()
      plot.Draw()
      cproj.SaveAs("{}{}{}".format(path,canvasname,".png"))
      cproj.Close()
      del plot


def CheckPdf(pdf,xvar,yvar,zvar,checkvar,path):
   if pdf.dependsOn(checkvar):
      print("Checking pdf {} against {}".format(pdf.GetName(),checkvar.GetName()))

      xcheckdefval=0
      ycheckdefval=0
      zcheckdefval=0
      nbinsx=0
      nbinsy=0
      nbinsz=0
      xcheckvals=[]
      ycheckvals=[]
      zcheckvals=[]
      if xvar is not None:
         xcheckdefval=xvar.getVal()
         varBinning=xvar.getBinning()
         varBinArray=varBinning.array()
         nbinsx=varBinning.numBins()
         for ix in range(0,nbinsx):
            lowedge=varBinArray[ix]
            highedge=varBinArray[ix+1]
            xcheckvals.append((highedge+lowedge)/2.)
      else: xcheckvals.append(xcheckdefval)
      if yvar is not None:
         ycheckdefval=yvar.getVal()
         varBinning=yvar.getBinning()
         varBinArray=varBinning.array()
         nbinsy=varBinning.numBins()
         for iy in range(0,nbinsy):
            lowedge=varBinArray[iy]
            highedge=varBinArray[iy+1]
            ycheckvals.append((highedge+lowedge)/2.)
      else: ycheckvals.append(ycheckdefval)
      if zvar is not None:
         zcheckdefval=zvar.getVal()
         varBinning=zvar.getBinning()
         varBinArray=varBinning.array()
         nbinsz=varBinning.numBins()
         for iz in range(0,nbinsz):
            lowedge=varBinArray[iz]
            highedge=varBinArray[iz+1]
            zcheckvals.append((highedge+lowedge)/2.)
      else: zcheckvals.append(zcheckdefval)

      for xval in xcheckvals:
         if xvar is not None: xvar.setVal(xval)
         for yval in ycheckvals:
            if yvar is not None: yvar.setVal(yval)
            for zval in zcheckvals:
               if zvar is not None: zvar.setVal(zval)
               plotname="c_CheckPdf_{}_{}".format(pdf.GetName(),checkvar.GetName())
               canvasname = plotname
               if xvar is not None: plotname="{}_X{}".format(plotname,FloatToString(xval))
               if yvar is not None: plotname="{}_Y{}".format(plotname,FloatToString(yval))
               if zvar is not None: plotname="{}_Z{}".format(plotname,FloatToString(zval))

               print("Adding {} to {}".format(plotname,canvasname))

               cproj = ROOT.TCanvas( canvasname, canvasname, 200, 200 )
               plot = checkvar.frame()
               pdf.plotOn(plot, ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(1))
               plot.GetXaxis().CenterTitle()
               plot.GetYaxis().SetTitleOffset(1.2)
               plot.GetYaxis().CenterTitle()
               plot.GetXaxis().SetTitle(checkvar.GetName())
               plot.GetYaxis().SetTitle(pdf.GetName())
               plot.GetXaxis().SetNdivisions(510)
               plot.SetTitle(plotname)
               plot.Draw()
               savename="{}{}{}".format(path,canvasname,".gif+")
               cproj.SaveAs(savename)
               cproj.Close()

               if zvar is not None: zvar.setVal(zcheckdefval)
            if yvar is not None: yvar.setVal(ycheckdefval)
         if xvar is not None: xvar.setVal(xcheckdefval)



def PlotPdf1D(pdf,norm,xvar,yvar,zvar,systarg,path,appendname=""):
   varlist=[]
   if xvar is not None: varlist.append(xvar)
   if yvar is not None: varlist.append(yvar)
   if zvar is not None: varlist.append(zvar)
   ndims=len(varlist)
   for ivar in range(0,ndims):
      varplot=varlist[ivar]
      if systarg is not None:
         canvasname = "c_{}_pdf_{}_KD{}".format(pdf.GetName(),systarg.GetName(),ivar+1)
      else:
         canvasname = "c_{}_pdf_Nominal_KD{}".format(pdf.GetName(),ivar+1)
      if appendname!="":
         canvasname = canvasname+"_"+appendname
      cproj = ROOT.TCanvas( canvasname, canvasname, 750, 700 )
      plot = varplot.frame()
      projvars = ROOT.RooArgSet()
      for jvar in range(0,ndims):
         if jvar==ivar: continue
         projvars.add(varlist[jvar])
      pdf.plotOn(plot, ROOT.RooFit.Project(projvars),ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.LineColor(ROOT.kBlack))
      if systarg is not None:
         defval = systarg.getVal()
         setvals = [ -1., 1. ]
         vallabels = [ "Down", "Up" ]
         for argval,arglabel in zip(setvals,vallabels):
            systarg.setVal(argval)
            applabel = arglabel
            lcolor=ROOT.kBlack
            if arglabel=="Down":
               lcolor=ROOT.kBlue
            else:
               lcolor=ROOT.kRed
            pdf.plotOn(plot, ROOT.RooFit.Project(projvars),ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.NumEvent),ROOT.RooFit.LineColor(lcolor),ROOT.RooFit.LineStyle(7))
         systarg.setVal(defval)
      plot.Draw()
      cproj.SaveAs("{}{}{}".format(path,canvasname,".png"))
      cproj.Close()
      del plot


class WidthDatacardMaker:
   def __init__(self, options, theInputCard, theEqnsMaker, theCategorizer, theSystematizer, iCat, theOutputDir):
      self.iCat = iCat
      self.theOutputDir = theOutputDir

      self.options = options
      self.theInputCard = theInputCard
      self.theEqnsMaker = theEqnsMaker
      self.theCategorizer = theCategorizer
      #self.iCatScheme = self.theCategorizer.iCatScheme
      self.theSystematizer = theSystematizer

      self.templateDir = self.options.templateDir
      self.customDataDir = self.options.customDataDir
      self.mH = self.options.mPOLE
      self.mLow = self.options.mLow
      self.mHigh = self.options.mHigh
      self.coordList = options.coordList

      self.sqrts = self.theInputCard.sqrts
      self.theSqrtsPeriod = self.theInputCard.theSqrtsPeriod
      self.channel = self.theInputCard.decayChan
      self.theChannelName = self.theInputCard.decayChanName
      self.catName = self.theCategorizer.catNameList[self.iCat]

      self.theSystVarsDict = self.theSystematizer.getVariableDict()

      self.workspace = ROOT.RooWorkspace("w", "w")
      self.workspace.importClassCode(ROOT.AsymPow.Class(),True)
      self.workspace.importClassCode(ROOT.AsymQuad.Class(),True)
      self.workspace.importClassCode(ROOT.RooFormulaVar.Class(), True)
      #self.workspace.importClassCode(ROOT.RooqqZZPdf_v2.Class(), True)
      self.workspace.importClassCode(ROOT.RooDoubleCB.Class(), True)
      self.workspace.importClassCode(ROOT.FastHistoFunc_f.Class(),True)
      self.workspace.importClassCode(ROOT.FastHisto2DFunc_f.Class(),True)
      self.workspace.importClassCode(ROOT.FastHisto3DFunc_f.Class(),True)
      self.workspace.importClassCode(ROOT.RooRealFlooredSumPdf.Class(),True)
      self.workspace.importClassCode(ROOT.VerticalInterpPdf.Class(),True)

      # Other input-independent RooRealVars are taken from the equations maker class
      eqnrrvars = self.theEqnsMaker.rrvars
      self.theLumi = eqnrrvars["lumi"]
      self.mass = eqnrrvars["mass"]

      # External mass shapes
      self.extMassShapesDir=self.options.extMassShapes
      self.hasExtMassShapes=(self.extMassShapesDir is not None)
      self.extShapeHandle=None
      if self.hasExtMassShapes:
         shapesFileNameMain = "Hto{}_{}_FinalMassShape_".format(self.theChannelName, self.catName)
         shapesFileName = "{0}/{1}{2}{3}".format(self.extMassShapesDir, shapesFileNameMain, "AllProcesses", ".root")
         self.extShapeHandle=ExternalShapeHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,shapesFileName,self.iCat)

      self.onevar = eqnrrvars["one"]

      # Determine all observables and rename them for mass range, final state, category and data period
      addVarMassName = "i{0}_e{1}_{2}_{3}".format(
         FloatToString(self.mLow), FloatToString(self.mHigh),
         self.catName, self.theSqrtsPeriod
      )
      addVarMassName.replace(".","p")
      addVarName = "{}_{}".format("mass", addVarMassName)

      self.KD1=None
      self.KD2=None
      self.KD3=None
      self.dataVars = ROOT.RooArgSet()
      KDsHaveMass=False
      for coord in self.coordList:
         if self.KD3 is not None:
            raise RuntimeError("There are >3 KDs in the list of coordinates, which is not cupported.")
         for key, value in eqnrrvars.iteritems():
            if key==coord:
               repName=addVarName
               if "mass" in key:
                  repName = addVarMassName
                  KDsHaveMass=True
               newVarName = "{}_{}".format(value.GetName(),repName)
               if not("CMS_{}".format(self.theChannelName) in newVarName):
                  newVarName = "CMS_{}_{}".format(self.theChannelName,newVarName)
               print("Renaming",value.GetName(),"to",newVarName)
               value.SetName(newVarName)

               self.dataVars.add(value)
               if self.KD1 is None:
                  self.KD1=value
               elif self.KD2 is None:
                  self.KD2=value
               elif self.KD3 is None:
                  self.KD3=value
               else:
                  sys.exit("Too many KDs!")

      if not KDsHaveMass:
         newVarName = "{}_{}".format(self.mass.GetName(),addVarMassName)
         if not("CMS_{}".format(self.theChannelName) in newVarName):
            newVarName = "CMS_{}_{}".format(self.theChannelName,newVarName)
         print("Renaming",self.mass.GetName(),"to",newVarName)
         self.mass.SetName(newVarName)

      self.theBunches = [] # To keep track of which files are open
      self.pdfList = []
      self.rateList = []
      self.extraRateList = []
      self.extraPdfList = []
      self.normList = []
      self.extraVars = [] # To keep track of which variables are created on the fly

      self.dataFileDir = "CMSdata"
      if (self.customDataDir != ''):
         self.dataFileDir = self.customDataDir
      self.dataTreeName = "data_obs"
      self.dataFileName = "{0}/hto{1}_{2}_{3}.root".format(
         self.dataFileDir, self.theChannelName, self.catName, self.theSqrtsPeriod
      )
      self.datacardName = "{0}/HCG/{3}/hto{1}_{2}.txt".format(
         self.theOutputDir, self.theChannelName, self.catName, self.theSqrtsPeriod
      )
      self.workspaceFileName = "{0}/HCG/{3}/hto{1}_{2}.input.root".format(
         self.theOutputDir, self.theChannelName, self.catName, self.theSqrtsPeriod
      )
      self.plotsPathName = "{0}/figs/{3}/hto{1}_{2}/".format(
         self.theOutputDir, self.theChannelName, self.catName, self.theSqrtsPeriod
      )

      globalCondDim=1
      condDimsAdjusted=False
      for proc in self.theInputCard.channels:
         procname = proc[0]
         procopts = proc[4]
         isConditional=False
         procnamefile = procname
         condDim=0
         for procopt in procopts:
            procoptl=procopt.lower()
            if "conditional" in procoptl:
               isConditional=True
               condNormVars = ROOT.RooArgSet()
               condDim = 2**int("kd1" in procoptl) * 3**int("kd2" in procoptl) * 5**int("kd3" in procoptl)
               if condDim%2==0 and globalCondDim%2!=0: globalCondDim = globalCondDim*2
               if condDim%3==0 and globalCondDim%3!=0: globalCondDim = globalCondDim*3
               if condDim%5==0 and globalCondDim%5!=0: globalCondDim = globalCondDim*5
      if globalCondDim==1:
         globalCondDim=0
      else:
         print("Some processes have conditional templates, so globalCondDim={}".format(globalCondDim))

      for proc in self.theInputCard.channels:
         procname = proc[0]
         proctype = proc[3]
         procopts = proc[4]
         isConditional = False
         isDataDriven = ("zjets" in procname.lower() or "zx" in procname.lower())
         removeIntrinsicLumi = False
         procnamefile = procname
         condNormVars=None # Unconditional variables which are not involved in conditional template construction
         condDim=0
         for procopt in procopts:
            procoptl=procopt.lower()
            if "conditional" in procoptl:
               isConditional=True
               condNormVars = ROOT.RooArgSet()
               if not("kd1" in procoptl) and self.KD1 is not None: condNormVars.add(self.KD1)
               if not("kd2" in procoptl) and self.KD2 is not None: condNormVars.add(self.KD2)
               if not("kd3" in procoptl) and self.KD3 is not None: condNormVars.add(self.KD3)
               condDim = 2**int("kd1" in procoptl) * 3**int("kd2" in procoptl) * 5**int("kd3" in procoptl)
               if condNormVars.getSize()==0:
                  raise RuntimeError("Process {} has all variables listed as conditional.".format(procname))
               elif condNormVars.getSize()==self.dataVars.getSize():
                  raise RuntimeError("Process {} has no variables listed as conditional even though the conditional option is specified.".format(procname))
            if "filenamealias" in procoptl:
               procnamefile = procopt.split('=')[1]
            if "datadriven" in procoptl:
               isDataDriven = True
            if "includeslumi" in procoptl:
               removeIntrinsicLumi = True

         bunchNominal = None
         bunchVariations = []
         procPdf = None # Using VerticalInterpPdf
         procRate = None # Using AsymQuad
         procNorm = None # procRate*theLumi

         quadNVars = []

         # Template file name core piece
         templateFileNameMain = "Hto{}_{}_FinalTemplates_".format(self.theChannelName, self.catName)

         # Nominal pdf and rate
         systName = "Nominal"
         templateFileName = "{0}/{1}{2}_{3}{4}".format(self.templateDir, templateFileNameMain, procnamefile, systName, ".root")
         if proctype==0 or isConditional:
            bunchNominal = SimpleTemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,proc,templateFileName,self.iCat,systName)
            bunchNominal.getTemplates()
         else:
            bunchNominal = BSITemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,proc,templateFileName,self.iCat,systName)
            bunchNominal.getTemplates(processName=procname)
         self.theBunches.append(bunchNominal)

         # Systematic "template" variations
         for syst in self.theInputCard.systematics:
            systType = syst[1]
            systOpts = syst[3]
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
                        templateFileName = "{0}/{1}{2}_{3}{4}".format(self.templateDir, templateFileNameMain, procnamefile, systName, ".root")
                        bunchVar = None
                        if proctype==0 or isConditional:
                           bunchVar = SimpleTemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,proc,templateFileName,self.iCat,systName)
                           bunchVar.getTemplates()
                        else:
                           bunchVar = BSITemplateHelper(self.options,self,self.theEqnsMaker,self.theCategorizer,proc,templateFileName,self.iCat,systName)
                           bunchVar.getTemplates(processName=procname)
                        tmplist.append(bunchVar)
                        self.theBunches.append(bunchVar)
                     if len(tmplist)==2: # Expect exactly two templates
                        if self.theSystVarsDict[systVariableName] is not None:
                           tmplist.append(self.theSystVarsDict[systVariableName]) # Append the systematic control RooRealVar as the last element
                           # Check for possible options and append to the list
                           tmpSystOptList=[]
                           for tmpoptpair in systOpts:
                              if procname in tmpoptpair[1] or "all" in tmpoptpair[1].lower():
                                 tmpSystOptList.append(tmpoptpair[0])
                           tmplist.append(tmpSystOptList)
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
                           systName = "{0}UpConst_{1}_{2}_{3}_{4}".format(systVariableName,procname,self.catName,self.theChannelName,self.theSqrtsPeriod)
                        else:
                           systName = "{0}DownConst_{1}_{2}_{3}_{4}".format(systVariableName,procname,self.catName,self.theChannelName,self.theSqrtsPeriod)
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
            procRateExtra = ROOT.AsymQuad(ratename, ratename, extraMorphRateList, extraMorphVarList, 1.0, 2)
            self.extraRateList.append(procRateExtra)

         # Construct the ultimate pdf and rate for the process
         if len(bunchVariations)>0:
            print("Bunch variations are found. Constructing VerticalInterpPdf pdf and AsymQuad rate...")
            morphPdfVarList = ROOT.RooArgList()
            morphPdfList = ROOT.RooArgList()
            morphRateVarList = ROOT.RooArgList()
            morphRateList = ROOT.RooArgList()
            morphPdfList.add(bunchNominal.getThePdf())
            morphRateList.add(bunchNominal.getTheRate())
            for systvar in bunchVariations:
               tmpBunchOpts = systvar[3]
               print("\t- Including bunch variation {}".format(systvar[2].GetName()))
               normOnly=False
               shapeOnly=False
               for tmpBunchOpt in tmpBunchOpts:
                  if "normonly" in tmpBunchOpt.lower():
                     normOnly=True
                  if "shapeonly" in tmpBunchOpt.lower():
                     shapeOnly=True
               if normOnly and shapeOnly:
                  raise RuntimeError("\t=> {} systematic in process {} cannot be both norm-only and shape-only.".format(systvar[2].GetName(), procname))
               elif normOnly:
                  print("\t=> {} is a norm-only systematic in process {}.".format(systvar[2].GetName(), procname))
               elif shapeOnly:
                  print("\t=> {} is a shape-only systematic in process {}.".format(systvar[2].GetName(), procname))
               for isyst in range(0,2):
                  if not normOnly:
                     morphPdfList.add(systvar[isyst].getThePdf())
                  if not shapeOnly:
                     morphRateList.add(systvar[isyst].getTheRate())
               if not normOnly:
                  morphPdfVarList.add(systvar[2])
               if not shapeOnly:
                  morphRateVarList.add(systvar[2])
            procPdf = ROOT.VerticalInterpPdf(procname, procname, morphPdfList, morphPdfVarList,1.0, 2)

            ratename = bunchNominal.getTheRate().GetName() + "_AsymQuad"
            ratename = ratename.replace("_Nominal","")
            procRate = ROOT.AsymQuad(ratename, ratename, morphRateList, morphRateVarList, 1.0, 2)
         else:
            print("Bunch variations do not exist. Constructing pdf and rate from nominal bunch...")
            procPdf = bunchNominal.getThePdf()
            procPdf.SetName(procname)
            procPdf.SetTitle(procname)
            procRate = bunchNominal.getTheRate()

         # Construct the product pdfs from conditional pdf x mass shape
         procExtPdf = None
         if isConditional:
            if self.hasExtMassShapes:
               condpdfname = bunchNominal.getThePdf().GetName() + "_ConditionalPdf"
               condpdfname = condpdfname.replace("_Nominal","")
               procCondPdf = procPdf; procPdf=None
               procCondPdf.SetName(condpdfname); procCondPdf.SetTitle(condpdfname)
               self.extraPdfList.append(procCondPdf)
               procExtPdf = self.extShapeHandle.getThePdf(proc)
               self.extraPdfList.append(procExtPdf)
               procPdf = ROOT.RooProdPdf(
                  procname,procname,
                  ROOT.RooArgSet(procExtPdf),
                  ROOT.RooFit.Conditional(ROOT.RooArgSet(procCondPdf),condNormVars)
               )
            else:
               raise RuntimeError("No external shape handle exists even though process {} is conditional".format(procname))

         self.pdfList.append(procPdf)
         self.rateList.append(procRate)

         normname = procname + "_norm"
         if procRateExtra is None:
            if isDataDriven:
               procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0,1e-15)", ROOT.RooArgList(procRate))
            elif removeIntrinsicLumi:
               procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1/{},1e-15)".format(FloatToString(self.theInputCard.lumi)), ROOT.RooArgList(procRate, self.theLumi))
            else:
               procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1,1e-15)", ROOT.RooArgList(procRate, self.theLumi))
         else:
            if isDataDriven:
               procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1,1e-15)", ROOT.RooArgList(procRate, procRateExtra))
            elif removeIntrinsicLumi:
               procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1*@2/{},1e-15)".format(FloatToString(self.theInputCard.lumi)), ROOT.RooArgList(procRate, procRateExtra, self.theLumi))
            else:
               procNorm = ROOT.RooFormulaVar(normname, "TMath::Max(@0*@1*@2,1e-15)", ROOT.RooArgList(procRate, procRateExtra, self.theLumi))
         self.normList.append(procNorm)

         #if procExtPdf is not None: procExtPdf.Print("v")
         print("Last check on pdf value:",procPdf.getVal(None))
         #if procExtPdf is not None: procExtPdf.Print("v")
         if self.KD1 is not None:
            print("\t- KD1 =",self.KD1.getVal())
         if self.KD2 is not None:
            print("\t- KD2 =",self.KD2.getVal())
         if self.KD3 is not None:
            print("\t- KD3 =",self.KD3.getVal())
         print("Last check on pdf norm:",procNorm.getVal())
         print("\t- Rate = ",procRate.getVal())
         print("\t- Lumi = ",self.theLumi.getVal())
         if procRateExtra is not None:
            print("\t- Extra rate = ",procRateExtra.getVal())

         # Fix conditional dimension binning: For a mass distribution, toy generation might create some issues
         if globalCondDim>0 and not condDimsAdjusted:
            condVars=[]
            if globalCondDim%2==0:
               condVars.append(self.KD1)
            if globalCondDim%3==0:
               condVars.append(self.KD2)
            if globalCondDim%5==0:
               condVars.append(self.KD3)
            for var in condVars:
               print("Adjusting the binning of {}".format(var.GetName()))
               varBinning=var.getBinning()
               varBinArray=varBinning.array()
               varNbins=varBinning.numBins()
               newBinning=ROOT.RooBinning(varBinning)
               ndiv=100
               for ix in range(0,varNbins):
                  lowedge=varBinArray[ix]
                  highedge=varBinArray[ix+1]
                  step=(highedge-lowedge)/float(ndiv)
                  for iy in range(1,ndiv):
                     newBinning.addBoundary(lowedge+step*float(iy))
               var.setBinning(newBinning)
            condDimsAdjusted=True

         # Import the pdf and rates for this process
         getattr(self.workspace, 'import')(procPdf,ROOT.RooFit.RecycleConflictNodes())
         getattr(self.workspace, 'import')(procNorm,ROOT.RooFit.RecycleConflictNodes())

      # Get the data
      self.GetData()

      # Check the pdfs
      print("Now checking the pdfs")
      self.CheckPdfs()

      # Plot the pdfs
      print("Now plotting the pdfs and rates")
      self.PlotPdfs()
      self.PlotRates()

      # Write the workspace
      print("Now writing the workspace")
      self.theWorkspaceFile = ROOT.TFile.Open(self.workspaceFileName,"recreate")
      self.theWorkspaceFile.WriteTObject(self.workspace)
      self.theWorkspaceFile.Close()

      #print("Testing the workspace")
      #ftmp = ROOT.TFile.Open(self.workspaceFileName,"read")
      #wtmp = ftmp.Get("w")
      #wtmp.Print("v")
      #ftmp.Close()

      # Write datacards
      print("Now writing the data cards")
      self.WriteDatacard()

      # Garbage collection
      for bunch in self.theBunches:
         bunch.close()
      if self.extShapeHandle is not None:
         self.extShapeHandle.close()


   def GetData(self):
      data_obs = ROOT.RooDataSet("dummy", "dummy", self.dataVars)
      self.theDataFile = ROOT.TFile.Open(self.dataFileName,"read")
      if not self.theDataFile:
         print("Data tree is not found!")
         self.theDataRDS = data_obs
      else:
         self.theDataTree = self.theDataFile.Get(self.dataTreeName)
         if not self.theDataTree:
            print("File \"",self.dataFileName,"\" or tree \"",self.dataTreeName, "\" is not found.")
            self.theDataRDS = data_obs
         else:
            del(data_obs)

            print("File \"",self.dataFileName,"\" and tree \"",self.dataTreeName,"\" are found. Extracting the data...")

            dataHasMass=False
            if self.theDataTree.GetBranchStatus("mass"):
               self.dataVars.add(self.mass) # If mass is already present in list of data variables, this line does nothing.
               dataHasMass=True

            datasetName = "data_obs_full"
            data_obs = ROOT.RooDataSet(datasetName, datasetName, self.dataVars)

            dataScalars = dict()
            for coord in self.coordList:
               dataScalars[coord] = array('f', [0])
            if dataHasMass:
               dataScalars["mass"] = array('f', [0])
            for name,var in dataScalars.iteritems():
               self.theDataTree.SetBranchAddress(name,var)
            for ev in range(0,self.theDataTree.GetEntries()):
               self.theDataTree.GetEntry(ev)
               for name,var in dataScalars.iteritems():
                  print("\t- Setting variable {} to {} in event {}".format(name,var[0],ev))
                  self.theEqnsMaker.rrvars[name].setVal(var[0])
               data_obs.add(self.dataVars)

            if dataHasMass:
               data_obs_rds = data_obs.reduce("{0}>={1} && {0}<{2}".format(self.mass.GetName(), FloatToString(self.mLow), FloatToString(self.mHigh)))
            else:
               data_obs_rds = data_obs
               print("Data does not have ",self.mass.GetName(),"as a variable.")
            self.theDataRDS = data_obs_rds
      self.theDataRDS.SetName("data_obs")
      self.theDataRDS.Print("v")
      getattr(self.workspace, 'import')(self.theDataRDS, ROOT.RooFit.Rename("data_obs"))


   def WriteDatacard(self):
      print("Being WriteDatacard...")

      self.theDatacardFile = open(self.datacardName, "wb")
      tmplist = self.workspaceFileName.split('/')
      nameWS = tmplist[len(tmplist)-1]

      nSigProcs = self.theInputCard.getNSigProcs()
      nBkgProcs = self.theInputCard.getNBkgProcs()

      self.theDatacardFile.write("##########################################################################################################################\n")
      self.theDatacardFile.write("## Combine manual:                       https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit     ##\n")
      self.theDatacardFile.write("## Non-standard combine use cases:       https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/SWGuideNonStandardCombineUses ##\n")
      self.theDatacardFile.write("## Latest Higgs combination conventions: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions   ##\n")
      self.theDatacardFile.write("## CMS StatCom twiki:                    https://twiki.cern.ch/twiki/bin/view/CMS/StatisticsCommittee                   ##\n")
      self.theDatacardFile.write("##########################################################################################################################\n")
      self.theDatacardFile.write("## Mass window [{0}, {1}]\n".format(FloatToString(self.mLow), FloatToString(self.mHigh)))
      if len(self.normList)==len(self.pdfList):
         for tmppdf, tmpnorm in zip(self.pdfList,self.normList):
            theProcExtRate=float(1)
            procFound=False
            for proc in self.theInputCard.channels:
               if proc[0]==tmppdf.GetName():
                  procFound=True
                  deflumi = proc[2]
                  defrate = proc[1]
                  if deflumi>0.:
                     defrate = defrate/deflumi*self.theInputCard.lumi
                  theProcExtRate = defrate
            if procFound:
               self.theDatacardFile.write("## {0} rate: {1} ({2} events @ {3} fb-1)\n".format(
                     tmppdf.GetName(), FloatToString(float(tmpnorm.getVal()*theProcExtRate/self.theLumi.getVal())),
                     FloatToString(float(tmpnorm.getVal()*theProcExtRate)), FloatToString(self.theLumi.getVal())
                  )
               )
         self.theDatacardFile.write("##########################################################################################################################\n")

      self.theDatacardFile.write("imax *\n")
      self.theDatacardFile.write("jmax *\n")
      self.theDatacardFile.write("kmax *\n")

      self.theDatacardFile.write("------------\n")
      self.theDatacardFile.write("shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
      self.theDatacardFile.write("------------\n")

      binname = "a{0:.0f}".format(self.channel)

      self.theDatacardFile.write("bin {} \n".format(binname))
      self.theDatacardFile.write("observation {0:.0f} \n".format(int(self.theDataRDS.numEntries())))

      self.theDatacardFile.write("------------\n")

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
            ctr_bkg += 1
            self.theDatacardFile.write("{0:.0f} ".format(ctr_bkg))
         else:
            ctr_sig += 1
            self.theDatacardFile.write("{0:.0f} ".format(ctr_sig))
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


   def CheckPdfs(self):
      if self.options.checkpdfs is None: return
      if not self.options.dopdfproj: return
      for tmppdf in self.pdfList:
         for pdfname in self.options.checkpdfs:
            if pdfname in tmppdf.GetName():
               CheckPdf(tmppdf,self.KD1,self.KD2,self.KD3,self.theEqnsMaker.rrvars["fai1"],self.plotsPathName)
               CheckPdf(tmppdf,self.KD1,self.KD2,self.KD3,self.theEqnsMaker.rrvars["GHratio"],self.plotsPathName)
               break


   def PlotPdfs(self):
      if not self.options.dopdfproj: return
      if len(self.normList)==len(self.pdfList):
         for tmppdf, tmpnorm in zip(self.pdfList,self.normList):
            theProcExtRate=float(1)
            procFound=False
            for proc in self.theInputCard.channels:
               if proc[0]==tmppdf.GetName():
                  procFound=True
                  deflumi = proc[2]
                  defrate = proc[1]
                  if deflumi>0.:
                     defrate = defrate/deflumi*self.theInputCard.lumi
                  theProcExtRate = defrate
            if procFound:
               tmprate = tmpnorm.getVal()*theProcExtRate/self.theLumi.getVal()
               PlotPdf1D(tmppdf,tmprate,self.KD1,self.KD2,self.KD3,None,self.plotsPathName)
               args = self.GetProcessSystVars(tmpnorm.GetName())
               for arg in args:
                  PlotPdf1D(tmppdf,tmprate,self.KD1,self.KD2,self.KD3,arg,self.plotsPathName)


   def PlotRates(self):
      for tmpnorm in self.normList:
         isBkg=False
         for proc in self.theInputCard.channels:
            if proc[0] in tmpnorm.GetName():
               procFound=True
               isBkg = (proc[3]==0)
         args = self.GetProcessSystVars(tmpnorm.GetName())
         if not isBkg:
            args.append(self.theEqnsMaker.rrvars["fai1"])
            args.append(self.theEqnsMaker.rrvars["GHratio"])
         if len(args)>0:
            PlotRate(tmpnorm,args,self.plotsPathName)


   def GetProcessSystVars(self,procidname):
      args = []
      for name,systvar in self.theSystVarsDict.iteritems():
         if systvar.hasClients():
            clientsIter = systvar.clientIterator()
            client=clientsIter.Next()
            while client:
               procidnameALT=procidname
               if "_norm" in procidname:
                  procidnameALT=procidname.replace("_norm", "")+"TotalRate"
               if client.GetName() in procidname or ("_norm" in procidname and procidnameALT in client.GetName()):
                  args.append(systvar)
                  break
               client=clientsIter.Next()
      return args



