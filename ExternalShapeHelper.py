#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
import ROOT
from array import array
from CategoryHelper import CategoryHelper
from ExtendedTemplate import ExtendedTemplate

class ExternalShapeHelper:
   def __init__(self, options, theMaker, theEqnsMaker, theCategorizer, shapesFileName, iCat):
      self.condDim = 0
      # sqrts and channel index from the datacard maker class
      self.sqrts = theMaker.sqrts
      self.theSqrtsPeriod = theMaker.theSqrtsPeriod
      self.channel = theMaker.channel
      self.theChannelName = theMaker.theChannelName
      self.theInputCard = theMaker.theInputCard

      # RooRealVars from the datacard maker class
      self.mass = theMaker.mass
      self.MH = theEqnsMaker.rrvars["MH"]

      self.iCatScheme = theCategorizer.iCatScheme
      self.catNameList = theCategorizer.catNameList
      self.nCategories = theCategorizer.nCategories
      self.iCat = iCat
      if self.iCat>=self.nCategories:
         sys.exit("self.iCat={} >= self.nCategories={}!".format(self.iCat,self.nCategories))

      self.shapesFileName = shapesFileName
      self.shapeSuffix = "{0}_{1}_{2}".format(self.catNameList[self.iCat],self.theChannelName,self.theSqrtsPeriod)

      self.shapeFile = None
      self.theWS = None
      self.thePdfs = dict()

      self.getShapes()


# Open the shape files
   def openFile(self):
      print "Opening file ",self.shapesFileName
      self.shapeFile = ROOT.TFile.Open(self.shapesFileName, "read")
      if self.shapeFile is None:
         raise RuntimeError("ExternalShapeHelper file {} is None!".format(self.shapesFileName))
      elif self.shapeFile.IsZombie():
         raise RuntimeError("ExternalShapeHelper could not open file {}!".format(self.shapesFileName))
# Close the shape files
   def close(self):
      if self.shapeFile is not None:
         if self.shapeFile.IsOpen():
            self.shapeFile.Close()


   def getThePdf(self, theProcess):
      procname = theProcess[0]
      if procname in self.thePdfs:
         return self.thePdfs[procname]
      else:
         raise RuntimeError("ExternalShapeHelper: External shape for process {} does not exist".format(procname))
         return None


# Get shapes for each category
   def getShapes(self):
      self.openFile()
      self.theWS=self.shapeFile.Get("w").Clone("WSinput_{}".format(self.shapeSuffix))

      MH_in=self.theWS.factory("MH")
      mass_in=self.theWS.factory("mass")
      w_in_vars = [ mass_in, MH_in ]
      w_out_vars = [ self.mass, self.MH ]

      for var, var_out in zip(w_in_vars,w_out_vars):
         if var.hasClients():
            clientsIter = var.clientIterator()
            client=clientsIter.Next()
            while client:
               client.redirectServers(ROOT.RooArgSet(var_out), False, False, False)
               client=clientsIter.Next()

      for proc in self.theInputCard.channels:
         procname = proc[0]
         proctype = proc[3]
         procopts = proc[4]
         procTplAlias = procname
         isConditional=False
         for procopt in procopts:
            procoptl=procopt.lower()
            if "conditional" in procoptl:
               isConditional=True
            if "templatenamealias" in procoptl:
               procTplAlias = procopt.split('=')[1]
         if isConditional:
            shapeName = "{}_MassShape".format(procTplAlias)
            #shapeName = "MassShapeModel"
            pdf=self.theWS.pdf(shapeName)
            pdf.SetName("{}_{}_ExtMassShape".format(procname, self.shapeSuffix))
            pdf.SetTitle("{}_{}_ExtMassShape".format(procname, self.shapeSuffix))
            print "Acquiring mass shape {} for process {}".format(shapeName, procname)
            self.thePdfs[procname]=pdf

