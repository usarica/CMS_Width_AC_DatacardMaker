#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array

class CategoryHelper:

   def __init__(
               self,
               iCatScheme#,
               #workspace
               ):
      self.iCatScheme = iCatScheme
      self.catNameList = []

      if(self.iCatScheme == 1): # icat==0: VBF, ==1: Non-VBF
         self.catNameList.append("Djet")
         self.catNameList.append("nonDjet")
      elif(self.iCatScheme == 2): # icat==0: VBF, ==1: VH, ==3: Untagged
         self.catNameList.append("VBFTagged")
         self.catNameList.append("VHTagged")
         self.catNameList.append("Untagged")
      else: # No categorization
         self.catNameList.append("")
      self.nCategories = len(self.catNameList)

      #NO NEED FOR THESE, THE CONSTRUCTION OF ASYMQUAD AND VERTICALHISTINTERPOLATORPDF HANDLE EVERYTHING CONSISTENTLY
      #self.CatUncScheme=0 # HARDCODED, TO BE DECIDED LATER ON
      #self.systVar_catsplit_gg = []
      #self.systVar_catsplit_VBF = []
      #self.systVar_catsplit_qqzz = []
      #self.systVar_catsplit_zjets = []

      ## No need to fix these since nCategories==1 does not enter the loop (i.e. no systVars are created)!
      #if(self.CatUncScheme==1): # Uncertainties on i->j contamination
         #for icat in range(0,self.nCategories):
            #for jcat in range(0,self.nCategories):
               #if(jcat<=icat):
                  #continue
               #self.systVar_catsplit_gg.append(workspace.factory("cat_{0:.0f}in{1:.0f}_syst_ggzz[-3,3]".format(icat,jcat)))
               #self.systVar_catsplit_VBF.append(workspace.factory("cat_{0:.0f}in{1:.0f}_syst_vbf_offshell[-3,3]".format(icat,jcat)))
               #self.systVar_catsplit_qqzz.append(workspace.factory("cat_{0:.0f}in{1:.0f}_syst_qqzz[-3,3]".format(icat,jcat)))
               #self.systVar_catsplit_zjets.append(workspace.factory("cat_{0:.0f}in{1:.0f}_syst_zjets[-3,3]".format(icat,jcat)))
      #else: # Uncertainties on i->un-tagged contamination
         #for icat in range(0,self.nCategories-1):
            #self.systVar_catsplit_gg.append(workspace.factory("cat_{0:.0f}_syst_ggzz[-3,3]".format(icat)))
            #self.systVar_catsplit_VBF.append(workspace.factory("cat_{0:.0f}_syst_vbf_offshell[-3,3]".format(icat)))
            #self.systVar_catsplit_qqzz.append(workspace.factory("cat_{0:.0f}_syst_qqzz[-3,3]".format(icat)))
            #self.systVar_catsplit_zjets.append(workspace.factory("cat_{0:.0f}_syst_zjets[-3,3]".format(icat)))





