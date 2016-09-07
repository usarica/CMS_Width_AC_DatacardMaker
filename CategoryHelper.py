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

   def __init__(self, iCatScheme, workspace):
      self.iCatScheme = iCatScheme
      self.catNameList = []

      if(self.iCatScheme == 1): # icat==0: VBF, ==1: Non-VBF
         self.catNameList.append("Djet")
         self.catNameList.append("nonDjet")
      else: # No categorization
         self.catNameList.append("")
      self.nCategories = len(self.catNameList)

      self.systVar_catsplit_gg = []
      self.systVar_catsplit_VBF = []

      for icat in range(0,self.nCategories-1):
        self.systVar_catsplit_gg.append(workspace.factory("Djetscale_ggzz[-3,3]"))
        self.systVar_catsplit_VBF.append(workspace.factory("Djetscale_vbf_offshell[-3,3]"))
        self.systVar_catsplit_qqzz.append(workspace.factory("Djetscale_qqzz[-3,3]"))
        self.systVar_catsplit_zjets.append(workspace.factory("Djetscale_zjets[-3,3]"))
        # No neec to fix these since nCategories==1 does not enter the loop (i.e. no systVars are created)!






