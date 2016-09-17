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



