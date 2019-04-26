#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from array import array

class CategoryHelper:
   def __init__(self, iCatScheme):
      self.iCatScheme = iCatScheme
      self.catNameList = []

   # NOTE: Add as many categorization schemes as necessary
      if(self.iCatScheme == 0): # No categorization
         self.catNameList.append("Inclusive")
      elif(self.iCatScheme == 1): # icat==0: VBF, ==1: Untagged
         self.catNameList.append("JJVBFTagged")
         self.catNameList.append("Untagged")
      elif(self.iCatScheme == 2): # icat==0: VBF, ==1: VH, ==3: Untagged
         self.catNameList.append("JJVBFTagged")
         self.catNameList.append("HadVHTagged")
         self.catNameList.append("Untagged")
      else: # No categorization
         raise RuntimeError(
            "CategoryHelper::init: Categorization scheme {} is not defined.\n".format(self.iCatScheme)
            +
            " - If you want to use the inclusive category, pass iCatScheme=0."
         )
      self.nCategories = len(self.catNameList)



