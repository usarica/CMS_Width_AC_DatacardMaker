#! /usr/bin/env python


class CategoryHelper:
   def __init__(self, iCatScheme):
      self._catDict={
         "inclusive": 0,
         "vbfcat": 1,
         "vbfvhcat": 2
      }

      try:
         iCatScheme = int(iCatScheme)
      except ValueError:
         if iCatScheme in self._catDict.keys():
            iCatScheme = self._catDict[iCatScheme]

      self.iCatScheme = iCatScheme
      self.catNameList = []

   # NOTE: Add as many categorization schemes as necessary
      if(self.iCatScheme == self._catDict["inclusive"]): # No categorization
         self.catNameList.append("Inclusive")
      elif(self.iCatScheme == self._catDict["vbfcat"]): # icat==0: VBF, ==1: Untagged
         self.catNameList.append("JJVBFTagged")
         self.catNameList.append("Untagged")
      elif(self.iCatScheme == self._catDict["vbfvhcat"]): # icat==0: VBF, ==1: VH, ==2: Untagged
         self.catNameList.append("JJVBFTagged")
         self.catNameList.append("HadVHTagged")
         self.catNameList.append("Untagged")
      else: # Undefined categorization
         errormsg="CategoryHelper::init: Categorization scheme {} is not defined. The following enumerators need to be used:".format(self.iCatScheme)
         for key,val in self._catDict.iteritems():
            errormsg = errormsg + "\n - {}: {}".format(key,val)
         raise RuntimeError(errormsg)
      self.nCategories = len(self.catNameList)



