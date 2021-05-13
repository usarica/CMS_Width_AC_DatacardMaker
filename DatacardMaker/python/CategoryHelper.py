#! /usr/bin/env python


class CategoryHelper:
   def __init__(self, iCatScheme):
      self._catDict={
         "inclusive": 0,
         "vbfcat": 1,
         "vbfvhcat": 2,
         "nj012_2l2nu": 3,
         "nj012boostedhadvh_2l2nu": 4,
         "run2legacy4l": 5
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
      elif(self.iCatScheme == self._catDict["run2legacy4l"]):
         self.catNameList.append("VHLepttagged")
         self.catNameList.append("VHHadrtagged")
         self.catNameList.append("VBFtagged")
         self.catNameList.append("VBF1jtagged")
         self.catNameList.append("Boosted")
         self.catNameList.append("Untagged")
      elif(self.iCatScheme == self._catDict["nj012_2l2nu"]):
         self.catNameList.append("Nj_eq_0")
         self.catNameList.append("Nj_eq_1")
         self.catNameList.append("Nj_geq_2_pTmiss_lt_200")
         self.catNameList.append("Nj_geq_2_pTmiss_ge_200")
      elif(self.iCatScheme == self._catDict["nj012_boostedhadvh_2l2nu"]):
         self.catNameList.append("Nj_eq_0")
         self.catNameList.append("Nj_eq_1")
         self.catNameList.append("Nj_geq_2_pTmiss_lt_200")
         self.catNameList.append("Nj_geq_2_pTmiss_ge_200")
         self.catNameList.append("BoostedHadVH")
      else: # Undefined categorization
         errormsg="CategoryHelper::init: Categorization scheme {} is not defined. The following enumerators need to be used:".format(self.iCatScheme)
         for key,val in self._catDict.iteritems():
            errormsg = errormsg + "\n - {}: {}".format(key,val)
         raise RuntimeError(errormsg)
      self.nCategories = len(self.catNameList)
