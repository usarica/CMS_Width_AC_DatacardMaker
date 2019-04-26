#! /usr/bin/env python
import os
import re
import math
import ROOT
from array import array

# Some helper functions that don't need to belong to the systematics class
def FloatToString(inputValue):
   return ('%.10f' % inputValue).rstrip('0').rstrip('.')
def FloatToStringScientific(inputValue,etype='e'):
   s = FloatToString(inputValue)
   q = s.replace(".","").lstrip('0')
   f = float(s)
   nq = len(q)-1
   strcmd = ("%%%c%i%c" % ('.',nq,etype))
   return (strcmd % f)
def GetDataPeriodString(sqrts,period):
   try: # Check whether python knows about 'basestring'
      basestring
   except NameError: # No, it doesn't (it's Python3); use 'str' instead
      basestring=str
   res=""
   if sqrts is None:
      raise RuntimeError("sqrts has to have a value")
   elif isinstance(sqrts,basestring):
      res = "{}TeV".format(sqrts)
   else:
      res = "{}TeV".format(FloatToString(sqrts))
   if period is not None:
      if isinstance(period,basestring):
         res = "{}_{}".format(res,period)
      else:
         res = "{}_{}".format(res,FloatToString(period))
   return res

