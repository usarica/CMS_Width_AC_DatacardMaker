#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT

folder1 = "/afs/cern.ch/work/u/usarica/public/CombineTemplates/8TeV/"
namefileNew = "HtoZZ4l_gg2VV_125p6_ModifiedTemplatesForCombine_D_Gamma_gg_r10_Nominal.root"
tnames = ["T_2D_1","T_2D_2","T_2D_4","","T_2D_qqZZ"]
onames = ["signal","background","interference","total","qqZZ_background"]
channels = ["2e2mu","4e","4mu"]
for i in range(len(channels)):
    fileNew = TFile.Open(folder1+channels[i]+"/220/"+namefileNew)
    print " "
    print channels[i]
    for j in range(len(onames)):
        if j is not 3 :
            print onames[j], ": %.3f" % (fileNew.Get(tnames[j]).Integral("width")*19.712)
        else :
            print onames[j], ": %.3f" % ((fileNew.Get(tnames[0]).Integral("width")+fileNew.Get(tnames[1]).Integral("width")+fileNew.Get(tnames[2]).Integral("width"))*19.712)

          
    
