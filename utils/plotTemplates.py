#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT

#file = TFile.Open("/afs/cern.ch/work/u/usarica/public/CombineTemplates/2e2mu/220/HtoZZ4l_gg2VV_125p6_TemplatesForCombine_D_Gamma_gg_r10_Nominal.root")
#h = file.Get("T_2D_1").ProjectionX()
#h.DrawClone()
#h = file.Get("T_2D_2").ProjectionX()
#h.DrawClone("SAME")
#h = file.Get("T_2D_4").ProjectionX()
#h.DrawClone("SAME")


channels = ["2e2mu","4e","4mu"]
namefileOld = "HtoZZ4l_gg2VV_125p6_TemplatesForCombine_D_Gamma_gg_r10_Nominal.root"
namefileNew = "HtoZZ4l_gg2VV_125p6_ModifiedTemplatesForCombine_D_Gamma_gg_r10_Nominal.root"

folder1 = "/afs/cern.ch/work/u/usarica/public/CombineTemplates/8TeV/"
folder2 = "/afs/cern.ch/work/u/usarica/public/CombineTemplates/forGiacomo_OldTemplates/WidthTemplates/8TeV/"

tnames = ["T_2D_1","T_2D_2","T_2D_4","T_2D_qqZZ"]

canvases = []
for i in range(len(tnames)):
#for i in range(1):
    #print i
    canvases.append(TCanvas(tnames[i],tnames[i],500,500))
    canvases[i].Divide(2,1)
    holds = []
    hnews = []
    holdsk = []
    hnewsk = []
    for j in range(len(channels)):
    #for j in range(1):
        fileOld = TFile.Open(folder2+channels[j]+"/220/"+namefileOld)
        fileNew = TFile.Open(folder1+channels[j]+"/220/"+namefileNew)
        nameToLoad = tnames[i]
        if i == 3 : #load Roberto for old
            #nameToLoad = "mZZ_qq"
            file = TFile.Open("templates2D/templ2D_"+channels[j]+"_8TeV_m4l.root")
            holds.append(file.Get("mZZ_qq").ProjectionX(channels[j]+"h0"))
        else : holds.append(fileOld.Get(nameToLoad).ProjectionX(channels[j]+"h0"))
        hnews.append(fileNew.Get(nameToLoad).ProjectionX(channels[j]+"hN"))
        #print j, " ", len(holds)
        print folder1+channels[j]+"/220/"+namefileNew
        holds[j].SetLineColor(j+1)
        hnews[j].SetLineColor(j+1)
        hnews[j].SetLineStyle(2)
        canvases[i].cd(1)
        samestr = "E SAME"
        if j == 0 : samestr = "E"
        hnews[j].Draw(samestr)
        #holds[j].Draw("SAME")

        
        holdsk.append(fileOld.Get(nameToLoad).ProjectionY(channels[j]+"k0"))
        hnewsk.append(fileNew.Get(nameToLoad).ProjectionY(channels[j]+"kn"))
        holdsk[j].SetLineColor(j+1)
        hnewsk[j].SetLineColor(j+1)
        hnewsk[j].SetLineStyle(2)
        canvases[i].cd(2)
        hnewsk[j].Draw(samestr)
        holdsk[j].Draw("SAME")
    canvases[i].SaveAs(tnames[i]+".png")
            

            
        
