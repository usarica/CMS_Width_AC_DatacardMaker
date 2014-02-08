#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT

#def transformGraphtoHisto():
print "start"
nBins = 100
startPoint = 120
endPoint = 1600

fileName = ["pdfUncertaintiesDown.root","pdfUncertaintiesUp.root"]
graphName = ["Graph","Graph"]
histName = ["HistDown","HistUp"]
hist = []
fitter = TF1("fitter","pol2")
fout = TFile("fout.root","RECREATE")

for j in range(len(fileName)):
    
    graph = TFile.Open(fileName[j]).Get(graphName[j])
    hist.append(TH1F(histName[j],histName[j],nBins,startPoint,endPoint))

    for i in range (0,nBins):
        #print i
        if hist[j].GetBinLowEdge(i+1)+hist[j].GetBinWidth(i+1) < 1000:
            fitter.SetRange(hist[j].GetBinLowEdge(i+1),hist[j].GetBinLowEdge(i+1)+hist[j].GetBinWidth(i+1))
            graph.Fit(fitter,"RN")
            hist[j].SetBinContent(i+1,fitter.Integral(hist[j].GetBinLowEdge(i+1),hist[j].GetBinLowEdge(i+1)+hist[j].GetBinWidth(i+1)))
        else : hist[j].SetBinContent(i+1,0.0000001)
            
    fout.cd()
    hist[j].Write()
fout.Close()
    
