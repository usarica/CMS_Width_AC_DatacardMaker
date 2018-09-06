#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSpline.h"
#include "TLegend.h"
#include "TAxis.h"


TGraph* bestFit(TTree* t, TString x, TString y){
  TString strDraw = y + ":" + x;
  cout << "Drawing " << strDraw << endl;
  cout << "Number of entries at minimum: " << t->GetEntries("deltaNLL==0") << endl;

  typedef float var_t;
  var_t xvar, yvar, deltaNLL=9999;
  t->ResetBranchAddresses();
  t->SetBranchAddress(x, &xvar);
  t->SetBranchAddress(y, &yvar);
  t->SetBranchAddress("deltaNLL", &deltaNLL);
  double xx[1]={ 0 }; double yy[1]={ 0 };
  for (int ev=0; ev<t->GetEntries(); ev++){
    t->GetEntry(ev);
    if (deltaNLL==0.){
      xx[0] = xvar;
      yy[0] = yvar;
      break;
    }
  }
  TGraph* gr0 = new TGraph(1, xx, yy);
  if (!gr0) return nullptr;
  else cout << "Best fit found, attempting..." << endl;
  gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
  cout << "Marker style of best fit is set" << endl;
  if (gr0->GetN() > 1) gr0->Set(1);
  cout << "Set of best fit is set" << endl;
  return gr0;
}


int countGridPointsFromTree(TTree *t, TString x, double xmin = -1, double xmax = -1) {
    if (xmin == xmax) {
        xmin = t->GetMinimum(x);
        xmax = t->GetMaximum(x);
    }
    t->Draw(Form("%s>>h1000(1000,%10g,%10g)", x.Data(),xmin-1e-4,xmax+1e-4), "deltaNLL > 0");
    TH1 *h1000 = (TH1*) gROOT->FindObject("h1000");
    int bins = 0;
    for (int i = 1, n = h1000->GetNbinsX(); i <= n; ++i) {
        if (h1000->GetBinContent(i) != 0) bins++;
    }
    h1000->Delete();
    return bins;
}
TH2* treeToHist2D(TTree *t, TString x, TString y, TString name, double xmin = -1, double xmax = -1, double ymin = -1, double ymax = -1) {
    if (xmin == xmax) {
        xmin = t->GetMinimum(x);
        xmax = t->GetMaximum(x);
    } 
    if (ymin == ymax) {
        ymin = t->GetMinimum(y);
        ymax = t->GetMaximum(y);
    }
    int xbins = countGridPointsFromTree(t,x,xmin,xmax);
    int ybins = countGridPointsFromTree(t,y,ymin,ymax);
    double dx = (xmax-xmin)/(xbins-1);
    double dy = (ymax-ymin)/(ybins-1);
    xmin -= 0.5*dx; xmax += 0.5*dx;
    ymin -= 0.5*dy; ymax += 0.5*dy;
    if (fabs(xmin) < 1e-5) xmin = 0;
    if (fabs(xmax) < 1e-5) xmax = 0;
    //std::cout << "In making " << name << ", guessed " << xbins << " bins for " << x << " from " << xmin << " to " << xmax << std::endl;
    //std::cout << "In making " << name << ", guessed " << ybins << " bins for " << y << " from " << ymin << " to " << ymax << std::endl;
    t->Draw(Form("2*deltaNLL:%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)", y.Data(), x.Data(), name.Data(), xbins, xmin, xmax, ybins, ymin, ymax), "deltaNLL != 0", "PROF");
    TH2 *prof = (TH2*) gROOT->FindObject(name+"_prof");
    TH2F *h2d = new TH2F(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
    for (int ix = 1; ix <= xbins; ++ix) {
        for (int iy = 1; iy <= ybins; ++iy) {
             h2d->SetBinContent(ix, iy, prof->GetBinContent(ix,iy));
        }
    }
    h2d->SetDirectory(0);
    return h2d;
}

TH2F* frameTH2F(TH2F *in){

        double frameValue = 1000;
        if (TString(in->GetName()).Contains("bayes")) frameValue = 0.0;

	Double_t xw = in->GetXaxis()->GetBinWidth(0);
	Double_t yw = in->GetYaxis()->GetBinWidth(0);

	Int_t nx = in->GetNbinsX();
	Int_t ny = in->GetNbinsY();

	Double_t x0 = in->GetXaxis()->GetXmin();
	Double_t x1 = in->GetXaxis()->GetXmax();

	Double_t y0 = in->GetYaxis()->GetXmin();
	Double_t y1 = in->GetYaxis()->GetXmax();

	TH2F *framed = new TH2F(
			Form("%s framed",in->GetName()),
			Form("%s framed",in->GetTitle()),
			nx + 2, x0-xw, x1+xw,
			ny + 2, y0-yw, y1+yw
			);

	//Copy over the contents
	for(int ix = 1; ix <= nx ; ix++){
		for(int iy = 1; iy <= ny ; iy++){
			framed->SetBinContent(1+ix, 1+iy, in->GetBinContent(ix,iy));
		}
	}
	//Frame with huge values
	nx = framed->GetNbinsX();
	ny = framed->GetNbinsY();
	for(int ix = 1; ix <= nx ; ix++){
		framed->SetBinContent(ix,  1, frameValue);
		framed->SetBinContent(ix, ny, frameValue);
	}
	for(int iy = 2; iy <= ny-1 ; iy++){
		framed->SetBinContent( 1, iy, frameValue);
		framed->SetBinContent(nx, iy, frameValue);
	}

	return framed;
}



TList* contourFromTH2(TH2 *h2in, double threshold) {
    std::cout << "Getting contour at threshold " << threshold << " from " << h2in->GetName() << std::endl;
    //http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
    Double_t contours[1];
    contours[0] = threshold;

    TH2F *h2 = frameTH2F((TH2F*)h2in);

    h2->SetContour(1, contours);

    // Draw contours as filled regions, and Save points
    h2->Draw("CONT Z LIST");
    gPad->Update(); // Needed to force the plotting and retrieve the contours in TGraphs

    // Get Contours
    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;

    if (conts == NULL || conts->GetSize() == 0){
        printf("*** No Contours Were Extracted!\n");
        return 0;
    }

    TList *ret = new TList();
    for(int i = 0; i < conts->GetSize(); i++){
        contLevel = (TList*)conts->At(i);
        printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
        for (int j = 0, n = contLevel->GetSize(); j < n; ++j) {
            TGraph *gr1 = (TGraph*) contLevel->At(j)->Clone();
            ret->Add(gr1);
            //break;
        }
    }
    return ret;
}

