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
#include "TMarker.h"
#include "contours.h"


using namespace std;


template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it>=val){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
    }
  }
  if (!inserted) valArray.push_back(val);
}

std::vector<TString> lsdir(TString const& indir){
  std::vector<TString> res;

  struct dirent* ep;
  DIR* dp = opendir(indir.Data());
  if (dp != NULL){
    while ((ep = readdir(dp))) res.push_back(ep->d_name);
    closedir(dp);
  }
  else cerr << "Couldn't open the directory" << endl;

  return res;
}


void fixHistogram(TH2F*& hist){
  const int nx = hist->GetNbinsX();
  const int ny = hist->GetNbinsY();
  // ix, 1
  for (int ix=2; ix<=nx-1; ix++){
    float bc = hist->GetBinContent(ix, 1);
    if (bc==0.){
      float bcp = hist->GetBinContent(ix-1, 1);
      float bcn = hist->GetBinContent(ix+1, 1);
      bc = (bcp+bcn)/2.;
      cout << "Fixing bin " << ix << " , " << 1 << endl;
      hist->SetBinContent(ix, 1, bc);
    }
  }
  // ix, ny
  for (int ix=2; ix<=nx-1; ix++){
    float bc = hist->GetBinContent(ix, ny);
    if (bc==0.){
      float bcp = hist->GetBinContent(ix-1, ny);
      float bcn = hist->GetBinContent(ix+1, ny);
      bc = (bcp+bcn)/2.;
      cout << "Fixing bin " << ix << " , " << ny << endl;
      hist->SetBinContent(ix, ny, bc);
    }
  }
  // 1, iy
  for (int iy=2; iy<=ny-1; iy++){
    float bc = hist->GetBinContent(1, iy);
    if (bc==0.){
      float bcp = hist->GetBinContent(1, iy-1);
      float bcn = hist->GetBinContent(1, iy+1);
      bc = (bcp+bcn)/2.;
      cout << "Fixing bin " << 1 << " , " << iy << endl;
      hist->SetBinContent(1, iy, bc);
    }
  }
  // nx, iy
  for (int iy=2; iy<=ny-1; iy++){
    float bc = hist->GetBinContent(nx, iy);
    if (bc==0.){
      float bcp = hist->GetBinContent(nx, iy-1);
      float bcn = hist->GetBinContent(nx, iy+1);
      bc = (bcp+bcn)/2.;
      cout << "Fixing bin " << nx << " , " << iy << endl;
      hist->SetBinContent(nx, iy, bc);
    }
  }
  // ix, iy
  for (int ix=2; ix<=nx-1; ix++){
    for (int iy=2; iy<=ny-1; iy++){
      float bc = hist->GetBinContent(ix, iy);
      if (bc==0.){
        float bcpy = hist->GetBinContent(ix, iy-1);
        float bcny = hist->GetBinContent(ix, iy+1);
        float bcpx = hist->GetBinContent(ix-1, iy);
        float bcnx = hist->GetBinContent(ix+1, iy);
        bc = (bcpx+bcnx+bcpy+bcny)/4.;
        cout << "Fixing bin " << ix << " , " << iy << endl;
        hist->SetBinContent(ix, iy, bc);
      }
    }
  }
}
TH2F* getHistogramFromTree(TTree* tree, TString const strxvar, TString const stryvar){
  typedef float var_t;

  TH2F* res=nullptr;
  vector<pair<var_t, var_t>> points;
  var_t xvar, yvar, deltaNLL=9999;
  std::vector<var_t> xvalList, yvalList;

  if (tree){
    tree->SetBranchAddress(strxvar, &xvar);
    tree->SetBranchAddress(stryvar, &yvar);
    tree->SetBranchAddress("deltaNLL", &deltaNLL);
    var_t minNLL=1e9;
    var_t minX, minY;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      deltaNLL = deltaNLL*2.;
      if (deltaNLL<0.) continue;
      minNLL = std::min(minNLL, deltaNLL);
    }
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      var_t xval = xvar;
      var_t yval = yvar;
      if (strxvar=="GGsm") xval *= 4.07;
      if (stryvar=="GGsm") yval *= 4.07;
      deltaNLL = deltaNLL*2.;
      if (minNLL!=deltaNLL){
        bool doAdd=true;
        doAdd &= !((strxvar=="GGsm" && xval>21.) || (stryvar=="GGsm" && yval>21.));
        if (doAdd){
          addByLowest(xvalList, xval, true);
          addByLowest(yvalList, yval, true);
        }
      }
    }
    std::cout << "Minimum deltaNLL = " << minNLL << std::endl;
    var_t xmin=xvalList.front();
    var_t xmax=xvalList.back();
    var_t ymin=yvalList.front();
    var_t ymax=yvalList.back();
    var_t xwidth=(xmax-xmin)/float(xvalList.size()-1);
    var_t ywidth=(ymax-ymin)/float(yvalList.size()-1);
    cout << xvalList.size() << endl;
    cout << yvalList.size() << endl;
    cout << "x = [ " << xmin << " , " << xmax << " , " << xwidth << " ]" << endl;
    cout << "y = [ " << ymin << " , " << ymax << " , " << ywidth << " ]" << endl;
    std::vector<var_t> xBoundList, yBoundList;
    for (unsigned int i=0; i<xvalList.size()-1; i++) addByLowest(xBoundList, float((xvalList.at(i)+xvalList.at(i+1))/2.), true);
    addByLowest(xBoundList, float(xvalList.at(0) - (xvalList.at(1)-xvalList.at(0))/2.), true);
    addByLowest(xBoundList, float(xvalList.at(xvalList.size()-1) + (xvalList.at(xvalList.size()-1)-xvalList.at(xvalList.size()-2))/2.), true);
    for (unsigned int i=0; i<yvalList.size()-1; i++) addByLowest(yBoundList, float((yvalList.at(i)+yvalList.at(i+1))/2.), true);
    addByLowest(yBoundList, float(yvalList.at(0) - (yvalList.at(1)-yvalList.at(0))/2.), true);
    addByLowest(yBoundList, float(yvalList.at(yvalList.size()-1) + (yvalList.at(yvalList.size()-1)-yvalList.at(yvalList.size()-2))/2.), true);
    res = new TH2F(Form("gr_%s_%s", strxvar.Data(), stryvar.Data()), "", xvalList.size(), xBoundList.data(), yvalList.size(), yBoundList.data());
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      var_t xval = xvar;
      var_t yval = yvar;
      if (strxvar=="GGsm") xval *= 4.07;
      if (stryvar=="GGsm") yval *= 4.07;
      //deltaNLL = deltaNLL*2.-minNLL;
      deltaNLL = deltaNLL*2.;
      if (deltaNLL<0.) continue;
      int ix=res->GetXaxis()->FindBin(xval);
      int iy=res->GetYaxis()->FindBin(yval);
      double bc=res->GetBinContent(ix, iy);
      if (bc==0.) res->SetBinContent(ix, iy, deltaNLL);
      //cout << "X,Y,L = " << xvar << "," << yvar << "," << deltaNLL << endl;
    }
    fixHistogram(res);
  }
  return res;
}
TString getVariableLabel(TString const strvar, TString const strachypo){
  if (strvar=="deltaNLL") return "-2 #Delta lnL";
  else if (strvar=="GGsm") return "#Gamma_{H} (MeV)";
  else if (strvar=="CMS_zz4l_fai1"){
    TString strai=strachypo.Data();
    if (strai.Contains("Lambda")) strai.Prepend("#");
    return Form("f_{%s} cos(#phi_{#lower[-0.2]{%s}})", strai.Data(), strai.Data());
  }
  else return strvar;
}

void plotScan2D(TString const indir, TString const strxvar, TString const stryvar, TString const strachypo=""){
  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.20);
  gROOT->ForceStyle();

  TChain* tree=new TChain("limit");
  vector<TString> lsdirs = lsdir("./");
  for (auto const& strdir:lsdirs){
    if (strdir.Contains(indir)){
      vector<TString> lscontent = lsdir(strdir);
      for (auto const& s:lscontent){ if (s.Contains(".root")) tree->Add(strdir+"/"+s); }
    }
  }
  TH2F* hh = getHistogramFromTree(tree, strxvar, stryvar);

  if (!hh) cout << "ERROR: No histogram!" << endl;
  else cout << "Histogram is acquired" << endl;

  hh->GetXaxis()->SetTitle(getVariableLabel(strxvar, strachypo));
  hh->GetYaxis()->SetTitle(getVariableLabel(stryvar, strachypo));

  TFile *ff = new TFile("__tmp__.root","RECREATE");
  ff->cd();
  TList* tlFF68 = contourFromTH2(hh,2.3);
  TList* tlFF95 = contourFromTH2(hh,5.99);
  TGraph* g68=nullptr;
  TGraph* g95=nullptr;
  if (tlFF68 && tlFF68->GetSize()>0) g68 = (TGraph*) tlFF68->At(0);
  if (tlFF95 && tlFF95->GetSize()>0) g95 = (TGraph*) tlFF95->At(0);
  if (g68){
    cout << "68% CL contour is present" << endl;
    g68->SetName("g68");
    g68->SetLineColor(kBlack);
    g68->SetLineWidth(2);
    g68->SetLineStyle(6);
    //  g68->SetFillColor(kAzure+1);
    //  g68->SetFillStyle(3001);
    cout << "68% CL contour is processed" << endl;
  }
  if (g95){
    cout << "95% CL contour is present" << endl;
    g95->SetName("g95");
    g95->SetLineColor(kBlack);
    g95->SetLineWidth(2);
    g95->SetLineStyle(1);
    //  g95->SetFillColor(kCyan);
    //  g95->SetFillStyle(3001);
    cout << "95% CL contour is processed" << endl;
  }

  TGraph* best = bestFit(tree,strxvar,stryvar);
  if (!best) cout << "ERROR: Best fit not found!" << endl;
  else cout << "Best fit is acquired" << endl;
  if (strxvar=="GGsm") best->GetX()[0] *= 4.07;
  best->SetMarkerStyle(5);
  best->SetMarkerSize(1.5);

  TString canvasname = Form("c_%s%s%s_%sVS%s", indir.Data(), (strachypo!="" ? "_" : ""), strachypo.Data(), stryvar.Data(), strxvar.Data());
  TCanvas* c = new TCanvas(canvasname, "", 1000, 800);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->cd();
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);

  double bestx,besty;
  best->GetPoint(0,bestx,besty);

  hh->GetXaxis()->SetLabelSize(0.04);
  hh->GetXaxis()->SetRangeUser(0., 20.);
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetLabelSize(0.04);
  hh->GetYaxis()->SetTitleOffset(1);
  hh->GetYaxis()->CenterTitle();

   
  TLegend *l = new TLegend(0.635,0.58,0.835,0.78);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(42);
  if (g95) l->AddEntry(g95,"95% CL","l");
  if (g68) l->AddEntry(g68,"68% CL","l");
  l->AddEntry(best,"Best Fit","p");
  c->cd();
  hh->Draw("colz");
  if (g95) g95->Draw("Csame");
  if (g68) g68->Draw("Csame");


  TMarker m;
  m.SetMarkerSize(2);
  m.SetMarkerColor(kRed+1);
  m.SetMarkerStyle(33);
  //  m.SetMarkerSize(1.5);
  //  m.SetMarkerColor(kBlack);
  //  m.SetMarkerStyle(5);
  m.DrawMarker(4.07, 0);
  l->AddEntry(&m, "SM", "p");
  best->Draw("Psame");
  l->Draw();

  TText* text=nullptr;
  TPaveText* ptc = new TPaveText(0.635, 0.84, 0.92, 0.94, "brNDC");
  ptc->SetBorderSize(0);
  ptc->SetFillStyle(0);
  ptc->SetTextAlign(12);
  ptc->SetTextFont(42);
  ptc->SetTextSize(0.04);
  text = ptc->AddText(0.02, 0.45, "#font[61]{CMS}");
  //text = ptc->AddText(0.12, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.06);
  ptc->Draw();

	TPaveText* pt = new TPaveText(0.15,0.955,0.92,1,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextAlign(12);
	pt->SetTextFont(42);
	pt->SetTextSize(0.04);
  float lumi2011=5.1;
  float lumi2012=19.7;
  float lumi2015=2.7;
  float lumi201617=77.5;
  //text = pt->AddText(0.68, 0.45, Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi201617+0.05));
  text = pt->AddText(0.2, 0.45, Form("#font[42]{%.1f fb^{-1} (7 TeV) + %.1f fb^{-1} (8 TeV) + %.1f fb^{-1} (13 TeV)}", lumi2011, lumi2012, lumi2015+lumi201617+0.05));
  text->SetTextSize(0.0315);
  pt->Draw();

  TPaveText* pt2 = new TPaveText(0.95, 0.33, 0.99, 0.99, "brNDC");
  //TPaveText* pt2 = new TPaveText(0.95, 0.77, 0.99, 0.99, "brNDC");
  pt2->SetBorderSize(0);
	pt2->SetFillStyle(0);
	pt2->SetTextAlign(21);
	pt2->SetTextFont(42);
	pt2->SetTextSize(0.04);
	text = pt2->AddText(0.01,0.3,"#font[42]{-2 #Delta ln L}");
	text->SetTextSize(0.056);
	text->SetTextAngle(90);
	pt2->Draw();

  c->SaveAs(indir+"/"+canvasname+".pdf");
  c->SaveAs(indir+"/"+canvasname+".png");

  delete pt2;
  delete pt;
  delete l;
  c->Close();
  delete hh;
}
