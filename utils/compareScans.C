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


using namespace std;


struct GraphStyle{
  TString label;
  int marker_style;
  float marker_size;
  EColor marker_color;
  int line_style;
  float line_width;
  EColor line_color;
  GraphStyle(
    TString label_,
    int marker_style_,
    float marker_size_,
    EColor marker_color_,
    int line_style_,
    float line_width_,
    EColor line_color_
  ) :
    label(label_),
    marker_style(marker_style_),
    marker_size(marker_size_),
    marker_color(marker_color_),
    line_style(line_style_),
    line_width(line_width_),
    line_color(line_color_)
  {}
  GraphStyle(
    TString label_,
    int line_style_,
    float line_width_,
    EColor line_color_
  ) :
    label(label_),
    marker_style(1),
    marker_size(1),
    marker_color(line_color_),
    line_style(line_style_),
    line_width(line_width_),
    line_color(line_color_)
  {}
}

template<typename T> TGraph* makeGraphFromPair(std::vector<std::pair<T, T>> points, TString name){
  if (points.empty()) return nullptr;
  const unsigned int nbins = points.size();
  double xy[2][nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xy[0][bin] = points[bin].first;
    xy[1][bin] = points[bin].second;
  }
  TGraph* tg = new TGraph(nbins, xy[0], xy[1]);
  tg->SetName(name);
  return tg;
}

template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T const& val, U const& index, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first>=val){
        inserted=true;
        if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
        break;
      }
    }
  }
  if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
}

template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::pair<T, U> const& val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first==val.first){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first>=val.first){
        inserted=true;
        if ((*it).second!=val.second) valArray.insert(it, val);
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

TGraph* getGraphFromTree(TTree* tree, TString const strxvar, TString const stryvar){
  typedef float var_t;

  TGraph* gr=nullptr;
  vector<pair<var_t, var_t>> points;
  var_t xvar, yvar;

  if (tree){
    tree->SetBranchAddress(strxvar, &xvar);
    tree->SetBranchAddress(stryvar, &yvar);
    var_t minY=1e9;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      var_t xval = xvar;
      var_t yval = yvar;
      if (strxvar=="deltaNLL") xval *= 2;
      else if (strxvar=="GGsm") xval *= 4.07;
      if (stryvar=="deltaNLL") yval *= 2;
      else if (stryvar=="GGsm") yval *= 4.07;
      pair<var_t, var_t> point(xval, yval);
      minY=std::min(minY, yval);
      addByLowest(points, point, true);
    }
    if (stryvar=="deltaNLL"){ for (auto& point:points) point.second -= minY; }
    TString grname=Form("%s%s%s", strxvar.Data(), ":", stryvar.Data());
    gr=makeGraphFromPair(points, grname);
  }

  return gr;
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

void compareScans(vector<pair<TString, GraphStyle> indir_label_pair_list, TString const strxvar, TString const stryvar="deltaNLL", TString const strachypo=""){
  const float excl[2] ={ 1, 3.841 };

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);

  typedef float var_t;
  var_t xvar, yvar;
  vector<TGraph*> grlist;
  for (auto const& indir_label_pair:indir_label_pair_list){
    TString const& indir = indir_label_pair.first;
    TChain* tree=new TChain("limit");
    vector<TString> lsdirs = lsdir("./");
    for (auto const& strdir:lsdirs){
      if (strdir.Contains(indir)){
        vector<TString> lscontent = lsdir(strdir);
        for (auto const& s:lscontent){ if (s.Contains(".root")) tree->Add(strdir+"/"+s); }
      }
    }
    TGraph* gr = getGraphFromTree(tree, strxvar, stryvar);
    delete tree;
    grlist.push_back(gr);
  }

  gr->SetLineWidth(2);
  gr->SetMarkerColor(0);
  gr->SetLineColor(kBlack);
  gr->SetLineStyle(1);
  gr->SetLineWidth(2);

  TString canvasname = Form("cCompare_%sVS%s", stryvar.Data(), strxvar.Data());
  TCanvas* c1 = new TCanvas(canvasname, "", 800, 800);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->cd();
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);

  TLegend* leg = new TLegend(0.2, 0.7, 0.5, 0.92);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  float minY=0, maxY=-1;
  bool first=true;
  for (unsigned int ig=0; ig<indir_label_pair_list.size(); ig++){
    TGraph*& gr = grlist.at(ig);
    if (!gr) continue;
    float minX=gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->FindBin(gr->GetX()[0]));
    float maxX=gr->GetXaxis()->GetBinUpEdge(gr->GetXaxis()->FindBin(gr->GetX()[gr->GetN()-1]));
    gr->GetXaxis()->SetTitle(getVariableLabel(strxvar, strachypo));
    gr->GetYaxis()->SetTitle(getVariableLabel(stryvar, strachypo));
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetYaxis()->SetTitleFont(42);
    gr->GetYaxis()->SetNdivisions(505);
    gr->GetXaxis()->SetNdivisions(505);
    gr->GetYaxis()->CenterTitle();
    gr->GetXaxis()->CenterTitle();
    if (stryvar=="deltaNLL"){
      if (strxvar=="GGsm"){
        gr->GetXaxis()->SetNdivisions(510);
        minX=gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->FindBin(0.));
        maxX=gr->GetXaxis()->GetBinUpEdge(gr->GetXaxis()->FindBin(50.));
        maxY=-1; for (int ip=0; ip<gr->GetN(); ip++){ if (gr->GetX()[ip]<=maxX && gr->GetX()[ip]>=minX) maxY=std::max(maxY, (float) gr->GetY()[ip]); }
      }
      else if (strxvar=="CMS_zz4l_fai1"){
        gr->GetXaxis()->SetNdivisions(510);
        minX=gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->FindBin(-0.15));
        maxX=gr->GetXaxis()->GetBinUpEdge(gr->GetXaxis()->FindBin(0.15));
        maxY=-1; for (int ip=0; ip<gr->GetN(); ip++){ if (gr->GetX()[ip]<=maxX && gr->GetX()[ip]>=minX) maxY=std::max(maxY, (float) gr->GetY()[ip]); }
      }
      gr->GetXaxis()->SetRangeUser(minX, maxX);
      gr->GetYaxis()->SetRangeUser(minY, maxY*1.2);
    }
    else{
      gr->GetXaxis()->SetNdivisions(510);
      maxY=-9999; for (int ip=0; ip<gr->GetN(); ip++){ maxY=std::max(maxY, (float) gr->GetY()[ip]); }
      minY=9999; for (int ip=0; ip<gr->GetN(); ip++){ minY=std::min(minY, (float) gr->GetY()[ip]); }
      gr->GetYaxis()->SetRangeUser(minY*(minY>0. ? 0.8 : 1.2), maxY*1.2);
    }
    gr->SetMarkerStyle(indir_label_pair_list.at(ig).second.marker_style);
    gr->SetMarkerSize(indir_label_pair_list.at(ig).second.marker_size);
    gr->SetMarkerColor(indir_label_pair_list.at(ig).second.marker_color);
    gr->SetLineStyle(indir_label_pair_list.at(ig).second.line_style);
    gr->SetLineWidth(indir_label_pair_list.at(ig).second.line_width);
    gr->SetLineColor(indir_label_pair_list.at(ig).second.line_color);
    if (first) gr->Draw("ac");
    else gr->Draw("csame");
    leg->AddEntry(gr, indir_label_pair_list.at(ig).second.label, "l");
  }
  leg->Draw();

  TPaveText* pt = new TPaveText(0.15, 0.955, 0.96, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  TText* text = pt->AddText(0.02,0.45,"#font[61]{CMS}");
  //text = pt->AddText(0.14, 0.42, "#font[52]{Unpublished}");
  //text->SetTextSize(0.0315);
  //text = pt->AddText(0.48, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
  //text->SetTextSize(0.0315);
  pt->Draw();

  TPaveText* oneSig=nullptr; TLine* l1=nullptr;
  TPaveText* twoSig=nullptr; TLine* l2=nullptr;

  if (stryvar=="deltaNLL"){
    if (maxY*1.2>excl[0]){
      const float sigtext_minX=0.85;
      const float sigtext_widthX=0.05;
      const float sigtext_maxX=sigtext_minX+sigtext_widthX;
      const float sigtext_minY=excl[0]+0.1;
      const float sigtext_widthY=0.05*maxY*1.2/1.8;
      const float sigtext_maxY=sigtext_minY+sigtext_widthY;

      oneSig = new TPaveText(minX*(1.-sigtext_minX)+maxX*sigtext_minX, sigtext_minY, minX*(1.-sigtext_maxX)+maxX*sigtext_maxX, sigtext_maxY, "nb");
      oneSig->SetFillColor(0);
      oneSig->SetTextFont(42);
      oneSig->SetTextSize(0.0315);
      oneSig->SetTextColor(kBlack);
      oneSig->SetBorderSize(0);
      oneSig->AddText("68% CL");
      oneSig->Draw();

      l1=new TLine();
      l1->SetLineStyle(9);
      l1->SetLineWidth(2);
      l1->SetLineColor(kBlack);
      l1->DrawLine(minX, excl[0], maxX, excl[0]);
      l1->Draw("same");
    }

    if (maxY*1.2>excl[1]){
      const float sigtext_minX=0.85;
      const float sigtext_widthX=0.05;
      const float sigtext_maxX=sigtext_minX+sigtext_widthX;
      const float sigtext_minY=excl[1]+0.1;
      const float sigtext_widthY=0.05*maxY*1.2/1.8;
      const float sigtext_maxY=sigtext_minY+sigtext_widthY;

      twoSig = new TPaveText(minX*(1.-sigtext_minX)+maxX*sigtext_minX, sigtext_minY, minX*(1.-sigtext_maxX)+maxX*sigtext_maxX, sigtext_maxY, "nb");
      twoSig->SetFillColor(0);
      twoSig->SetTextFont(42);
      twoSig->SetTextSize(0.0315);
      twoSig->SetTextColor(kBlack);
      twoSig->SetBorderSize(0);
      twoSig->AddText("95% CL");
      twoSig->Draw();

      l2=new TLine();
      l2->SetLineStyle(9);
      l2->SetLineWidth(2);
      l2->SetLineColor(kBlack);
      l2->DrawLine(minX, excl[1], maxX, excl[1]);
      l2->Draw("same");
    }
  }

  c1->SaveAs(indir+"/"+canvasname+".pdf");
  c1->SaveAs(indir+"/"+canvasname+".png");

  delete l2; delete twoSig;
  delete l1; delete oneSig;
  delete pt;
  delete leg;
  c1->Close();
  delete gr;
}
