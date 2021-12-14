#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <sys/types.h>
#include <dirent.h>
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TGraph.h"
#include "TLine.h"
#include "TSpline.h"
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
};

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
  var_t nll, nll0;

  if (tree){
    tree->SetBranchAddress(strxvar, &xvar);
    tree->SetBranchAddress(stryvar, &yvar);
    if (stryvar=="deltaNLL"){
      tree->SetBranchAddress("nll", &nll);
      tree->SetBranchAddress("nll0", &nll0);
    }
    var_t minY=1e9;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      var_t xval = xvar;
      var_t yval = yvar;
      if (strxvar=="deltaNLL") xval *= 2;
      else if (strxvar=="GGsm") xval *= 4.07;
      if (stryvar=="deltaNLL"){
        yval += nll+nll0;
        yval *= 2;
      }
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
  else if (strvar=="R" || strvar=="r_offshell") return "#mu^{off-shell}";
  else if (strvar=="RF" || strvar=="rf_offshell") return "#mu_{F}^{off-shell}";
  else if (strvar=="RV" || strvar=="rv_offshell") return "#mu_{V}^{off-shell}";
  else if (strvar=="CMS_zz4l_fai1"){
    TString strai=strachypo.Data();
    if (strachypo=="L1") strai="Lambda1";
    if (strai.Contains("Lambda")) strai.Prepend("#");
    //return Form("f_{%s} cos(#phi_{#lower[-0.2]{%s}})", strai.Data(), strai.Data());
    return Form("f_{%s}", strai.Data());
  }
  else return strvar;
}

void compareScans(std::vector< std::pair<TString, GraphStyle> > const& indir_label_pair_list, TString strxvar, TString stryvar="deltaNLL", TString strachypo=""){
  // Magic numbers
  constexpr double npixels_stdframe_xy = 800;
  constexpr double relmargin_frame_left = 0.20;
  constexpr double relmargin_frame_right = 0.05;
  constexpr double relmargin_frame_CMS = 0.07;
  constexpr double relmargin_frame_XTitle = 0.15;
  constexpr double relsize_CMSlogo = 0.98;
  constexpr double relsize_CMSlogo_sqrts = 0.8;
  constexpr double relsize_XYTitle = 0.9;
  constexpr double relsize_XYLabel = 0.8;
  constexpr double offset_xlabel = 0.004;
  constexpr double offset_ylabel = 0.007;
  constexpr double offset_xtitle = 1.09;
  constexpr double offset_ytitle = 1.5;

  const double npixels_CMSlogo = npixels_stdframe_xy*relmargin_frame_CMS*relsize_CMSlogo;
  const double npixels_CMSlogo_sqrts = npixels_CMSlogo*relsize_CMSlogo_sqrts;
  const double npixels_XYTitle = npixels_CMSlogo*relsize_XYTitle;
  const double npixels_XYLabel = npixels_CMSlogo*relsize_XYLabel;

  const double npixels_x = int(
    npixels_stdframe_xy*(
      1.
      + relmargin_frame_left
      + relmargin_frame_right
      ) + 0.5
    );
  const double npixels_y = int(
    npixels_stdframe_xy*(
      relmargin_frame_CMS
      + 1.
      + relmargin_frame_XTitle
      ) + 0.5
    );

  const float excl[2] ={ 1, 3.841 };

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);

  typedef float var_t;
  var_t xvar, yvar;
  std::vector<TGraph*> grlist;
  std::unordered_map<TGraph*, TString> grlabels;
  {
    int igr = 0;
    for (auto const& indir_label_pair:indir_label_pair_list){
      TString const& indir = indir_label_pair.first;
      auto const& grstyle = indir_label_pair.second;

      TChain* tree = new TChain("limit");
      std::vector<TString> lsdirs = lsdir("./");
      for (auto const& strdir:lsdirs){
        if (strdir.Contains(indir)){
          vector<TString> lscontent = lsdir(strdir);
          for (auto const& s:lscontent){ if (s.Contains(".root")) tree->Add(strdir+"/"+s); }
        }
      }

      TGraph* gr = getGraphFromTree(tree, strxvar, stryvar);
      delete tree;

      if (!gr) continue;

      grlist.push_back(gr);

      gr->SetMarkerSize(grstyle.marker_size);
      gr->SetMarkerColor(grstyle.marker_color);
      gr->SetMarkerStyle(grstyle.marker_style);
      gr->SetLineColor(grstyle.line_color);
      gr->SetLineStyle(grstyle.line_style);
      gr->SetLineWidth(grstyle.line_width);
      gr->SetName(Form("gr%i", igr));
      gr->SetTitle("");
      grlabels[gr] = grstyle.label;
    }
  }

  TString canvasname = Form("cCompare_%sVS%s_%s", stryvar.Data(), strxvar.Data(), strachypo.Data());
  TCanvas* c1 = new TCanvas(canvasname, "", npixels_x, npixels_y);
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
  c1->SetLeftMargin(relmargin_frame_left/(1.+relmargin_frame_left+relmargin_frame_right));
  c1->SetRightMargin(relmargin_frame_right/(1.+relmargin_frame_left+relmargin_frame_right));
  c1->SetTopMargin(npixels_stdframe_xy*relmargin_frame_CMS/npixels_y);
  c1->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_XTitle/npixels_y);

  double leg_xmin = 0.2;
  double leg_xmax = 0.5;
  double leg_ymin = 0.7;
  double leg_ymax = 0.92;
  TLegend* leg = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(npixels_XYTitle);
  leg->SetTextAlign(12);

  float minX=0, maxX=0;
  float minY=0, maxY=-1;
  bool first=true;
  for (unsigned int ig=0; ig<indir_label_pair_list.size(); ig++){
    TGraph*& gr = grlist.at(ig);
    minX=gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->FindBin(gr->GetX()[0]));
    maxX=gr->GetXaxis()->GetBinUpEdge(gr->GetXaxis()->FindBin(gr->GetX()[gr->GetN()-1]));
    gr->GetXaxis()->SetTitle(getVariableLabel(strxvar, strachypo));
    gr->GetYaxis()->SetTitle(getVariableLabel(stryvar, strachypo));
    gr->GetXaxis()->SetLabelFont(43);
    gr->GetXaxis()->SetLabelOffset(offset_xlabel);
    gr->GetXaxis()->SetLabelSize(npixels_XYLabel);
    gr->GetXaxis()->SetTitleFont(43);
    gr->GetXaxis()->SetTitleSize(npixels_XYTitle);
    gr->GetXaxis()->SetTitleOffset(offset_xtitle);
    gr->GetYaxis()->SetLabelFont(43);
    gr->GetYaxis()->SetLabelOffset(offset_ylabel);
    gr->GetYaxis()->SetLabelSize(npixels_XYLabel);
    gr->GetYaxis()->SetTitleFont(43);
    gr->GetYaxis()->SetTitleSize(npixels_XYTitle);
    gr->GetYaxis()->SetTitleOffset(offset_ytitle);
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();
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

  constexpr bool markPreliminary = false;
  double lumi=138;
  TText* text;
  TPaveText pt(
    npixels_stdframe_xy*relmargin_frame_left/npixels_x,
    1.-(npixels_stdframe_xy*relmargin_frame_CMS-1)/npixels_y,
    1.-npixels_stdframe_xy*relmargin_frame_right/npixels_x,
    1,
    "brNDC"
  );
  pt.SetBorderSize(0);
  pt.SetFillStyle(0);
  pt.SetTextAlign(22);
  pt.SetTextFont(43);
  text = pt.AddText(0.001, 0.5, "CMS");
  text->SetTextFont(63);
  text->SetTextSize(npixels_CMSlogo);
  text->SetTextAlign(12);
  if (markPreliminary){
    text = pt.AddText(npixels_CMSlogo*2.2/npixels_stdframe_xy, 0.45, "Preliminary");
    text->SetTextFont(53);
    text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
    text->SetTextAlign(12);
  }
  int theSqrts=13;
  TString cErgTev = Form("#leq%.0f fb^{-1} (13 TeV)", lumi);
  text = pt.AddText(0.999, 0.45, cErgTev);
  text->SetTextFont(43);
  text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
  text->SetTextAlign(32);
  pt.Draw();

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

  c1->SaveAs(canvasname+".pdf");
  c1->SaveAs(canvasname+".png");

  delete l2; delete twoSig;
  delete l1; delete oneSig;
  delete leg;
  c1->Close();
  for (auto& gr:grlist) delete gr;
}
