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
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSpline.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TPaletteAxis.h"
#include "contours.h"


using namespace std;


class TPaletteAxisPatched : public TPaletteAxis{

  TPaletteAxisPatched() : TPaletteAxis(){}
  TPaletteAxisPatched(Double_t x1, Double_t y1, Double_t x2, Double_t  y2, TH1 *h) : TPaletteAxis(x1, y1, x2, y2, h){}
  TPaletteAxisPatched(const TPaletteAxisPatched& palette) : TPaletteAxis(palette){}
  TPaletteAxisPatched(const TPaletteAxis& palette) : TPaletteAxis(palette){}
  virtual ~TPaletteAxisPatched(){}
  virtual void Paint(Option_t* _opt_=""){
    ConvertNDCtoPad();

    SetFillStyle(1001);
    Double_t ymin = fY1;
    Double_t ymax = fY2;
    Double_t xmin = fX1;
    Double_t xmax = fX2;
    Double_t wmin = fH->GetMinimum();
    Double_t wmax = fH->GetMaximum();
    Double_t wlmin = wmin;
    Double_t wlmax = wmax;
    Double_t b1, b2, w1, w2, zc;
    Bool_t   kHorizontal = false;

    if ((wlmax - wlmin) <= 0) {
      Double_t mz = wlmin * 0.1;
      if (mz == 0) mz = 0.1;
      wlmin = wlmin - mz;
      wlmax = wlmax + mz;
      wmin  = wlmin;
      wmax  = wlmax;
    }

    if (GetX2NDC()-GetX1NDC() > GetY2NDC()-GetY1NDC()) kHorizontal = true;

    if (gPad->GetLogz()) {
      if (wmin <= 0 && wmax > 0) wmin = TMath::Min((Double_t) 1, (Double_t)0.001 * wmax);
      wlmin = TMath::Log10(wmin);
      wlmax = TMath::Log10(wmax);
    }
    Double_t ws    = wlmax - wlmin;
    Int_t ncolors = gStyle->GetNumberOfColors();
    Int_t ndivz = fH->GetContour();
    if (ndivz == 0) return;
    ndivz = TMath::Abs(ndivz);
    Int_t theColor, color;
    // import Attributes already here since we might need them for CJUST
    if (fH->GetDimension() == 2) fAxis.ImportAxisAttributes(fH->GetZaxis());
    // case option "CJUST": put labels directly at color boundaries
    TLatex *label = NULL;
    TLine *line = NULL;
    Double_t prevlab = 0;
    TString opt(fH->GetDrawOption());
    if (opt.Contains("CJUST", TString::kIgnoreCase)) {
      label = new TLatex();
      label->SetTextFont(fAxis.GetLabelFont());
      label->SetTextColor(fAxis.GetLabelColor());
      if (kHorizontal) label->SetTextAlign(kHAlignCenter+kVAlignTop);
      else             label->SetTextAlign(kHAlignLeft+kVAlignCenter);
      line = new TLine();
      line->SetLineColor(fAxis.GetLineColor());
      if (kHorizontal) line->PaintLine(xmin, ymin, xmax, ymin);
      else             line->PaintLine(xmax, ymin, xmax, ymax);
    }
    Double_t scale = ndivz / (wlmax - wlmin);
    for (Int_t i = 0; i < ndivz; i++) {

      zc = fH->GetContourLevel(i);
      if (fH->TestBit(TH1::kUserContour) && gPad->GetLogz())
        zc = TMath::Log10(zc);
      w1 = zc;
      if (w1 < wlmin) w1 = wlmin;

      w2 = wlmax;
      if (i < ndivz - 1) {
        zc = fH->GetContourLevel(i + 1);
        if (fH->TestBit(TH1::kUserContour) && gPad->GetLogz())
          zc = TMath::Log10(zc);
        w2 = zc;
      }

      if (w2 <= wlmin) continue;
      if (kHorizontal) {
        b1 = xmin + (w1 - wlmin) * (xmax - xmin) / ws;
        b2 = xmin + (w2 - wlmin) * (xmax - xmin) / ws;
      }
      else {
        b1 = ymin + (w1 - wlmin) * (ymax - ymin) / ws;
        b2 = ymin + (w2 - wlmin) * (ymax - ymin) / ws;
      }

      if (fH->TestBit(TH1::kUserContour)) {
        color = i;
      }
      else {
        color = Int_t(0.01 + (w1 - wlmin) * scale);
      }

      theColor = Int_t((color + 0.99) * Double_t(ncolors) / Double_t(ndivz));
      SetFillColor(gStyle->GetColorPalette(theColor));
      TAttFill::Modify();
      if (kHorizontal) gPad->PaintBox(b1, ymin, b2, ymax);
      else             gPad->PaintBox(xmin, b1, xmax, b2);
      // case option "CJUST": put labels directly
      if (label) {
        Double_t lof = fAxis.GetLabelOffset()*(gPad->GetUxmax()-gPad->GetUxmin());
        // the following assumes option "S"
        Double_t tlength = fAxis.GetTickSize() * (gPad->GetUxmax()-gPad->GetUxmin());
        Double_t lsize = fAxis.GetLabelSize();
        Double_t lsize_user = lsize*(gPad->GetUymax()-gPad->GetUymin());
        Double_t zlab = fH->GetContourLevel(i);
        if (gPad->GetLogz()&& !fH->TestBit(TH1::kUserContour)) {
          zlab = TMath::Power(10, zlab);
        }
        // make sure labels dont overlap
        if (i == 0 || (b1 - prevlab) > 1.5*lsize_user) {
          if (kHorizontal) label->PaintLatex(b1, ymin - lof, 0, lsize, Form("%g", zlab));
          else             label->PaintLatex(xmax + lof, b1, 0, lsize, Form("%g", zlab));
          prevlab = b1;
        }
        if (kHorizontal) line->PaintLine(b2, ymin+tlength, b2, ymin);
        else             line->PaintLine(xmax-tlength, b1, xmax, b1);
        if (i == ndivz-1) {
          // label + tick at top of axis
          if ((b2 - prevlab > 1.5*lsize_user)) {
            if (kHorizontal) label->PaintLatex(b2, ymin - lof, 0, lsize, Form("%g", fH->GetMaximum()));
            else             label->PaintLatex(xmax + lof, b2, 0, lsize, Form("%g", fH->GetMaximum()));
          }
          if (kHorizontal) line->PaintLine(b1, ymin+tlength, b1, ymin);
          else             line->PaintLine(xmax-tlength, b2, xmax, b2);
        }
      }
    }

    // Take primary divisions only
    Int_t ndiv = fH->GetZaxis()->GetNdivisions();
    Bool_t isOptimized = ndiv>0;
    Int_t absDiv = abs(ndiv);
    Int_t maxD = absDiv/1000000;
    ndiv = absDiv%100 + maxD*1000000;
    if (!isOptimized) ndiv  = -ndiv;

    char chopt[6] = "S   ";
    chopt[1] = 0;
    strncat(chopt, "+L", 3);
    if (ndiv < 0) {
      ndiv = TMath::Abs(ndiv);
      strncat(chopt, "N", 2);
    }
    if (gPad->GetLogz()) {
      wmin = TMath::Power(10., wlmin);
      wmax = TMath::Power(10., wlmax);
      strncat(chopt, "G", 2);
    }
    if (label) {
      // case option "CJUST", cleanup
      delete label;
      delete line;
    }
    else {
      // default
      if (kHorizontal) fAxis.PaintAxis(xmin, ymin, xmax, ymin, wmin, wmax, ndiv, chopt);
      else             fAxis.PaintAxis(xmax, ymin, xmax, ymax, wmin, wmax, ndiv, chopt);
    }
  }

  ClassDef(TPaletteAxisPatched, 1)
};

ClassImp(TPaletteAxisPatched);



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
  cout << "Running histogram fix..." << endl;
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
    cout << "Looping over the tree to find minimum -2dNLL..." << endl;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      deltaNLL = deltaNLL*2.;
      if (deltaNLL<0.) continue;
      minNLL = std::min(minNLL, deltaNLL);
    }
    std::cout << "Minimum deltaNLL = " << minNLL << std::endl;

    cout << "Looping over the tree to find the list of points..." << endl;
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
        doAdd &= !((strxvar=="rf_offshell" && xval>3.) || (stryvar=="rf_offshell" && yval>3.));
        if (doAdd){
          addByLowest(xvalList, xval, true);
          addByLowest(yvalList, yval, true);
        }
      }
    }
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
    cout << "Looping over the tree to fill the histogram..." << endl;
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

void plotScan2D(TString const indir, TString const strxvar, TString const stryvar, TString const strachypo=""){
  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.20);
  gROOT->ForceStyle();

  TChain* tree = new TChain("limit");
  {
    int nfiles = 0;
    std::vector<TString> strappendlist{ "", "_UCSD", "_uscms", "_centeredScan", "_centered", "_more", "_morelow", "_uscms_more", "_centeredScan_more", "_centered_more" };
    for (auto const& strappend:strappendlist) nfiles += tree->Add(indir + strappend + "/*.root");
    cout << "Added " << nfiles << " files for " << indir << "..." << endl;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(strxvar, 1);
    tree->SetBranchStatus(stryvar, 1);
    tree->SetBranchStatus("deltaNLL", 1);
  }

  TFile* ff = TFile::Open("__tmp__.root", "read");
  if (!ff) ff = TFile::Open("__tmp__.root", "recreate");
  bool doReuse = !ff->IsWritable();
  ff->cd();

  TH2F* hh = nullptr;
  TGraph* g68 = nullptr;
  TGraph* g95 = nullptr;
  TList* tlFF68 = nullptr;
  TList* tlFF95 = nullptr;
  TGraph* best = nullptr;
  constexpr double thr_dNLL_68 = 2.30; // 1 sigma
  constexpr double thr_dNLL_95 = 5.99; // 95%
  constexpr double thr_dNLL_3s = 11.83; // 3 sigma
  constexpr double thr_dNLL_4s = 19.33; // 4 sigma
  constexpr double thr_dNLL_5s = 28.74; // 5 sigma

  if (!doReuse){
    hh = getHistogramFromTree(tree, strxvar, stryvar);
    hh->SetName("hist");
    if (!hh) cout << "ERROR: No histogram!" << endl;

    hh->GetXaxis()->SetTitle(getVariableLabel(strxvar, strachypo));
    hh->GetYaxis()->SetTitle(getVariableLabel(stryvar, strachypo));

    tlFF68 = contourFromTH2(hh, thr_dNLL_68);
    tlFF95 = contourFromTH2(hh, thr_dNLL_95);
    if (tlFF68 && tlFF68->GetSize()>0) g68 = (TGraph*) tlFF68->At(0);
    if (tlFF95 && tlFF95->GetSize()>0) g95 = (TGraph*) tlFF95->At(0);
    if (g68){
      cout << "68% CL contour is present" << endl;
      g68->SetName("g68");
    }
    if (g95){
      cout << "95% CL contour is present" << endl;
      g95->SetName("g95");
    }

    best = bestFit(tree, strxvar, stryvar);
    best->SetName("best");
    if (!best) cout << "ERROR: Best fit not found!" << endl;
    else{
      if (strxvar=="GGsm") best->GetX()[0] *= 4.07;
      if (stryvar=="GGsm") best->GetY()[0] *= 4.07;
    }

    ff->WriteTObject(hh);
    ff->WriteTObject(best);
    ff->WriteTObject(g68);
    ff->WriteTObject(g95);
  }
  else{
    hh = (TH2F*) ff->Get("hist");
    best = (TGraph*) ff->Get("best");
    g68 = (TGraph*) ff->Get("g68");
    g95 = (TGraph*) ff->Get("g95");
  }

  if (best){
    double bestx, besty;
    best->GetPoint(0, bestx, besty);

    cout << "Best fit: " << bestx << ", " << besty << endl;
    best->SetMarkerStyle(5);
    best->SetMarkerSize(3);
  }
  if (g68){
    g68->SetLineColor(kBlack);
    g68->SetLineWidth(2);
    g68->SetLineStyle(6);
    cout << "68% CL contour is processed" << endl;
  }
  if (g95){
    g95->SetLineColor(kBlack);
    g95->SetLineWidth(2);
    g95->SetLineStyle(1);
    cout << "95% CL contour is processed" << endl;
  }

  double zmax=-9e9, zmin=0e0;
  for (int ix=1; ix<=hh->GetNbinsX(); ix++){
    for (int iy=1; iy<=hh->GetNbinsY(); iy++){
      double bc = hh->GetBinContent(ix, iy);
      zmax = std::max(zmax, bc);
      zmin = std::min(zmin, bc);
    }
  }
  cout << "Min -2dNLL = " << zmin << endl;
  cout << "Max -2dNLL = " << zmax << endl;
  hh->SetMinimum(zmin);
  hh->SetMaximum(zmax);

  constexpr unsigned int ncolors = 500;
  std::vector<int> colors; colors.reserve(ncolors);
  std::vector<double> Red_abs{ 0, 255, 255, 0, 255, 224 };
  std::vector<double> Green_abs{ 255, 0, 222, 230, 0, 224 };
  std::vector<double> Blue_abs{ 0, 255, 3, 230, 0, 224 };
  //std::vector<double> Red_abs{ 0, 255, 255, 255, 0, 224 };
  //std::vector<double> Green_abs{ 255, 0, 222, 0, 230, 224 };
  //std::vector<double> Blue_abs{ 0, 255, 3, 0, 230, 224 };
  std::vector<double> stops_dNLL{ 0.0, thr_dNLL_68, thr_dNLL_95, thr_dNLL_3s, thr_dNLL_4s, thr_dNLL_5s };
  std::vector<double> Red;
  std::vector<double> Green;
  std::vector<double> Blue;
  std::vector<double> Length;
  for (unsigned int ic=0; ic<stops_dNLL.size(); ic++){
    double srat = stops_dNLL.at(ic) / zmax;
    double rr = Red_abs.at(ic)/255.;
    double gg = Green_abs.at(ic)/255.;
    double bb = Blue_abs.at(ic)/255.;

    if (ic==0 || srat<=1.0){
      Red.push_back(rr);
      Green.push_back(gg);
      Blue.push_back(bb);
      Length.push_back(srat);
      if (srat==1.0) break;
    }
    else{
      double srat_prev = stops_dNLL.at(ic-1) / zmax;

      double rr_prev = Red_abs.at(ic-1)/255.;
      double gg_prev = Green_abs.at(ic-1)/255.;
      double bb_prev = Blue_abs.at(ic-1)/255.;

      rr = (rr-rr_prev)/(srat-srat_prev)*(1.-srat_prev)+rr_prev;
      gg = (gg-gg_prev)/(srat-srat_prev)*(1.-srat_prev)+gg_prev;
      bb = (bb-bb_prev)/(srat-srat_prev)*(1.-srat_prev)+bb_prev;

      Red.push_back(rr);
      Green.push_back(gg);
      Blue.push_back(bb);
      Length.push_back(1.0);

      break;
    }
  }
  for (unsigned int ic=0; ic<Length.size(); ic++){
    cout << "(R, G, B, L) = (" << Red.at(ic) << ", " << Green.at(ic) << ", " << Blue.at(ic) << ", " << Length.at(ic) << ")" << endl;
  }
  int FI = TColor::CreateGradientColorTable(Length.size(), Length.data(), Red.data(), Green.data(), Blue.data(), ncolors);
  if (FI<0) cout << "Failed to set the color palette." << endl;
  else{
    for (unsigned int ic=0; ic<ncolors; ic++) colors.push_back(FI+ic);
    gStyle->SetPalette(ncolors, colors.data());
    gStyle->SetNumberContours(ncolors);
  }

  TString canvasname = Form("cCompare_%sVS%s_%s", strxvar.Data(), stryvar.Data(), strachypo.Data());
  TCanvas* c = new TCanvas(canvasname, "", 2000, 1600);
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

  hh->GetXaxis()->SetLabelSize(0.04);
  //hh->GetXaxis()->SetRangeUser(0., 20.);
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetLabelSize(0.04);
  hh->GetYaxis()->SetTitleOffset(1);
  hh->GetYaxis()->CenterTitle();

   
  TLegend *l = new TLegend(0.475,0.72,0.835,0.92);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(42);
  l->SetTextSize(0.04);
  if (g95) l->AddEntry(g95,"-2 #Delta ln L = 5.99","l");
  if (g68) l->AddEntry(g68,"-2 #Delta ln L = 2.30","l");
  l->AddEntry(best,"Best fit","p");
  c->cd();
  hh->Draw("colz");
  if (g95) g95->Draw("csame");
  if (g68) g68->Draw("csame");

  TMarker m;
  m.SetMarkerSize(4);
  m.SetMarkerColor(kBlue+1);
  m.SetMarkerStyle(33);
  //  m.SetMarkerSize(1.5);
  //  m.SetMarkerColor(kBlack);
  //  m.SetMarkerStyle(5);
  m.DrawMarker(1, 1);
  l->AddEntry(&m, "SM", "p");
  best->Draw("Psame");
  l->Draw();

  TText* text=nullptr;
  TPaveText* ptc = new TPaveText(0.22, 0.85, 0.37, 0.95, "brNDC");
  ptc->SetBorderSize(0);
  ptc->SetFillStyle(0);
  ptc->SetTextAlign(12);
  ptc->SetTextFont(42);
  ptc->SetTextSize(0.04);
  text = ptc->AddText(0.02, 0.5, "#font[61]{CMS}");
  //text = ptc->AddText(0.12, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.06);
  ptc->Draw();

	TPaveText* pt = new TPaveText(0.15,0.955,0.8,1,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextAlign(32);
	pt->SetTextFont(42);
	pt->SetTextSize(0.04);
  float lumi2011=5.1;
  float lumi2012=19.7;
  float lumi2015=2.7;
  float lumi201617=77.5;
  float lumi20161718=138;
  //text = pt->AddText(0.68, 0.45, Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi201617+0.05));
  //text = pt->AddText(0.2, 0.45, Form("#font[42]{%.1f fb^{-1} (7 TeV) + %.1f fb^{-1} (8 TeV) + %.1f fb^{-1} (13 TeV)}", lumi2011, lumi2012, lumi2015+lumi201617+0.05));
  text = pt->AddText(0.999, 0.45, Form("#font[42]{#leq%.0f fb^{-1} (13 TeV)}", lumi20161718));
  text->SetTextSize(0.04);
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

  c->Update();
  c->Modified();

  c->SaveAs(canvasname+".pdf");
  c->SaveAs(canvasname+".png");
  c->SaveAs(canvasname+".root");
  c->SaveAs(canvasname+".C");

  delete pt2;
  delete pt;
  delete ptc;
  delete l;
  c->Close();

  delete tlFF95;
  delete tlFF68;
  delete hh;

  ff->Close();
  delete tree;
}
