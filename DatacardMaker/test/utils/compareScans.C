#include <iostream>
#include <fstream>
#include <sstream>
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
  int marker_color;
  int line_style;
  float line_width;
  int line_color;
  std::vector<double> knockout_x;
  GraphStyle(
    TString label_,
    int marker_style_,
    float marker_size_,
    int marker_color_,
    int line_style_,
    float line_width_,
    int line_color_
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
    int line_color_
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

TGraph* getGraphFromTree(TTree* tree, TString const strxvar, TString const stryvar, std::vector<double> const& knockout_x){
  typedef float var_t;

  TGraph* gr=nullptr;
  vector<pair<var_t, var_t>> points;
  var_t xvar, yvar;
  double nll, nll0;

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
        //yval += nll+nll0;
        yval *= 2;
        if (yval<0.) continue;
      }
      else if (stryvar=="GGsm") yval *= 4.07;

      bool doKnockout = false;
      for (double const& vk:knockout_x){
        if (std::abs((static_cast<double>(xval)-vk))<=((static_cast<double>(xval)+vk)/2.)*1e-4){
          doKnockout = true;
          break;
        }
      }
      if (doKnockout) continue;

      pair<var_t, var_t> point(xval, yval);
      minY=std::min(minY, yval);
      addByLowest(points, point, true);
    }
    if (stryvar=="deltaNLL"){ for (auto& point:points) point.second -= minY; }
    TString grname=Form("%s_%s", strxvar.Data(), stryvar.Data());
    gr=makeGraphFromPair(points, grname);
  }

  return gr;
}

TString getVariableLabel(TString const strvar, TString const strachypo){
  if (strvar=="deltaNLL") return "-2 #Delta lnL";
  else if (strvar=="GGsm") return "#Gamma_{H} (MeV)";
  else if (strvar=="rfv_offshell") return "#mu_{i}^{off-shell}";
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

std::vector<double> findIntersections(TGraph* gr, double const& excl){
  constexpr int niter = 100;
  std::vector<double> res;
  if (gr){
    TSpline3 sp("", gr, "b2e2", 0, 0);
    for (int ip=0; ip<gr->GetN()-1; ip++){
      if (
        (gr->GetY()[ip+1]>=excl && (gr->GetY()[ip]<excl || (ip==0 && gr->GetY()[ip]==excl)))
        ||
        (gr->GetY()[ip+1]<=excl && (gr->GetY()[ip]>excl || (ip==0 && gr->GetY()[ip]==excl)))
        ){
        double diffy=-1;
        double locx=-9e9;
        for (int j=0; j<=niter; j++){
          double xx=gr->GetX()[ip] + double(j)/double(niter)*(gr->GetX()[ip+1]-gr->GetX()[ip]);
          double yy=sp.Eval(xx);
          if (diffy<0. || diffy>std::abs(yy-excl)){
            locx = xx;
            diffy = std::abs(yy-excl);
          }
        }
        if (diffy>=0.) res.push_back(locx);
      }
    }
  }
  return res;
}

void print_hepdata_script(
  TString scriptname,
  TString const& stropt,
  std::vector< std::pair<TString, TGraph*> > const& label_gr_pairs,
  TString const& xtitle, TString const& ytitle,
  double const& xmin, double const& xmax
){
  ofstream tout(scriptname.Data());

  tout << R"V0G0N(
from hepdata_lib import Submission
from hepdata_lib import Table
from hepdata_lib import Variable

submission = Submission()
)V0G0N" << endl;

  for (unsigned int igr=0; igr<label_gr_pairs.size(); igr++){
    auto const& grlabel = label_gr_pairs.at(igr).first;
    auto const& gr = label_gr_pairs.at(igr).second;

    TString xlabel = gr->GetXaxis()->GetTitle();
    TString ylabel = gr->GetYaxis()->GetTitle();

    tout << Form("table = Table(\"%s vs %s (%s)\")", xtitle.Data(), ytitle.Data(), grlabel.Data()) << endl;
    tout << Form("table.description = \"%s vs %s (%s)\"", xtitle.Data(), ytitle.Data(), grlabel.Data()) << endl;
    tout << R"V0G0N(
table.location = "Table from CMS-HIG-21-013"
table.keywords["reactions"] = [ "HIGGS --> Z0 Z0" ]
)V0G0N" << endl;

    tout << Form("xvar = Variable(\"%s\", is_independent=True, is_binned=False, units=\"%s\")", xtitle.Data(), (xtitle.Contains("GGsm") ? "GeV" : "")) << endl;
    tout << Form("yvar = Variable(\"%s\", is_independent=False, is_binned=False, units=\"%s\")", ytitle.Data(), (ytitle.Contains("GGsm") ? "GeV" : "")) << endl;
    tout << Form("yvar.add_qualifier(\"%s\",\"\")", ylabel.Data()) << endl;
    tout << "yvar.add_qualifier(\"SQRT(S)\", 13, \"TeV\")" << endl;

    stringstream ssx, ssy;
    bool first = true;
    for (int ip=0; ip<gr->GetN(); ip++){
      if (gr->GetX()[ip]>=xmin && gr->GetX()[ip]<=xmax){
        if (!first){
          ssx << ", ";
          ssy << ", ";
        }
        ssx << gr->GetX()[ip];
        ssy << gr->GetY()[ip];
        first = false;
      }
    }
    tout << "xvar.values = [ " << ssx.str() << " ]" << endl;
    tout << "yvar.values = [ " << ssy.str() << " ]" << endl;

    tout << R"V0G0N(
table.add_variable(xvar)
table.add_variable(yvar)
table.keywords["cmenergies"] = [13000]
submission.add_table(table)
)V0G0N" << endl;
  }

  tout << Form("submission.create_files(\"hepdata_%s\",remove_old=True)", stropt.Data()) << endl;

  tout.close();
}

void compareScans_single(std::vector< std::pair<TString, GraphStyle> > const& indir_label_pair_list, TString strxvar, TString stryvar="deltaNLL", TString strachypo="", TString strappend="", TString strcase="", bool print_hepdata_rcd=false){
  // Magic numbers
  constexpr double npixels_stdframe_xy = 800;
  constexpr double relmargin_frame_left = 0.20;
  constexpr double relmargin_frame_right = 0.05;
  constexpr double relmargin_frame_CMS = 0.07;
  constexpr double relmargin_frame_XTitle = 0.15;
  constexpr double relsize_longy_extra = 0.3;
  constexpr double relsize_CMSlogo = 0.98;
  constexpr double relsize_CMSlogo_sqrts = 0.8;
  constexpr double relsize_XYTitle = 0.9;
  constexpr double relsize_XYLabel = 0.8;
  constexpr double offset_xlabel = 0.004;
  constexpr double offset_ylabel = 0.007;
  constexpr double offset_xtitle = 1.03;
  constexpr double offset_ytitle = 1.3;

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
      + 1. + (strcase!="GGsm_long" ? 0. : relsize_longy_extra)
      + relmargin_frame_XTitle
      ) + 0.5
    );

  const float excl[2] ={ 1, 3.841 };

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);

  TGraph* gr_obs = nullptr;
  TGraph* gr_exp = nullptr;

  typedef float var_t;
  var_t xvar, yvar;
  std::vector<TGraph*> grlist;
  std::unordered_map<TGraph*, TString> grindirs;
  std::unordered_map<TGraph*, GraphStyle const*> grstyles;
  unsigned int nleg = 0;
  {
    int igr = 0;
    for (auto& indir_label_pair:indir_label_pair_list){
      TString const& indir = indir_label_pair.first;
      auto& grstyle = indir_label_pair.second;

      TChain* tree = new TChain("limit");
      int ntrees = 0;
      ntrees += tree->Add(indir + "/*.root");
      ntrees += tree->Add(indir + "_centered/*.root");
      ntrees += tree->Add(indir + "_morelow/*.root");
      ntrees += tree->Add(indir + "_UCSD/*.root");
      cout << "Added " << ntrees << " from " << indir << " derivatives..." << endl;

      TString strxvar_eff = strxvar;
      if (strxvar=="rfv_offshell") strxvar_eff = (igr%2==0 ? "rf_offshell" : "rv_offshell");
      TGraph* gr = getGraphFromTree(tree, strxvar_eff, stryvar, grstyle.knockout_x);
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
      grstyles[gr] = &grstyle;
      grindirs[gr] = indir;

      if (grstyle.label!="") nleg++;

      cout << "Done with " << gr->GetName() << " with " << gr->GetN() << " points..." << endl;

      if (!gr_obs && grstyle.line_style==1){
        gr_obs = (TGraph*) gr->Clone("gr_obs");
        gr_obs->SetLineColor(kGray);
      }
      if (!gr_exp && grstyle.line_style!=1){
        gr_exp = (TGraph*) gr->Clone("gr_exp");
        gr_exp->SetLineColor(kGray);
      }

      igr++;
    }
  }

  if (strappend!="") strappend = strappend + "_";
  TString canvasname_core = Form("%s%sVS%s_%s", strappend.Data(), stryvar.Data(), strxvar.Data(), strachypo.Data());
  TString exclname = Form("limits_%s%s", canvasname_core.Data(), ".txt");

  ofstream tout(exclname.Data());
  for (auto const& pp:grindirs){
    auto const& gr = pp.first;
    auto const& indir = pp.second;

    cout << indir << endl;
    tout << indir << endl;
    {
      double llAtZero=-99;
      double xAtMin=-99;
      double llmin=9e9;
      for (int ix=0; ix<gr->GetN(); ix++){
        double const& yy = gr->GetY()[ix];
        double const& xx = gr->GetX()[ix];
        if (yy<llmin){
          llmin = yy;
          xAtMin = xx;
        }
        if (xx==0. && stryvar=="deltaNLL"){
          llAtZero = yy;
        }
      }
      
      cout << "Minimum: " << xAtMin << " @ -2dNLL = " << llmin << endl;
      tout << "Minimum: " << xAtMin << " @ -2dNLL = " << llmin << endl;
      if (llAtZero>=0.){
        cout << "-2dNLL = " << llAtZero << " @ " << strxvar << " = 0" << endl;
        tout << "-2dNLL = " << llAtZero << " @ " << strxvar << " = 0" << endl;
      }
    }
    for (unsigned char ie=0; ie<2; ie++){
      std::vector<double> intersections = findIntersections(gr, excl[ie]);

      cout << (ie==0 ? "68%:" : "95%:");
      tout << (ie==0 ? "68%:" : "95%:");
      for (auto const& px:intersections){
        cout << " " << px;
        tout << " " << px;
      }
      cout << endl;
      tout << endl;
    }
    cout << endl;
    tout << endl;
  }
  tout.close();

  TString canvasname = Form("cCompare_%s", canvasname_core.Data());
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

  double leg_xmin = (relmargin_frame_left+0.03)/(1.+relmargin_frame_left+relmargin_frame_right);
  if (strcase=="fa2"){
    leg_xmin = (relmargin_frame_left+0.45)/(1.+relmargin_frame_left+relmargin_frame_right);
  }
  else if (strcase=="fa3"){
    leg_xmin = (relmargin_frame_left+0.50)/(1.+relmargin_frame_left+relmargin_frame_right);
  }
  else if (strcase=="fL1"){
    leg_xmin = (relmargin_frame_left+0.4)/(1.+relmargin_frame_left+relmargin_frame_right);
  }
  double leg_xmax = leg_xmin + 0.52/(1.+relmargin_frame_left+relmargin_frame_right);
  double leg_ymax = static_cast<double>(npixels_y - npixels_stdframe_xy*(relmargin_frame_CMS+0.03))/static_cast<double>(npixels_y);
  double leg_ymin = leg_ymax - static_cast<double>(nleg)*1.5*npixels_XYTitle/npixels_y;
  TLegend* leg = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(npixels_XYTitle);
  leg->SetTextAlign(12);

  if (strcase=="fa3"){
    double leg_dx = leg_xmax - leg_xmin;
    leg_xmin = (relmargin_frame_left+0.12)/(1.+relmargin_frame_left+relmargin_frame_right);
    leg_xmax = leg_xmin + leg_dx;
    leg_ymin = leg_ymax - 2.*1.5*npixels_XYTitle/npixels_y;
  }
  else{
    leg_ymax = leg_ymin - npixels_XYTitle/npixels_y;
    leg_ymin = leg_ymax - 2.*1.5*npixels_XYTitle/npixels_y;
  }
  TLegend* leg2 = nullptr;
  if (gr_obs && gr_exp){
    leg2 = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
    leg2->SetFillColor(0);
    leg2->SetLineColor(0);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(43);
    leg2->SetTextSize(npixels_XYTitle);
    leg2->SetTextAlign(12);
    leg2->AddEntry(gr_obs, "Observed", "l");
    leg2->AddEntry(gr_exp, "Expected", "l");
  }


  int ndiv_x = 510;
  double minX=0, maxX=0;
  double minY=1e9, maxY=-1e9;
  bool first=true;
  for (TGraph*& gr:grlist){
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

    if (first){
      minX = gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->FindBin(gr->GetX()[0]));
      maxX = gr->GetXaxis()->GetBinUpEdge(gr->GetXaxis()->FindBin(gr->GetX()[gr->GetN()-1]));
    }
    else{
      minX = std::max(minX, (double) gr->GetXaxis()->GetBinLowEdge(gr->GetXaxis()->FindBin(gr->GetX()[0])));
      maxX = std::min(maxX, (double) gr->GetXaxis()->GetBinUpEdge(gr->GetXaxis()->FindBin(gr->GetX()[gr->GetN()-1])));
    }

    for (int ip=0; ip<gr->GetN(); ip++){
      if (gr->GetX()[ip]<minX) continue;
      if (gr->GetX()[ip]>maxX) continue;
      minY = std::min(minY, gr->GetY()[ip]);
      maxY = std::max(maxY, gr->GetY()[ip]);
      //cout << gr->GetName() << " (x, y) = " << gr->GetX()[ip] << ", " << gr->GetY()[ip] << endl;
    }
    first=false;
  }

  // Special x min/max
  constexpr double adj_maxy = 1.2;
  if (stryvar=="deltaNLL"){
    if (strxvar=="r_offshell"){
      minX = 0; maxX = 4; ndiv_x = 505;
      minY = 0; maxY = 15./adj_maxy;
    }
    else if (strxvar=="rfv_offshell" || strxvar=="rf_offshell" || strxvar=="rv_offshell"){
      minX = 0; maxX = 5; ndiv_x = 510;
      minY = 0; maxY = 15./adj_maxy;
    }
    else if (strxvar=="GGsm"){
      minX = 0; maxX = 15; ndiv_x = 505;
      minY = 0; maxY = 15./adj_maxy;
    }
    else if (strxvar=="CMS_zz4l_fai1"){
      if (strachypo=="a2"){
        minX = -0.012; maxX = 0.012; ndiv_x = 505;
        minY = 0; maxY = 15./adj_maxy;
      }
      else if (strachypo=="a3"){
        minX = -0.006; maxX = 0.006; ndiv_x = 503;
        minY = 0; maxY = 17./adj_maxy;
      }
      else if (strachypo=="L1"){
        minX = -0.0012; maxX = 0.0012; ndiv_x = 1003;
        minY = 0; maxY = 15./adj_maxy;
      }
    }
  }

  if (print_hepdata_rcd){
    gSystem->mkdir("hepdata", true);
    TString hepdataname = Form("hepdata/%s%s", canvasname_core.Data(), ".py");

    std::vector< std::pair<TString, TGraph*> > label_gr_pairs; label_gr_pairs.reserve(grlist.size());
    {
      unsigned int igr = 0;
      for (TGraph*& gr:grlist){
        auto const& grstyle = *(grstyles[gr]);

        TString tmplabel = grstyle.label;
        if (tmplabel=="") tmplabel = grstyles[grlist.at(igr - grlist.size()/2)]->label + " expected";
        else tmplabel += " observed";

        label_gr_pairs.emplace_back(tmplabel, gr);

        igr++;
      }
    }

    TString tmpxvar = strxvar;
    if (tmpxvar=="CMS_zz4l_fai1") tmpxvar = Form("f%s", strachypo.Data());

    print_hepdata_script(
      hepdataname,
      strcase,
      label_gr_pairs,
      tmpxvar, (stryvar=="deltaNLL" ? "-2dNLL" : stryvar),
      minX, maxX
    );
  }

  first = true;
  for (TGraph*& gr:grlist){
    auto const& grstyle = *(grstyles[gr]);

    gr->GetXaxis()->SetRangeUser(minX, maxX);
    gr->GetXaxis()->SetNdivisions(ndiv_x);
    gr->GetYaxis()->SetRangeUser(minY*(minY>0. ? 0.8 : 1.2), maxY*adj_maxy);

    if (first) gr->Draw("ac");
    else gr->Draw("csame");
    if (grstyle.label!="") leg->AddEntry(gr, grstyle.label, "l");
    first=false;
  }
  leg->Draw();
  if (leg2) leg2->Draw();

  constexpr bool markPreliminary = false;
  constexpr bool markSupplementary = false;
  double lumi=(strxvar=="GGsm" || strxvar=="CMS_zz4l_fai1" ? 140. : 138.);
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
  else if (markSupplementary){
    text = pt.AddText(npixels_CMSlogo*2.2/npixels_stdframe_xy, 0.45, "Supplementary");
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
    float sigtext_minX=0.77;
    if (strcase=="r_offshell" || strcase=="rv_offshell") sigtext_minX = 0.10;
    else if (strcase=="GGsm_long") sigtext_minX = 0.12;
    else if (strcase=="rfv_offshell_paper") sigtext_minX = 0.07;
    else if (strcase=="fL1") sigtext_minX = 0.05;
    const float sigtext_widthX=0.10;
    const float sigtext_maxX=sigtext_minX+sigtext_widthX;
    if (maxY*1.2>excl[0]){
      const float sigtext_minY=excl[0] + (maxY-minY)/npixels_stdframe_xy*npixels_XYLabel*0.1;
      const float sigtext_widthY=(maxY-minY)/npixels_stdframe_xy*npixels_XYLabel;
      const float sigtext_maxY=sigtext_minY+sigtext_widthY;

      oneSig = new TPaveText(minX*(1.-sigtext_minX)+maxX*sigtext_minX, sigtext_minY, minX*(1.-sigtext_maxX)+maxX*sigtext_maxX, sigtext_maxY, "nb");
      oneSig->SetFillColor(0);
      oneSig->SetTextFont(43);
      oneSig->SetTextAlign(12);
      oneSig->SetTextSize(npixels_XYLabel);
      oneSig->SetTextColor(kBlack);
      oneSig->SetBorderSize(0);
      oneSig->AddText("68% CL");
      oneSig->Draw();

      l1 = new TLine();
      l1->SetLineStyle(9);
      l1->SetLineWidth(2);
      l1->SetLineColor(kOrange-3);
      l1->DrawLine(minX, excl[0], maxX, excl[0]);
      l1->Draw("same");
    }

    if (maxY*1.2>excl[1]){
      const float sigtext_maxY=excl[1] - (maxY-minY)/npixels_stdframe_xy*npixels_XYLabel*0.1;
      const float sigtext_widthY=(maxY-minY)/npixels_stdframe_xy*npixels_XYLabel;
      const float sigtext_minY=sigtext_maxY-sigtext_widthY;

      twoSig = new TPaveText(minX*(1.-sigtext_minX)+maxX*sigtext_minX, sigtext_minY, minX*(1.-sigtext_maxX)+maxX*sigtext_maxX, sigtext_maxY, "nb");
      twoSig->SetFillColor(0);
      twoSig->SetTextFont(43);
      twoSig->SetTextAlign(12);
      twoSig->SetTextSize(npixels_XYLabel);
      twoSig->SetTextColor(kBlack);
      twoSig->SetBorderSize(0);
      twoSig->AddText("95% CL");
      twoSig->Draw();

      l2 = new TLine();
      l2->SetLineStyle(9);
      l2->SetLineWidth(2);
      l2->SetLineColor(kOrange-3);
      l2->DrawLine(minX, excl[1], maxX, excl[1]);
      l2->Draw("same");
    }
  }

  c1->SaveAs(canvasname+".pdf");
  c1->SaveAs(canvasname+".png");

  delete l2; delete twoSig;
  delete l1; delete oneSig;
  delete leg2;
  delete leg;
  c1->Close();
  for (auto& gr:grlist) delete gr;
  delete gr_exp;
  delete gr_obs;
}

void compareScans(TString strvar, bool print_hepdata_rcd=false){
  std::vector< std::pair<TString, GraphStyle> > inputs;

  TString stryvar = "deltaNLL";
  TString strappend = "";
  TString strachypo = "SM";
  TString strxvar = strvar;
  if (strvar=="r_offshell_paper"){
    strxvar = "r_offshell";
    strappend = "Paper";
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFixed_wCRZW_Obs_210812", {
          "R_{V,F}^{off-shell}=1",
          1, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFloated_wCRZW_Obs_210812", {
          "R_{V,F}^{off-shell} (u)",
          1, 2, static_cast<int>(kViolet)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFixed_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFloated_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kViolet)
        }
        )
    );
  }
  else if (strvar=="rfv_offshell_paper"){
    strxvar = "rfv_offshell";
    strappend = "Paper";
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RF_RVFloated_wCRZW_Obs_210812", {
          "#mu_{F}^{off-shell}",
          1, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RV_RFFloated_wCRZW_Obs_210812", {
          "#mu_{V}^{off-shell}",
          1, 2, static_cast<int>(kViolet)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RF_RVFloated_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RV_RFFloated_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kViolet)
        }
        )
    );
  }
  else if (strvar=="rf_offshell"){
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RF_RVFloated_wCRZW_Obs_210812", {
          "2l2#nu+4l",
          1, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_RF_RVFloated_wCRZW_Obs_210809", {
          "Only 2l2#nu",
          1, 2, static_cast<int>(kGreen+2)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RF_RVFloated_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_RF_RVFloated_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kGreen+2)
        }
        )
    );
  }
  else if (strvar=="rv_offshell"){
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RV_RFFloated_wCRZW_Obs_210812", {
          "2l2#nu+4l",
          1, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_RV_RFFloated_wCRZW_Obs_210809", {
          "Only 2l2#nu",
          1, 2, static_cast<int>(kGreen+2)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_RV_RFFloated_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_RV_RFFloated_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kGreen+2)
        }
        )
    );
  }
  else if (strvar=="r_offshell"){
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFixed_wCRZW_Obs_210812", {
          "R_{V,F}^{off-shell}=1 (2l2#nu+4l)",
          1, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFloated_wCRZW_Obs_210812", {
          "R_{V,F}^{off-shell} (u, 2l2#nu+4l)",
          1, 2, static_cast<int>(kViolet)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_R_RVFixed_wCRZW_Obs_210809", {
          "R_{V,F}^{off-shell}=1 (2l2#nu)",
          1, 2, static_cast<int>(kRed)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_R_RVFloated_wCRZW_Obs_210809", {
          "R_{V,F}^{off-shell} (u, 2l2#nu)",
          1, 2, static_cast<int>(kGreen+2)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFixed_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/OffshellOnly_R_RVFloated_wCRZW_Exp_210812", {
          "",
          7, 2, static_cast<int>(kViolet)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_R_RVFixed_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kRed)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/Offshell2L2NuOnly_R_RVFloated_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kGreen+2)
        }
        )
    );
  }
  else if (strvar=="GGsm"){
    strachypo = "all";
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_wCRZW_Obs_210809", {
          "SM-like (f_{ai}=0)",
          1, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_fa2Floated_wCRZW_Obs_210809", {
          "f_{a2} (u)",
          1, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_fa3Floated_wCRZW_Obs_210809", {
          "f_{a3} (u)",
          1, 2, static_cast<int>(kRed)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_fL1Floated_wCRZW_Obs_210809", {
          "f_{#Lambda1} (u)",
          1, 2, static_cast<int>(kViolet)
        }
        )
    );

    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_fa2Floated_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_fa3Floated_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kRed)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_fL1Floated_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kViolet)
        }
        )
    );
  }
  else if (strvar=="GGsm_long"){
    strxvar = "GGsm";
    strachypo = "SM";
    strappend = "2l2nu_vs_4l";
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_wCRZW_Obs_210809", {
          "2l2#nu+4l off-shell + 4l on-shell",
          1, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/PartialComb_NoOffshell4L_GGsm_Obs_220108", {
          "2l2#nu off-shell + 4l on-shell",
          1, 2, static_cast<int>(kGreen+2)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/PartialComb_NoOffshell2L2Nu_GGsm_Obs_220108", {
          "4l off-shell + 4l on-shell",
          1, 2, static_cast<int>(kViolet)
        }
        )
    );

    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/FullComb_GGsm_wCRZW_Exp_210809", {
          "",
          7, 2, static_cast<int>(kBlack)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/PartialComb_NoOffshell4L_GGsm_Exp_220108", {
          "",
          7, 2, static_cast<int>(kGreen+2)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        "AllResults/PartialComb_NoOffshell2L2Nu_GGsm_Exp_220108", {
          "",
          7, 2, static_cast<int>(kViolet)
        }
        )
    );
  }
  else if (strvar.BeginsWith("f")){
    strxvar = "CMS_zz4l_fai1";
    strachypo = strvar;
    strachypo = strachypo(1, strachypo.Length());
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        Form("AllResults/FullComb_f%s_GGsmFixed_wCRZW_Obs_210809", strachypo.Data()), {
          "#Gamma_{H}=4.1 MeV",
          1, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        Form("AllResults/FullComb_f%s_GGsmFloated_wCRZW_Obs_210809", strachypo.Data()), {
          "#Gamma_{H} (u)",
          1, 2, static_cast<int>(kViolet)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        Form("AllResults/FullComb_f%s_Onshell4LOnly_Obs_210912", strachypo.Data()), {
          "On-shell 4l",
          1, 2, static_cast<int>(kGreen+2)
        }
        )
    );

    inputs.push_back(
      std::pair<TString, GraphStyle>(
        Form("AllResults/FullComb_f%s_GGsmFixed_wCRZW_Exp_210809", strachypo.Data()), {
          "",
          7, 2, static_cast<int>(kBlue)
        }
        )
    );
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        Form("AllResults/FullComb_f%s_GGsmFloated_wCRZW_Exp_210809", strachypo.Data()), {
          "",
          7, 2, static_cast<int>(kViolet)
        }
        )
    );
    if (strvar=="fL1") inputs.back().second.knockout_x.push_back(0.00114533);
    inputs.push_back(
      std::pair<TString, GraphStyle>(
        Form("AllResults/FullComb_f%s_Onshell4LOnly_Exp_210912", strachypo.Data()), {
          "",
          7, 2, static_cast<int>(kGreen+2)
        }
        )
    );
  }


  compareScans_single(inputs, strxvar, stryvar, strachypo, strappend, strvar, print_hepdata_rcd);
}
