#include <unistd.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <functional>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include "TIterator.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TSystem.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "RooCmdArg.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooRealIntegral.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "QuantFuncMathCore.h"
#include <HiggsAnalysis/CombinedLimit/interface/AsymPow.h>


namespace std{

  template<> struct hash<TString>{
    typedef TString argument_type;
    typedef size_t result_type;
    result_type operator()(argument_type const& arg) const{ return hash<string>{}(arg.Data()); }
  };

}

using namespace RooFit;
using namespace std;

struct process_spec{
  RooAbsPdf* pdf;
  RooAbsReal* norm;
  double rate;
  TString name;
  std::vector<RooAbsReal*> rateModifiers;

  process_spec() : pdf(0), norm(0), rate(0){}
  process_spec(TString name_, double rate_) : pdf(nullptr), norm(nullptr), rate(rate_), name(name_){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()){}
  process_spec(const process_spec& other) : pdf(other.pdf), norm(other.norm), rate(other.rate), name(other.name){}

  void setRateModifiers(std::vector<RooAbsReal*> const& vlist){ rateModifiers=vlist; }
  double getExtraNorm();
};
double process_spec::getExtraNorm(){
  double res = 1;
  for (auto const& var:rateModifiers) res *= var->getVal();
  return res;
}


template <typename T> void divideBinWidth(T* histo);
template<> void divideBinWidth<TH1F>(TH1F* histo);
template<> void divideBinWidth<TH2F>(TH2F* histo);
template<> void divideBinWidth<TH3F>(TH3F* histo);
template <typename T> void multiplyBinWidth(T* histo);
template<> void multiplyBinWidth<TH1F>(TH1F* histo);
template<> void multiplyBinWidth<TH2F>(TH2F* histo);
template<> void multiplyBinWidth<TH3F>(TH3F* histo);
template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, bool useWidth, double* error=nullptr);
template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error=nullptr);
template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error=nullptr);
TH1F* getHistogramSlice(TH1F const* histo, TString newname);
TH1F* getHistogramSlice(TH2F const* histo, unsigned char XDirection, int iy, int jy, TString newname="");
TH1F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, int iy, int jy, int iz, int jz, TString newname=""); // "y" and "z" are cylical, so if Xdirection==1 (Y), "y"=Z and "z"=X
TH2F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, unsigned char YDirection, int iz, int jz, TString newname="");


void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
void splitOption(const TString rawoption, TString& wish, TString& value, char delimiter);
void splitOptionRecursive(const TString rawoption, vector<TString>& splitoptions, char delimiter);


unsigned int extractTemplates(process_spec& proc, RooDataSet* data, unordered_map<TString, TH1F>& procshape_1D, unordered_map<TString, TH2F>& procshape_2D, unordered_map<TString, TH3F>& procshape_3D, TString newname=""){
  vector<RooRealVar*> deps;
  RooArgSet* depList = proc.pdf->getDependents(data);
  TIterator* coefIter = depList->createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*) coefIter->Next())) deps.push_back((RooRealVar*) coef);
  delete coefIter;
  delete depList;
  const unsigned int ndims=deps.size();

  vector<RooRealVar*> pars;
  RooRealVar* fai1=0;
  RooRealVar* GGsm=0;
  RooArgSet* parList = proc.pdf->getParameters(data);
  coefIter = parList->createIterator();
  coef=0;
  while ((coef = (RooAbsArg*) coefIter->Next())){
    pars.push_back((RooRealVar*) coef);
    if (TString(pars.at(pars.size()-1)->GetName()).Contains("fai1")) fai1 = pars.at(pars.size()-1);
    else if (TString(pars.at(pars.size()-1)->GetName()).Contains("GGsm")) GGsm = pars.at(pars.size()-1);
  }
  delete coefIter;
  delete parList;

  RooBinning xbinning;
  RooBinning ybinning;
  RooBinning zbinning;
  RooCmdArg xcmd;
  if (deps.at(0)->getBins()>100){
    xbinning=RooBinning(deps.at(0)->getBins()/100, deps.at(0)->getMin(), deps.at(0)->getMax(), "plotbinningX");
    xcmd = Binning(xbinning);
  }
  RooCmdArg ycmd;
  RooCmdArg zcmd;
  if (deps.size()>1){
    if (deps.at(1)->getBins()>100){
      ybinning=RooBinning(deps.at(1)->getBins()/100, deps.at(1)->getMin(), deps.at(1)->getMax(), "plotbinningY");
      ycmd = YVar(*(deps.at(1)), Binning(ybinning));
    }
    else ycmd = YVar(*(deps.at(1)));
  }
  if (deps.size()>2){
    if (deps.at(2)->getBins()>100){
      zbinning=RooBinning(deps.at(2)->getBins()/100, deps.at(2)->getMin(), deps.at(2)->getMax(), "plotbinningZ");
      zcmd = ZVar(*(deps.at(2)), Binning(zbinning));
    }
    else zcmd = ZVar(*(deps.at(2)));
  }

  TString procname=proc.name;
  if (newname!="") procname=newname;
  TString tplname = "T_";
  tplname += procname;

  if (ndims==3){
    TH3F* tpl=(TH3F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    multiplyBinWidth(tpl);
    double normval = proc.rate * proc.getExtraNorm(); if (proc.norm) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
    tpl->SetName(tplname);
    tpl->SetTitle("");
    tpl->Scale(scale);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), 1, tpl->GetNbinsZ(), false, nullptr) << "?=" << normval << endl;

    bool isAdded=false;
    for (auto it=procshape_3D.begin(); it!=procshape_3D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      procshape_3D[procname]=TH3F(*tpl);
      procshape_3D[procname].SetName(tplname);
    }
    cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_3D[procname]), 1, procshape_3D[procname].GetNbinsX(), 1, procshape_3D[procname].GetNbinsY(), 1, procshape_3D[procname].GetNbinsZ(), false, nullptr) << "?=" << normval << endl;
    delete tpl;
  }
  else if (ndims==2){
    TH2F* tpl=(TH2F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    multiplyBinWidth(tpl);
    double normval = proc.rate * proc.getExtraNorm(); if (proc.norm) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
    tpl->SetName(tplname);
    tpl->SetTitle("");
    tpl->Scale(scale);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), false, nullptr) << "?=" << normval << endl;

    bool isAdded=false;
    for (auto it=procshape_2D.begin(); it!=procshape_2D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      procshape_2D[procname]=TH2F(*tpl);
      procshape_2D[procname].SetName(tplname);
    }
    cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_2D[procname]), 1, procshape_2D[procname].GetNbinsX(), 1, procshape_2D[procname].GetNbinsY(), false, nullptr) << "?=" << normval << endl;
    delete tpl;
  }
  else if (ndims==1){
    TH1F* tpl=(TH1F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    multiplyBinWidth(tpl);
    double normval = proc.rate * proc.getExtraNorm(); if (proc.norm) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
    tpl->SetName(tplname);
    tpl->SetTitle("");
    tpl->Scale(scale);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), false, nullptr) << "?=" << normval << endl;

    bool isAdded=false;
    for (auto it=procshape_1D.begin(); it!=procshape_1D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      procshape_1D[procname]=TH1F(*tpl);
      procshape_1D[procname].SetName(tplname);
    }
    cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_1D[procname]), 1, procshape_1D[procname].GetNbinsX(), false, nullptr) << "?=" << normval << endl;
    delete tpl;
  }
  cout << "Work on template " << tplname << " is complete!" << endl;
  return ndims;
}
unsigned int extractDataTemplates(process_spec& proc, RooDataSet* data, unordered_map<TString, TH1F>& procshape_1D, unordered_map<TString, TH2F>& procshape_2D, unordered_map<TString, TH3F>& procshape_3D, TString newname=""){
  vector<RooRealVar*> deps;
  RooArgSet* depList = proc.pdf->getDependents(data);
  TIterator* coefIter = depList->createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*) coefIter->Next())) deps.push_back((RooRealVar*) coef);
  delete coefIter;
  delete depList;
  const unsigned int ndims=deps.size();

  RooCmdArg xcmd;
  if (deps.at(0)->getBins()>100){
    RooBinning xbinning=RooBinning(deps.at(0)->getBins()/100, deps.at(0)->getMin(), deps.at(0)->getMax(), "plotbinningX");
    xcmd = Binning(xbinning);
  }
  RooCmdArg ycmd;
  RooCmdArg zcmd;
  if (deps.size()>1){
    if (deps.at(1)->getBins()>100){
      RooBinning ybinning=RooBinning(deps.at(1)->getBins()/100, deps.at(1)->getMin(), deps.at(1)->getMax(), "plotbinningY");
      ycmd = YVar(*(deps.at(1)), Binning(ybinning));
    }
    else ycmd = YVar(*(deps.at(1)));
  }
  if (deps.size()>2){
    if (deps.at(2)->getBins()>100){
      RooBinning zbinning=RooBinning(deps.at(2)->getBins()/100, deps.at(2)->getMin(), deps.at(2)->getMax(), "plotbinningZ");
      zcmd = ZVar(*(deps.at(2)), Binning(zbinning));
    }
    else zcmd = ZVar(*(deps.at(2)));
  }

  TString procname="data";
  if (newname!="") procname=newname;
  TString tplname = "T_";
  tplname += procname;

  if (ndims==3){
    TH3F* tpl=(TH3F*) data->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    tpl->SetName(tplname);
    tpl->SetTitle("");
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), 1, tpl->GetNbinsZ(), false, nullptr) << endl;

    bool isAdded=false;
    for (auto it=procshape_3D.begin(); it!=procshape_3D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      cout << "\t- Creating " << tplname << endl;
      procshape_3D[procname]=TH3F(*tpl);
      procshape_3D[procname].SetName(tplname);
    }
    cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_3D[procname]), 1, procshape_3D[procname].GetNbinsX(), 1, procshape_3D[procname].GetNbinsY(), 1, procshape_3D[procname].GetNbinsZ(), false, nullptr) << endl;
    delete tpl;
  }
  else if (ndims==2){
    TH2F* tpl=(TH2F*) data->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), false, nullptr) << endl;

    bool isAdded=false;
    for (auto it=procshape_2D.begin(); it!=procshape_2D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      cout << "\t- Creating " << tplname << endl;
      procshape_2D[procname]=TH2F(*tpl);
      procshape_2D[procname].SetName(tplname);
    }
    cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_2D[procname]), 1, procshape_2D[procname].GetNbinsX(), 1, procshape_2D[procname].GetNbinsY(), false, nullptr) << endl;
    delete tpl;
  }
  else if (ndims==1){
    TH1F* tpl=(TH1F*) data->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), false, nullptr) << endl;

    bool isAdded=false;
    for (auto it=procshape_1D.begin(); it!=procshape_1D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      cout << "\t- Creating " << tplname << endl;
      procshape_1D[procname]=TH1F(*tpl);
      procshape_1D[procname].SetName(tplname);
    }
    cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_1D[procname]), 1, procshape_1D[procname].GetNbinsX(), false, nullptr) << endl;
    delete tpl;
  }
  cout << "Work on data " << tplname << " is complete!" << endl;
  return ndims;
}


float getACMuF(TString cinputdir, float fai1){
  float res=1;
  if (cinputdir.Contains("/a3")) res = 1;
  else if (cinputdir.Contains("/a2")) res = sqrt((1.-fabs(fai1))*fabs(fai1))*TMath::Sign(1., fai1)*(1065.55457-2.*290.58626)/290.58626+1.;
  else if (cinputdir.Contains("/L1")) res = sqrt((1.-fabs(fai1))*fabs(fai1))*TMath::Sign(1., fai1)*(9.618383-2.*290.58626)/290.58626+1.;
  cout << "getACMuF result = " << 1./res << endl;
  return 1./res;
}
float getACMuV(TString cinputdir, float fai1){
  float res=getACMuF(cinputdir, fai1);
  float fai1_intcoef=sqrt((1.-fabs(fai1))*fabs(fai1))*TMath::Sign(1., fai1);
  float fai1_pureSMcoef=(1.-fabs(fai1));
  float fai1_pureBSMcoef=fabs(fai1);
  if (cinputdir.Contains("/a3")) res /= fai1_pureSMcoef + fai1_pureBSMcoef*pow(2.55052/0.297979018705, 2);
  else if (cinputdir.Contains("/a2")) res /= fai1_pureSMcoef + fai1_pureBSMcoef*pow(1.65684/0.27196538, 2) + fai1_intcoef*(2207.73/968.674-2.)*pow(1.65684/0.27196538, 1);
  else if (cinputdir.Contains("/L1")) res /= fai1_pureSMcoef + fai1_pureBSMcoef*pow(12100.42/2158.21307286, 2) + fai1_intcoef*(2861.213/968.674-2.)*pow(12100.42/2158.21307286, 1);
  cout << "getACMuV result = " << res << endl;
  return res;
}
void setControlVariableValue(std::unordered_map<TString, RooRealVar*> const& controlVars, TString const& key, double value){
  auto it_var = controlVars.find(key);
  if (it_var!=controlVars.end() && it_var->second) it_var->second->setVal(value);
}

TString getFractionString(float fai1){
  TString res = Form("%.5f", fai1);
  while (res.EndsWith("0")) res.Resize(res.Length()-1);
  if (res.EndsWith(".")) res.Resize(res.Length()-1);
  return res;
}


TGraphAsymmErrors* getDataGraph(TH1F* hdata){
  TGraphAsymmErrors* tgdata = nullptr;
  if (hdata){
    int ndata = 0;
    double xx_data[999];
    double xu_data[999];
    double xd_data[999];
    double yy_data[999];
    double yu_data[999];
    double yd_data[999];
    double rr_data[999];
    double ru_data[999];
    double rd_data[999];
    const double quant = (1.0 - 0.6827) / 2.0;
    double integral_data = 1;

    for (int bin = 1; bin <= hdata->GetNbinsX(); bin++){
      double bincenter = hdata->GetBinCenter(bin);
      double bincontent = hdata->GetBinContent(bin);

      if (bincontent > 0.){
        xx_data[ndata] = bincenter;
        yy_data[ndata] = bincontent / integral_data;
        xu_data[ndata] = 0;
        xd_data[ndata] = 0;
        yu_data[ndata] = (ROOT::Math::chisquared_quantile_c(quant, 2 * (bincontent + 1)) / 2. - bincontent) / integral_data;
        yd_data[ndata] = ((bincontent == 0) ? 0 : (bincontent - ROOT::Math::chisquared_quantile_c(1 - quant, 2 * bincontent) / 2.)) / integral_data;
        ndata++;
      }
    }
    cout << "Number of graph points: " << ndata << endl;
    tgdata = new TGraphAsymmErrors(ndata, xx_data, yy_data, xd_data, xu_data, yd_data, yu_data);
    tgdata->SetName("tgdata");
    tgdata->SetMarkerSize(1.2);
    tgdata->SetMarkerStyle(20);
    tgdata->SetMarkerColor(kBlack);
    tgdata->SetLineColor(kBlack);
    tgdata->SetLineWidth(1);
  }
  else cout << "Data histogram is null." << endl;
  //if (tgdata){
  //  cout << tgdata->GetName() << " graph is not null!" << endl;
  //  cout << "\t- Np = " << tgdata->GetN() << endl;
  //  for (int ip=0; ip<tgdata->GetN(); ip++) cout << "Point " << ip << ": ( " << tgdata->GetX()[ip] << ", " << tgdata->GetY()[ip] << " )" << endl;
  //}
  return tgdata;
}


bool FileExists(const char* fname){
  if (!fname) return false;
  struct stat sb;
  return (stat(fname, &sb) == 0 && S_ISREG(sb.st_mode));
}

void plotStacked(TString cinputdir, int onORoffshell=0, bool isEnriched=true, bool markPreliminary=false, bool isBlind=false, TString strFitResultFile=""){
  gStyle->SetOptStat(0);

  TString cinputdir_lower = cinputdir; cinputdir_lower.ToLower();
  bool const is_2l2nu = cinputdir_lower.Contains("2l2nu");
  bool const is_3l1nu = cinputdir_lower.Contains("3l1nu");
  bool const is_4l = !is_2l2nu && !is_3l1nu;

  std::vector<TString> const sqrtsnames{ "13TeV_2016", "13TeV_2017", "13TeV_2018" };
  const unsigned int nsqrts = sqrtsnames.size();
  std::vector<TString> catnames;
  if (is_4l) catnames = std::vector<TString>{ "Untagged", "JJVBFTagged", "HadVHTagged" };
  else if (is_2l2nu) catnames = std::vector<TString>{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2_pTmiss_lt_200", "Nj_geq_2_pTmiss_ge_200" };
  else if (is_3l1nu) catnames = std::vector<TString>{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  else return;
  const unsigned int ncats = catnames.size();

  if (onORoffshell==0 && !is_4l){
    cerr << "\t- Cannot have an on-shell category for non-4l final states! Skipping..." << endl;
    return;
  }

  std::vector<TString> channames;
  if (is_2l2nu) channames = std::vector<TString>{ "2e2nu", "2mu2nu" };
  else if (is_3l1nu) channames = std::vector<TString>{ "2l1e", "2l1mu" };
  else /*if (is_4l)*/ channames = std::vector<TString>{ "4mu", "4e", "2e2mu" };
  const unsigned int nchans = channames.size();

  bool const isPostFit = (strFitResultFile!="");

  TDirectory* curdir = gDirectory;

  // Determine alternative model parameters
  TString failabel="f_{ai}";
  TString ailabel="";
  TString aiKDlabel="";
  TString aihypo="";
  if (cinputdir.Contains("/a3")){ failabel="f_{a3}"; ailabel="a_{3}"; aiKDlabel=aihypo="a3"; }
  else if (cinputdir.Contains("/a2")){ failabel="f_{a2}"; ailabel="a_{2}"; aiKDlabel=aihypo="a2"; }
  else if (cinputdir.Contains("/L1")){ failabel="f_{#Lambda1}"; ailabel="#Lambda_{1}"; aiKDlabel="#Lambda1"; aihypo="L1"; }
  else if (cinputdir.Contains("/SM")){ aiKDlabel="a2"; }

  // Determine alternative model parameters
  constexpr double GHref = 4.07;
  unordered_map<TString, double> val_fai1ALT; val_fai1ALT["RV"]=1; val_fai1ALT["RF"]=1; val_fai1ALT["fai1"]=0;
  unordered_map<TString, double> val_GGsmALT; val_GGsmALT["RV"]=1; val_GGsmALT["RF"]=1; val_GGsmALT["GGsm"]=1;
  /*
  if (onORoffshell==1){ // Offshell
    if (cinputdir.Contains("/a3")){ val_fai1ALT["RV"]=0.22; val_fai1ALT["RF"]=1.06; val_fai1ALT["fai1"]=0.01; }
    else if (cinputdir.Contains("/a2")){ val_fai1ALT["RV"]=0.3; val_fai1ALT["RF"]=1.3; val_fai1ALT["fai1"]=-0.015; }
    else if (cinputdir.Contains("/L1")){ val_fai1ALT["RV"]=0.1; val_fai1ALT["RF"]=0.85; val_fai1ALT["fai1"]=0.015; }
    val_GGsmALT["RV"]=0.29; val_GGsmALT["RF"]=0.88; val_GGsmALT["GGsm"]=14./GHref;
  }
  else{ // Onshell
    if (cinputdir.Contains("/a3")){ val_fai1ALT["RV"]=0.22; val_fai1ALT["RF"]=1.06; val_fai1ALT["fai1"]=0.01; }
    else if (cinputdir.Contains("/a2")){ val_fai1ALT["RV"]=0.3; val_fai1ALT["RF"]=1.3; val_fai1ALT["fai1"]=-0.015; }
    else if (cinputdir.Contains("/L1")){ val_fai1ALT["RV"]=0.1; val_fai1ALT["RF"]=0.85; val_fai1ALT["fai1"]=0.015; }
    else{ val_GGsmALT["RV"]=0.29; val_GGsmALT["RF"]=0.88; val_GGsmALT["GGsm"]=1; }
  }
  */
  if (onORoffshell==1){ // Offshell
    if (is_4l || is_2l2nu){
      if (cinputdir.Contains("/a3")){ val_fai1ALT["fai1"]=0.05; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/a2")){ val_fai1ALT["fai1"]=0.05; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/L1")){ val_fai1ALT["fai1"]=0.05; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      val_GGsmALT["RV"]=1; val_GGsmALT["RF"]=1; val_GGsmALT["GGsm"]=20./GHref;
    }
    else if (is_3l1nu){
      if (cinputdir.Contains("/a3")){ val_fai1ALT["fai1"]=1; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/a2")){ val_fai1ALT["fai1"]=1; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/L1")){ val_fai1ALT["fai1"]=1; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      val_GGsmALT["RV"]=1; val_GGsmALT["RF"]=1; val_GGsmALT["GGsm"]=2000./GHref;
    }
  }
  else{ // Onshell
    if (cinputdir.Contains("/a3")){ val_fai1ALT["fai1"]=0.5; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
    else if (cinputdir.Contains("/a2")){ val_fai1ALT["fai1"]=-0.5; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
    else if (cinputdir.Contains("/L1")){ val_fai1ALT["fai1"]=0.5; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
    else{ val_GGsmALT["RV"]=0.29; val_GGsmALT["RF"]=0.88; val_GGsmALT["GGsm"]=1; }
  }

  std::unordered_map<TString, double> postfit_vals_map;
  if (isPostFit){
    TFile* finput_postfit = TFile::Open(strFitResultFile, "read");
    RooFitResult* fitResult = dynamic_cast<RooFitResult*>(finput_postfit->Get("fit_mdf"));
    RooArgList const& finalFloatPars = fitResult->floatParsFinal();
    TIterator* it = nullptr;
    it = finalFloatPars.createIterator();
    RooAbsArg* var;
    while ((var = (RooAbsArg*) it->Next())){
      RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
      if (rvar){
        TString strvar = rvar->GetName();
        double val = rvar->getVal();
        postfit_vals_map[strvar] = val;
      }
    }
    delete it;
    finput_postfit->Close();
    curdir->cd();
  }

  for (unsigned int icat=0; icat<ncats; icat++){
    TString const& catname = catnames.at(icat);
    cout << "Acquiring shapes for category " << catname << endl;

    unsigned int ndims=0;
    unordered_map<TString, TH1F> procshape_1D;
    unordered_map<TString, TH2F> procshape_2D;
    unordered_map<TString, TH3F> procshape_3D;

    // Get process order/labels
    std::vector<TString> proc_order, proc_label;
    std::vector<int> proc_color;
    std::vector<int> proc_code;
    if (is_4l){
      if (onORoffshell==1){ // Offshell
        proc_order=std::vector<TString>{ "Zjets", "qqZZ", "VVZZ_offshell", "ggZZ_offshell", "total_GGsmALT" };
        if (!cinputdir.Contains("/SM")) proc_order.push_back("total_fai1ALT");
      }
      else{
        proc_order=std::vector<TString>{ "zjets", "bkg_zz", "VVZZ", "ggH", "total_fai1ALT" };
      }
    }
    else if (is_2l2nu){
      proc_order=std::vector<TString>{ "tZX", "qqZZ_offshell", "qqWZ_offshell", "NRB_2l2nu", "InstrMET", "VVVV_onshell", "VVVV_offshell", "ggZZ_offshell", "total_GGsmALT" };
      if (!cinputdir.Contains("/SM")) proc_order.push_back("total_fai1ALT");
    }
    else if (is_3l1nu){
      proc_order=std::vector<TString>{ "tWX", "tZX", "ttbar_2l2nu", "qqZG", "DY_2l", "qqZZ", "qqWZ", "VVVV_onshell", "VVVV_offshell", "total_GGsmALT" };
      if (!cinputdir.Contains("/SM")) proc_order.push_back("total_fai1ALT");
    }
    proc_order.push_back("data");
    for (auto const& p:proc_order){
      if (is_4l){
        if (p=="bkg_qqzz" || p=="qqZZ"){ proc_color.push_back(int(TColor::GetColor("#99ccff"))); proc_label.push_back("q#bar{q}#rightarrow4l bkg."); proc_code.push_back(0); }
        else if (p=="bkg_gg"){ proc_color.push_back(int(kBlue)); proc_label.push_back("gg#rightarrow4l bkg."); proc_code.push_back(0); }
        else if (p=="zjets" || p=="Zjets"){ proc_color.push_back(int(TColor::GetColor("#669966"))); proc_label.push_back("Z+X"); proc_code.push_back(0); }
        else if (p=="bkg_vv"){ proc_color.push_back(int(kPink+9)); proc_label.push_back("EW bkg."); proc_code.push_back(0); }
        else if (p=="bkg_zz"){ proc_color.push_back(int(TColor::GetColor("#99ccff"))); proc_label.push_back("ZZ/Z#gamma*/#gamma*#gamma*#rightarrow4l bkg."); proc_code.push_back(0); }
        else if (p=="ggZZ_offshell" || p=="ggH"){
          //proc_color.push_back(int(kOrange-2));
          proc_color.push_back(int(TColor::GetColor("#ffdcdc")));
          if (p=="ggH"){ proc_label.push_back("Total SM"); proc_code.push_back(1); }
          else if (p=="ggZZ_offshell"){ proc_label.push_back("gg#rightarrow4l SM s+b+i"); proc_code.push_back(2); }
        }
        else if (p=="VVZZ" || p=="VVZZ_offshell"){
          proc_color.push_back(int(TColor::GetColor("#ff9b9b")));
          if (p=="VVZZ"){ proc_label.push_back("VBF+VH SM"); proc_code.push_back(1); }
          else if (p=="VVZZ_offshell"){ proc_label.push_back("EW SM s+b+i"); proc_code.push_back(2); }
          //if (p=="VVZZ"){ proc_label.push_back("EW sig."); proc_code.push_back(1); }
          //else if (p=="VVZZ_offshell"){ proc_label.push_back("EW SM total"); proc_code.push_back(2); }
        }
        else if (p.Contains("ALT")){
          if (onORoffshell==1){ // Offshell
            if (p.Contains("GGsmALT")){ proc_label.push_back(Form("Total (f_{ai}=0, #Gamma_{H}=%s MeV)", getFractionString(val_GGsmALT["GGsm"]*GHref).Data())); proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2); }
            else if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total (%s=%s, #Gamma_{H}=#Gamma_{H}^{SM})", failabel.Data(), getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-2); }
          }
          else{ // Onshell
            if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total, %s=%s", failabel.Data(), getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-1); }
          }
        }
        else if (p=="data"){
          proc_color.push_back(int(kBlack));
          proc_code.push_back(-99);
          proc_label.push_back("Observed");
        }
        else cerr << p << " is unmatched for labeling!" << endl;
      }
      else if (is_2l2nu){
        if (p=="qqZZ_offshell"){ proc_color.push_back(int(TColor::GetColor("#99ccff"))); proc_label.push_back("q#bar{q}#rightarrowZZ bkg."); proc_code.push_back(0); }
        else if (p=="qqWZ_offshell"){ proc_color.push_back(int(kBlue)); proc_label.push_back("q#bar{q}#rightarrowWZ bkg."); proc_code.push_back(0); }
        else if (p=="InstrMET"){ proc_color.push_back(int(TColor::GetColor("#669966"))); proc_label.push_back("Instr. p_{T}^{miss}"); proc_code.push_back(0); }
        else if (p=="NRB_2l2nu"){ proc_color.push_back(int(kGray+1)); proc_label.push_back("Nonresonant bkg."); proc_code.push_back(0); }
        else if (p=="tZX"){ proc_color.push_back(int(kOrange-6)); proc_label.push_back("tZ+X bkg."); proc_code.push_back(0); }
        else if (p=="ggZZ_offshell"){
          proc_color.push_back(int(TColor::GetColor("#ffdcdc")));
          proc_label.push_back("gg SM total"); proc_code.push_back(2);
        }
        else if (p=="VVVV_offshell"){
          proc_color.push_back(int(TColor::GetColor("#ff9b9b")));
          proc_label.push_back("EW SM total (off-shell)"); proc_code.push_back(1);
        }
        else if (p=="VVVV_onshell"){
          proc_color.push_back(int(kPink+9));
          proc_label.push_back("EW SM total (on-shell)"); proc_code.push_back(1);
        }
        else if (p.Contains("ALT")){
          if (onORoffshell==1){ // Offshell
            if (p.Contains("GGsmALT")){ proc_label.push_back(Form("Total (f_{ai}=0, #Gamma_{H}=%s MeV)", getFractionString(val_GGsmALT["GGsm"]*GHref).Data())); proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2); }
            else if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total (f_{ai}=%s, #Gamma_{H}=#Gamma_{H}^{SM})", getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-2); }
          }
          else{ // Onshell
            if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total, %s=%s", failabel.Data(), getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-1); }
          }
        }
        else if (p=="data"){
          proc_color.push_back(int(kBlack));
          proc_code.push_back(-99);
          proc_label.push_back("Observed");
        }
        else cerr << p << " is unmatched for labeling!" << endl;
      }
      else if (is_3l1nu){
        if (p=="qqZZ"){ proc_color.push_back(int(TColor::GetColor("#99ccff"))); proc_label.push_back("q#bar{q}#rightarrowZZ bkg."); proc_code.push_back(0); }
        else if (p=="qqWZ"){ proc_color.push_back(int(kBlue)); proc_label.push_back("q#bar{q}#rightarrowWZ bkg."); proc_code.push_back(0); }
        else if (p=="qqZG"){ proc_color.push_back(int(TColor::GetColor("#f1c232"))); proc_label.push_back("q#bar{q}#rightarrowZ#gamma bkg."); proc_code.push_back(0); }
        else if (p=="DY_2l"){ proc_color.push_back(int(TColor::GetColor("#669966"))); proc_label.push_back("Drell-Yan bkg."); proc_code.push_back(0); }
        else if (p=="ttbar_2l2nu"){ proc_color.push_back(int(kGray+1)); proc_label.push_back("t#bar{t} bkg."); proc_code.push_back(0); }
        else if (p=="tZX"){ proc_color.push_back(int(kOrange-6)); proc_label.push_back("tZ+X bkg."); proc_code.push_back(0); }
        else if (p=="tWX"){ proc_color.push_back(int(TColor::GetColor("#674ea7"))); proc_label.push_back("tW+X bkg."); proc_code.push_back(0); }
        else if (p=="VVVV_offshell"){
          proc_color.push_back(int(TColor::GetColor("#ff9b9b")));
          proc_label.push_back("EW SM total (off-shell)"); proc_code.push_back(1);
        }
        else if (p=="VVVV_onshell"){
          proc_color.push_back(int(kPink+9));
          proc_label.push_back("EW SM total (on-shell)"); proc_code.push_back(1);
        }
        else if (p.Contains("ALT")){
          if (onORoffshell==1){ // Offshell
            if (p.Contains("GGsmALT")){ proc_label.push_back(Form("Total (f_{ai}=0, #Gamma_{H}=%s GeV)", getFractionString(val_GGsmALT["GGsm"]*GHref/1000.).Data())); proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2); }
            else if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total (f_{ai}=%s, #Gamma_{H}=#Gamma_{H}^{SM})", getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-2); }
          }
          else{ // Onshell
            if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total, %s=%s", failabel.Data(), getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-1); }
          }
        }
        else if (p=="data"){
          proc_color.push_back(int(kBlack));
          proc_code.push_back(-99);
          proc_label.push_back("Observed");
        }
        else cerr << p << " is unmatched for labeling!" << endl;
      }
    }

    for (unsigned int ic=0; ic<nchans; ic++){
      TString const& channame = channames.at(ic);
      cout << "\t- Channel:  " << channame << endl;
      for (unsigned int is=0; is<nsqrts; is++){
        TString const& sqrtsname = sqrtsnames.at(is);
        cout << "\t- Data period:  " << sqrtsname << endl;
        TString cinput_main = cinputdir + "/" + (onORoffshell ? "Offshell_" : "Onshell_") + sqrtsname + "/hto" + channame + "_" + catname;
        TString cinput_file = cinput_main + ".input.root";
        TString cinput_dc = cinput_main + ".txt";
        if (!FileExists(cinput_file)){ cout << "File " << cinput_file << " does not exist. Skipping this channel..." << endl; continue; }
        else if (!FileExists(cinput_dc)){ cout << "File " << cinput_dc << " does not exist. Skipping this channel..." << endl; continue; }

        vector<TString> procname;
        vector<double> procrate;
        vector<RooRealVar*> lnNmodvars;
        vector<RooConstVar*> kappavars;
        unordered_map<TString, std::vector<RooAbsReal*>> procname_rateModifier_map;

        // Get process names
        ifstream tin(cinput_dc.Data());
        char line[512];
        while (string(line).find("process")==string::npos) tin.getline(line, 512);
        char* chars_array = strtok(line, " ");
        chars_array = strtok(NULL, " ");
        while (chars_array){
          procname.push_back(TString(chars_array));
          chars_array = strtok(NULL, " ");
          procname_rateModifier_map[procname.back()] = std::vector<RooAbsReal*>();
        }
        // Get process rates
        while (string(line).find("rate")==string::npos) tin.getline(line, 512);
        chars_array = strtok(line, " ");
        chars_array = strtok(NULL, " ");
        while (chars_array){
          procrate.push_back(atof(chars_array));
          chars_array = strtok(NULL, " ");
        }
        // Get lnN systematics
        while (!tin.eof()){
          tin.getline(line, 512);

          string strline = line;
          std::replace(strline.begin(), strline.end(), ',', ' ');
          strline.erase(std::remove(strline.begin(), strline.end(), '['), strline.end());
          strline.erase(std::remove(strline.begin(), strline.end(), ']'), strline.end());

          bool isLog = strline.find("lnN")!=string::npos;
          if (isLog){
            vector<string> systdist;
            splitOptionRecursive(strline, systdist, ' ');
            string systname = systdist.at(0);
            string systtype = systdist.at(1);
            string accumulate="";
            cout << "Processing systematic " << systname << endl;
            RooRealVar* systvar = new RooRealVar(systname.data(), "", 0, -7, 7); lnNmodvars.push_back(systvar);
            for (unsigned int ip=0; ip<procname.size(); ip++){
              string systline = systdist.at(ip+2);
              TString const& pname = procname.at(ip);
              auto& rateModifiers = procname_rateModifier_map.find(pname)->second;
              if (systline.find("-")==string::npos && systline!=""){
                double kdn=1, kup=1;
                if (systline.find("/")!=std::string::npos){
                  string ssf, ssl;
                  splitOption(systline, ssf, ssl, '/');
                  kdn = stod(ssf);
                  kup = stod(ssl);
                }
                else{
                  kup = stod(systline);
                  kdn = 1./kup;
                }
                RooConstVar* var_kdn = new RooConstVar(Form("%s_%s_kdn", systname.data(), pname.Data()), "", kdn); kappavars.push_back(var_kdn);
                RooConstVar* var_kup = new RooConstVar(Form("%s_%s_kup", systname.data(), pname.Data()), "", kup); kappavars.push_back(var_kup);
                AsymPow* tmpvar = new AsymPow(Form("%s_%s", systname.data(), pname.Data()), "", *var_kdn, *var_kup, *systvar);
                rateModifiers.push_back((RooAbsReal*) tmpvar);
              }
            }
          }
        }
        tin.close();
        assert(procrate.size()==procname.size());

        // Construct the processes
        unordered_map<TString, process_spec> procSpecs;
        {
          std::vector<double>::iterator it_rate=procrate.begin();
          for (TString const& pname:procname){
            cout << "Process " << pname << " is found in the datacard" << endl;
            procSpecs[pname] = process_spec(pname, *it_rate);
            auto const& rateModifiers = procname_rateModifier_map.find(pname)->second;
            procSpecs[pname].setRateModifiers(rateModifiers);
            it_rate++;
          }
        }

        // Open the input
        TFile* finput = TFile::Open(cinput_file, "read");
        RooWorkspace* ws = (RooWorkspace*) finput->Get("w");

        // Set postfit values first if they are present.
        for (auto const& pp:postfit_vals_map){
          RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
          if (tmpvar){
            cout << "\t- Setting " << pp.first << " = " << pp.second << " in the workspace." << endl;
            tmpvar->setVal(pp.second);
          }

          for (auto const& lnNmodvar:lnNmodvars){
            if (TString(lnNmodvar->GetName()) == pp.first){
              cout << "\t- Setting " << pp.first << " = " << pp.second << " in the list of lnN nuisances." << endl;
              lnNmodvar->setVal(pp.second);
            }
          }
        }
        unordered_map<TString, RooRealVar*> controlVars;
        TString varsToCheck[10]={
          "R", "RF", "RF_13TeV", "RV", "RV_13TeV", "R_13TeV", "CMS_zz4l_fai1", "GGsm", "kbkg_gg", "kbkg_VBF"
        };
        for (unsigned int v=0; v<10; v++){
          controlVars[varsToCheck[v]] = ws->var(varsToCheck[v]);
          if (controlVars[varsToCheck[v]]){
            if (varsToCheck[v]!="CMS_zz4l_fai1"){
              controlVars[varsToCheck[v]]->setVal(1);
              if (varsToCheck[v]=="GGsm") controlVars[varsToCheck[v]]->removeRange();
            }
            else controlVars[varsToCheck[v]]->setVal(0);
          }
          else cerr << varsToCheck[v] << " could not be found!" << endl;
        }
        RooDataSet* data = (RooDataSet*) ws->data("data_obs");

        // Get process pdfs
        for (auto const& pname:procname){
          if (pname=="data") continue;
          curdir->cd();
          cout << "Extracting the pdf and norm from process " << pname << endl;
          RooAbsPdf* pdf = ws->pdf(pname);
          RooAbsReal* norm = (RooAbsReal*) ws->factory(pname+"_norm");
          if (!norm) cout << "Warning: " << pname << "_norm is not found." << endl;
          if (!pdf) cerr << pname << " pdf could not be found." << endl;
          procSpecs[pname].pdf=pdf;
          procSpecs[pname].norm=norm;

          if (!cinputdir.Contains("/SM") && !onORoffshell){ // AC onshell
            if (pname.Contains("ggZZ")){
              controlVars["kbkg_gg"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "ggH");
              controlVars["kbkg_gg"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_gg");
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              controlVars["kbkg_gg"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
              controlVars["kbkg_gg"]->setVal(1);
            }
            else if (pname.Contains("ttH") || pname.Contains("bbH")){
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "ggH");

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RF"]->setVal(1);
            }
            else if (pname.Contains("VBF")){
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVZZ");
              controlVars["kbkg_VBF"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_vv");
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
              controlVars["kbkg_VBF"]->setVal(1);
            }
            else if (pname.Contains("ZH") || pname.Contains("WH")){
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVZZ");
              controlVars["kbkg_VBF"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_vv");
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
              controlVars["kbkg_VBF"]->setVal(1);
            }
            else if (pname.Contains("bkg_qqzz")){
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
            }
          }
          else if (!onORoffshell && (pname.Contains("ttH") || pname.Contains("bbH"))){ // On-shell SM
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "ggH");
          }
          else if (!onORoffshell && (pname.Contains("ZH") || pname.Contains("WH") || pname.Contains("VBF"))){ // On-shell SM
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVZZ");
          }
          //else if (!onORoffshell && (pname.Contains("bkg_zzz") || pname.Contains("bkg_wzz") || pname.Contains("bkg_vbs"))){ // On-shell SM
          else if (!onORoffshell && (pname.Contains("bkg_zzz") || pname.Contains("bkg_wzz") || pname.Contains("bkg_vbs") || pname.Contains("bkg_gg") || pname.Contains("bkg_qqzz"))){ // On-shell SM
            //ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_vv");
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
          }
          else if (onORoffshell && (pname.Contains("ggZZ_offshell") || pname.Contains("VVZZ_offshell") || pname.Contains("VVVV_offshell"))){ // Off-shell process
            setControlVariableValue(controlVars, "GGsm", val_GGsmALT["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_GGsmALT["RV"]);
            setControlVariableValue(controlVars, "RF", val_GGsmALT["RF"]);
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_GGsmALT");
            setControlVariableValue(controlVars, "GGsm", 1);
            setControlVariableValue(controlVars, "RV", 1);
            setControlVariableValue(controlVars, "RF", 1);

            if (!cinputdir.Contains("/SM")){
              setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_fai1ALT["CMS_zz4l_fai1"]);
              setControlVariableValue(controlVars, "RV", val_fai1ALT["RV"]);
              setControlVariableValue(controlVars, "RF", val_fai1ALT["RF"]);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
              setControlVariableValue(controlVars, "CMS_zz4l_fai1", 0);
              setControlVariableValue(controlVars, "RV", 1);
              setControlVariableValue(controlVars, "RF", 1);
            }
          }
          else if (onORoffshell && (pname.Contains("ggZZ_onshell") || pname.Contains("VVZZ_onshell") || pname.Contains("VVVV_onshell"))){ // On-shell process in off-shell data cards
            setControlVariableValue(controlVars, "RV", val_GGsmALT["RV"]);
            setControlVariableValue(controlVars, "RF", val_GGsmALT["RF"]);
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_GGsmALT");
            setControlVariableValue(controlVars, "RV", 1);
            setControlVariableValue(controlVars, "RF", 1);

            if (!cinputdir.Contains("/SM")){
              setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_fai1ALT["CMS_zz4l_fai1"]);
              setControlVariableValue(controlVars, "RV", val_fai1ALT["RV"]);
              setControlVariableValue(controlVars, "RF", val_fai1ALT["RF"]);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
              setControlVariableValue(controlVars, "CMS_zz4l_fai1", 0);
              setControlVariableValue(controlVars, "RV", 1);
              setControlVariableValue(controlVars, "RF", 1);
            }
          }

          ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D);
        }
        extractDataTemplates(procSpecs[procname.front()], data, procshape_1D, procshape_2D, procshape_3D, "data");

        finput->Close();
        curdir->cd();

        // Clean up rate modifiers
        for (auto& pp:procname_rateModifier_map){
          for (auto& ppp:pp.second) delete ppp;
        }
        for (auto& vv:kappavars) delete vv;
        for (auto& vv:lnNmodvars) delete vv;
      }
    }

    vector<TString> varlabels;
    // KD1
    if (ndims>=1){
      if (is_4l){
        if (onORoffshell==1 || cinputdir.Contains("/SM")) varlabels.push_back("m_{4l} (GeV)");
        else{ // On-shell AC
          if (catname=="Untagged") varlabels.push_back("D_{bkg}");
          else if (catname=="JJVBFTagged") varlabels.push_back("D_{bkg,m4l}^{VBF+dec}");
          else if (catname=="HadVHTagged") varlabels.push_back("D_{bkg,m4l}^{VH+dec}");
        }
      }
      else if (is_2l2nu){
        varlabels.push_back("m_{T}^{ZZ} (GeV)");
      }
      else if (is_3l1nu){
        varlabels.push_back("m_{T}^{WZ} (GeV)");
      }
    }
    // KD2
    if (ndims>=2){
      if (is_4l){
        if (onORoffshell==1 || cinputdir.Contains("/SM")){
          if (catname=="Untagged") varlabels.push_back("D_{bkg}^{kin}");
          else if (catname=="JJVBFTagged") varlabels.push_back("D_{bkg}^{VBF+dec}");
          else if (catname=="HadVHTagged") varlabels.push_back("D_{bkg}^{VH+dec}");
        }
        else{ // On-shell AC
          if (catname=="Untagged"){
            if (cinputdir.Contains("/a3")) varlabels.push_back("D_{0-}^{dec}");
            else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{0h+}^{dec}");
            else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{#Lambda1}^{dec}");
          }
          else if (catname=="JJVBFTagged"){
            if (cinputdir.Contains("/a3")) varlabels.push_back("D_{0-}^{VBF+dec}");
            else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{0h+}^{VBF+dec}");
            else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{#Lambda1}^{VBF+dec}");
          }
          else if (catname=="HadVHTagged"){
            if (cinputdir.Contains("/a3")) varlabels.push_back("D_{0-}^{VH+dec}");
            else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{0h+}^{VH+dec}");
            else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{#Lambda1}^{VH+dec}");
          }
        }
      }
      else if (is_2l2nu){
        varlabels.push_back(((catname.Contains("Nj_geq_2") && is_2l2nu) ? "D_{2jet}^{VBF}" : "p_{T}^{miss} (GeV)"));
      }
      else if (is_3l1nu){
        cerr << "Second dimension is not implemented in 3l1nu." << endl;
        exit(1);
      }
    }
    // KD3
    if (ndims>=3){
      if (is_4l){
        if (onORoffshell==1){
          if (cinputdir.Contains("/SM")){
            if (catname=="Untagged") varlabels.push_back("D_{bsi}^{gg,dec}");
            else if (catname=="JJVBFTagged") varlabels.push_back("D_{bsi}^{VBF+dec}");
            else if (catname=="HadVHTagged") varlabels.push_back("D_{bsi}^{VH+dec}");
          }
          else{
            if (catname=="Untagged"){
              if (cinputdir.Contains("/a3")) varlabels.push_back("D_{0-}^{dec}");
              else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{0h+}^{dec}");
              else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{#Lambda1}^{dec}");
            }
            else if (catname=="JJVBFTagged"){
              if (cinputdir.Contains("/a3")) varlabels.push_back("D_{0-}^{VBF+dec}");
              else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{0h+}^{VBF+dec}");
              else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{#Lambda1}^{VBF+dec}");
            }
            else if (catname=="HadVHTagged"){
              if (cinputdir.Contains("/a3")) varlabels.push_back("D_{0-}^{VH+dec}");
              else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{0h+}^{VH+dec}");
              else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{#Lambda1}^{VH+dec}");
            }
          }
        }
        else{ // On-shell AC
          if (catname=="Untagged"){
            if (cinputdir.Contains("/a3")) varlabels.push_back("D_{CP}^{dec}");
            else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{int}^{dec}");
            else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{0h+}^{dec}");
          }
          else if (catname=="JJVBFTagged"){
            if (cinputdir.Contains("/a3")) varlabels.push_back("D_{CP}^{VBF}");
            else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{int}^{VBF}");
            else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{0h+}^{VBF+dec}");
          }
          else if (catname=="HadVHTagged"){
            if (cinputdir.Contains("/a3")) varlabels.push_back("D_{CP}^{VH}");
            else if (cinputdir.Contains("/a2")) varlabels.push_back("D_{int}^{VH}");
            else if (cinputdir.Contains("/L1")) varlabels.push_back("D_{0h+}^{VH+dec}");
          }
        }
      }
      else if (is_2l2nu){
        varlabels.push_back(Form("D_{2jet}^{VBF,%s}", aiKDlabel.Data()));
      }
      else if (is_3l1nu){
        cerr << "Third dimension is not implemented in 3l1nu." << endl;
        exit(1);
      }
    }

    if (ndims==0) continue;
    cout << "\t- Variable labels: ";
    for (auto& vl:varlabels) cout << vl << " ";
    cout << endl;

    curdir->cd();
    int massDim=-1;
    int kdDim=-1;
    for (unsigned int idim=0; idim<ndims; idim++){
      if (varlabels.at(idim).Contains("m_{4l}") || varlabels.at(idim).Contains("m_{T}^{ZZ}") || varlabels.at(idim).Contains("m_{T}^{WZ}")) massDim=idim;
      else if (
        (is_4l && varlabels.at(idim).Contains("D_{bkg}"))
        ||
        (
          is_2l2nu
          && 
          (varlabels.at(idim).Contains("p_{T}^{miss}") || varlabels.at(idim)=="D_{2jet}^{VBF}")
          )
        ) kdDim=idim;
    }
    cout << "Mass dim: " << massDim << endl;
    cout << "KD dim: " << kdDim << endl;

    std::vector< std::vector<float> > cutvals(ndims, std::vector<float>(ndims, -1));
    unordered_map<TString, std::vector<TH1F*>> procdist; // The key is for the process label, and the vector is for the dimensions
    for (unsigned int idim=0; idim<ndims; idim++){
      TString dimname = Form("_KD%i", idim+1);
      cout << "Preparing projections on dimension " << idim << "..." << endl;
      if (ndims==3){
        for (auto it=procshape_3D.begin(); it!=procshape_3D.end(); it++){
          TH3F* hist = &(it->second);
          TAxis* yaxis=nullptr;
          TAxis* zaxis=nullptr;
          switch (idim){
          case 0:
            yaxis=hist->GetYaxis();
            zaxis=hist->GetZaxis();
            break;
          case 1:
            yaxis=hist->GetZaxis();
            zaxis=hist->GetXaxis();
            break;
          case 2:
            yaxis=hist->GetXaxis();
            zaxis=hist->GetYaxis();
            break;
          }
          int iy=1, iz=1;
          int jy=yaxis->GetNbins(), jz=zaxis->GetNbins();
          if (isEnriched){
            const float valMassCut=(is_4l ? 340. : (is_2l2nu ? 450. : (is_3l1nu ? 400. : -1.)));
            if (onORoffshell==1 && (int) idim!=massDim){ // Cut on mass
              if (massDim==0){
                if (idim==1){ iz=zaxis->FindBin(valMassCut); }
                else if (idim==2){ iy=yaxis->FindBin(valMassCut); }
              }
              else if (massDim==1){
                if (idim==0){ iy=yaxis->FindBin(valMassCut); }
                else if (idim==2){ iz=zaxis->FindBin(valMassCut); }
              }
              else if (massDim==2){
                if (idim==0){ iz=zaxis->FindBin(valMassCut); }
                else if (idim==1){ iy=yaxis->FindBin(valMassCut); }
              }
              cutvals.at(idim).at(massDim) = valMassCut;
            }

            const float valKDCut=(onORoffshell==1 ? (is_2l2nu ? 0.8 : 0.6) : 0.5);
            if ((int) idim!=kdDim){ // Cut on D_bkg
              if (kdDim==0){
                if (idim==1){ iz=zaxis->FindBin(valKDCut); }
                else if (idim==2){ iy=yaxis->FindBin(valKDCut); }
              }
              else if (kdDim==1){
                if (idim==0){ iy=yaxis->FindBin(valKDCut); }
                else if (idim==2){ iz=zaxis->FindBin(valKDCut); }
              }
              else if (kdDim==2){
                if (idim==0){ iz=zaxis->FindBin(valKDCut); }
                else if (idim==1){ iy=yaxis->FindBin(valKDCut); }
              }
              cutvals.at(idim).at(kdDim) = valKDCut;
            }
          }
          cout << "Cutting y-axis in range " << iy << ", " << jy << endl;
          cout << "Cutting z-axis in range " << iz << ", " << jz << endl;
          TH1F* htmp=getHistogramSlice(
            hist, idim,
            iy, jy, iz, jz,
            it->first + dimname + "_hist"
          );
          htmp->GetXaxis()->SetTitle(varlabels.at(idim));
          cout << "\t- Constructed histogram " << htmp->GetName() << endl;
          cout << "\t\t- Initial histogram integral: " << getHistogramIntegralAndError(hist, 1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1, hist->GetNbinsZ(), false, nullptr) << endl;
          if (procdist.find(it->first)==procdist.end()) procdist[it->first] = std::vector<TH1F*>();
          procdist[it->first].push_back(htmp);
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first].back(), 1, procdist[it->first].back()->GetNbinsX(), false, nullptr) << endl;
        }
      }
      else if (ndims==2){
        for (auto it=procshape_2D.begin(); it!=procshape_2D.end(); it++){
          TH2F* hist = &(it->second);
          TAxis* yaxis=nullptr;
          switch (idim){
          case 0:
            yaxis=hist->GetYaxis();
            break;
          case 1:
            yaxis=hist->GetXaxis();
            break;
          }
          int iy=1;
          int jy=yaxis->GetNbins();
          if (isEnriched){
            const float valMassCut=(is_4l ? 340. : (is_2l2nu ? 450. : (is_3l1nu ? 400. : -1.)));
            if (onORoffshell==1 && (int) idim!=massDim){ // Cut on mass
              iy = yaxis->FindBin(valMassCut);
              cutvals.at(idim).at(massDim) = valMassCut;
            }
            const float valKDCut=(onORoffshell==1 ? (varlabels.at(kdDim).Contains("p_{T}^{miss}") ? 200. : (is_2l2nu ? 0.8 : 0.6)) : 0.5);
            if ((int) idim!=kdDim){ // Cut on D_bkg
              iy = yaxis->FindBin(valKDCut);
              cutvals.at(idim).at(kdDim) = valKDCut;
            }
          }
          cout << "Cutting y-axis in range " << iy << ", " << jy << endl;
          TH1F* htmp=getHistogramSlice(
            hist, idim,
            iy, jy,
            it->first + dimname + "_hist"
          );
          htmp->GetXaxis()->SetTitle(varlabels.at(idim));
          cout << "\t- Constructed histogram " << htmp->GetName() << endl;
          cout << "\t\t- Initial histogram integral: " << getHistogramIntegralAndError(hist, 1, hist->GetNbinsX(), 1, hist->GetNbinsY(), false, nullptr) << endl;
          if (procdist.find(it->first)==procdist.end()) procdist[it->first] = std::vector<TH1F*>();
          procdist[it->first].push_back(htmp);
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first].back(), 1, procdist[it->first].back()->GetNbinsX(), false, nullptr) << endl;
        }
      }
      else{
        for (auto it=procshape_1D.begin(); it!=procshape_1D.end(); it++){
          TH1F* hist = &(it->second);
          TH1F* htmp = getHistogramSlice(
            hist, it->first + dimname + "_hist"
          );
          htmp->GetXaxis()->SetTitle(varlabels.at(idim));
          cout << "\t- Constructed histogram " << htmp->GetName() << endl;
          cout << "\t\t- Initial histogram integral: " << getHistogramIntegralAndError(hist, 1, hist->GetNbinsX(), false, nullptr) << endl;
          if (procdist.find(it->first)==procdist.end()) procdist[it->first] = std::vector<TH1F*>();
          procdist[it->first].push_back(htmp);
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first].back(), 1, procdist[it->first].back()->GetNbinsX(), false, nullptr) << endl;
        }
      }
    }

    for (auto& pp:procdist) cout << "Process " << pp.first << " has " << pp.second.size() << " projections." << endl;
    for (unsigned short idim=0; idim<cutvals.size(); idim++){
      cout << "Cut values defined for the projection on dimension " << idim << ": ";
      for (unsigned short jdim=0; jdim<cutvals.at(idim).size(); jdim++) cout << (jdim>0 ? ", " : "") << cutvals.at(idim).at(jdim);
      cout << endl;
    }

    // Make a canvas for the legend
    {
      unsigned int nleft=0;
      for (unsigned int ip=0; ip<proc_order.size(); ip++){
        TString const& procname = proc_order.at(ip);
        auto const& proccode = proc_code.at(ip);
        if (proccode!=0 || (is_3l1nu && procname=="ttbar_2l2nu")) nleft++;
      }
      unsigned int nright = proc_order.size()-nleft;

      TCanvas canvas(
        TString((isEnriched ? "c_SignalEnriched_" : "c_")) + TString((markPreliminary ? "Preliminary_" : "")) + (onORoffshell ? "Offshell_" : "Onshell_") + (aihypo=="" ? "SM" : aihypo.Data()) + "_legend",
        "", 8, 30, 800, 500
      );
      canvas.cd();
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(0.17);
      canvas.SetRightMargin(0.05);
      canvas.SetTopMargin(0.07);
      canvas.SetBottomMargin(0.13);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);

      float leg_xmin=0.55;
      float leg_ymin=0.10; if (nright<nleft) leg_ymin += (0.9-0.1)/std::max(nleft, nright)*(nleft-nright);
      float leg_xmax=0.90;
      float leg_ymax=0.90;
      TLegend legend_right(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
      legend_right.SetBorderSize(0);
      legend_right.SetTextFont(42);
      legend_right.SetTextSize(0.045);
      legend_right.SetLineColor(1);
      legend_right.SetLineStyle(1);
      legend_right.SetLineWidth(1);
      legend_right.SetFillColor(0);
      legend_right.SetFillStyle(0);

      leg_xmin=0.10;
      leg_xmax=0.45;
      leg_ymin=0.10; if (nright>nleft) leg_ymin += (0.9-0.1)/std::max(nleft, nright)*(nright-nleft);
      TLegend legend_left(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
      legend_left.SetBorderSize(0);
      legend_left.SetTextFont(42);
      legend_left.SetTextSize(0.045);
      legend_left.SetLineColor(1);
      legend_left.SetLineStyle(1);
      legend_left.SetLineWidth(1);
      legend_left.SetFillColor(0);
      legend_left.SetFillStyle(0);

      TGraphAsymmErrors* tgdata=nullptr;
      for (unsigned int ip=0; ip<proc_order.size(); ip++){
        TString const& procname = proc_order.at(ip);
        TString const& proclabel = proc_label.at(ip);

        TH1F*& prochist = procdist[procname].front();
        if (!prochist) cout << procname << " histogram is null!" << endl;
        else cout << procname << " histogram is present." << endl;
        if (procname.Contains("ALT")){
          if (procname.Contains("GGsm")) prochist->SetLineStyle(2);
          else if (procname.Contains("fai1")) prochist->SetLineStyle(7);
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColor(proc_color[ip]);
        }
        else if (proclabel=="Total SM"){
          prochist->SetMarkerColor(kRed);
          prochist->SetLineColor(kRed);
        }
        else if (proclabel=="gg#rightarrow4l SM s+b+i"){
          prochist->SetMarkerColor(kRed);
          prochist->SetLineColor(kRed);
          prochist->SetFillColor(kRed-7);
          prochist->SetFillStyle(3345);
        }
        else if (proclabel=="VBF+VH SM"){
          prochist->SetMarkerColor(kBlue);
          prochist->SetLineColor(kBlue);
        }
        else if (proclabel=="EW SM s+b+i"){
          prochist->SetMarkerColor(kBlue);
          prochist->SetLineColor(kBlue);
          prochist->SetFillColor(kBlue-7);
          //prochist->SetFillStyle(3354);
          prochist->SetFillStyle(3344);
        }
        else if (proc_code.at(ip)==0){
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColorAlpha(kBlack, 0.5);
          prochist->SetFillColor(proc_color[ip]);
          prochist->SetFillStyle(1001);
        }
        else if (proc_code.at(ip)==-99){ // Data
          prochist->SetMarkerColor(proc_color[ip]);
        }
        else{
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColorAlpha(kRed+1, 0.5);
          prochist->SetFillColor(proc_color[ip]);
          prochist->SetFillStyle(1001);
        }
        prochist->SetLineWidth(1);

        if (procname=="data") tgdata=getDataGraph(prochist);
      }

      for (unsigned int ip=proc_order.size(); ip>0; ip--){
        TString const& procname = proc_order.at(ip-1);
        TString const& proclabel = proc_label.at(ip-1);
        auto const& proccode = proc_code.at(ip-1);
        cout << "Adding process " << procname << " to legend..." << endl;
        TH1F*& prochist = procdist[procname].front();
        if (!prochist) cout << procname << " histogram is null!" << endl;

        TLegend* legend_chosen = nullptr;
        if (proccode!=0 || (is_3l1nu && procname=="ttbar_2l2nu")){
          legend_chosen = &legend_left;
          cout << "\t- Plotting on the left..." << endl;
        }
        else{
          legend_chosen = &legend_right;
          cout << "\t- Plotting on the right..." << endl;
        }

        if (proclabel=="Total SM" || procname.Contains("ALT")) legend_chosen->AddEntry(prochist, proc_label.at(ip-1), "l");
        else if (procname!="data") legend_chosen->AddEntry(prochist, proc_label.at(ip-1), "f");
        else if (tgdata) legend_chosen->AddEntry(tgdata, proc_label.at(ip-1), "e1p");
      }

      legend_left.Draw();
      legend_right.Draw();
      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(TString(canvas.GetName())+".pdf");

      delete tgdata;
    }


    for (unsigned int idim=0; idim<ndims; idim++){
      TString dimname = Form("_KD%i", idim+1);
      cout << "Ploting dimension " << idim << "..." << endl;

      for (unsigned int ip=0; ip<proc_order.size(); ip++){
        if (ip>0){
          if (!proc_order.at(ip).Contains("ALT") && proc_order.at(ip)!="data") procdist[proc_order.at(ip)].at(idim)->Add(procdist[proc_order.at(ip-1)].at(idim));
          else if (proc_order.at(ip)!="data"){
            for (unsigned int jp=ip-1; jp>0; jp--){
              if (proc_code[jp]==0){
                procdist[proc_order.at(ip)].at(idim)->Add(procdist[proc_order.at(jp)].at(idim)); // Add the first bkg process and break
                break;
              }
            }
          }
        }
        cout << proc_order.at(ip) << " integral: " << getHistogramIntegralAndError(procdist[proc_order.at(ip)].at(idim), 1, procdist[proc_order.at(ip)].at(idim)->GetNbinsX(), false, nullptr) << endl;
      }

      // Divide by bin width when plotting off-shell mass
      //if (onORoffshell==1 && (int) idim==massDim){ for (auto proc:proc_order) divideBinWidth(procdist[proc]); }

      // Draw
      cout << "Creating the canvas..." << endl;
      TCanvas canvas(
        TString((isEnriched ? "c_SignalEnriched_" : "c_")) + TString((markPreliminary ? "Preliminary_" : "")) + (onORoffshell ? "Offshell_" : "Onshell_") + catname + "_" + (aihypo=="" ? "SM" : aihypo.Data()) + dimname + (!isPostFit ? "" : "_postfit"),
        "", 8, 30, 800, 1100
      );
      canvas.cd();
      /*
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(0.17);
      canvas.SetRightMargin(0.05);
      canvas.SetTopMargin(0.07);
      canvas.SetBottomMargin(0.13);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      */

      cout << "\t- Creating the pads..." << endl;
      std::vector<TPad*> pads;
      canvas.cd();
      pads.push_back(
        new TPad(
          Form("%s_top", canvas.GetName()), "",
          0, 2./11. + 0.13 * 8./11. + 0.5/11., 1, 1
        )
      );
      canvas.cd();
      pads.push_back(
        new TPad(
          Form("%s_bot", canvas.GetName()), "",
          0, 0, 1, 2./11. + 0.13 * 8./11. + 0.5/11.
        )
      );
      {
        unsigned int ipad=0;
        for (auto& pad:pads){
          pad->cd();
          pad->SetFillColor(0);
          pad->SetBorderMode(0);
          pad->SetBorderSize(2);
          pad->SetTickx(1);
          pad->SetTicky(1);
          pad->SetLeftMargin(0.17);
          pad->SetRightMargin(0.05);
          if (ipad==0){
            pad->SetTopMargin(0.07*800./(1100.*(1.-(2./11. + 0.13 * 8./11. + 0.5/11.))));
            pad->SetBottomMargin(0.5/11.);
          }
          else{
            pad->SetTopMargin(0.5/11.);
            pad->SetBottomMargin(0.13*800./(1100.*(2./11. + 0.13 * 8./11. + 0.5/11.)));
          }
          pad->SetFrameFillStyle(0);
          pad->SetFrameBorderMode(0);
          pad->SetFrameFillStyle(0);
          pad->SetFrameBorderMode(0);
          if (varlabels.at(idim).Contains("m_{T}^{ZZ}") || varlabels.at(idim).Contains("m_{T}^{WZ}") || varlabels.at(idim).Contains("p_{T}^{miss}")) pad->SetLogx(true);
          canvas.cd();
          ipad++;
        }
      }
      // Draw in reverse order
      pads.back()->Draw();
      pads.front()->Draw();

      TText* text;
      TPaveText pt(0.15, 0.93, 0.85, 1, "brNDC");
      pt.SetBorderSize(0);
      pt.SetFillStyle(0);
      pt.SetTextAlign(12);
      pt.SetTextFont(42);
      pt.SetTextSize(0.045);
      text = pt.AddText(0.025, 0.45, "#font[61]{CMS}");
      text->SetTextSize(0.044);
      if (markPreliminary){
        text = pt.AddText(0.165, 0.42, "#font[52]{Preliminary}");
        text->SetTextSize(0.0315);
      }
      else if (onORoffshell==1 && !(isEnriched && ((cinputdir.Contains("/SM") && idim!=1) || (cinputdir.Contains("/a3") && idim==1)))){
        text = pt.AddText(0.165, 0.42, "#font[52]{Supplementary}");
        text->SetTextSize(0.0315);
      }
      int theSqrts=13;
      TString cErgTev = Form("#font[42]{137.2 fb^{-1} (%i TeV)}", theSqrts);
      text = pt.AddText(0.78, 0.45, cErgTev);
      text->SetTextSize(0.0315);

      TString strCatLabel;
      if (catname=="Untagged") strCatLabel="Untagged";
      else if (catname=="JJVBFTagged") strCatLabel="VBF-tagged";
      else if (catname=="HadVHTagged") strCatLabel="VH-tagged";
      else if (catname=="Nj_eq_0") strCatLabel="N_{j}=0";
      else if (catname=="Nj_eq_1") strCatLabel="N_{j}=1";
      else if (catname.Contains("Nj_geq_2")) strCatLabel="N_{j}#geq2";
      else if (catname=="BoostedHadVH") strCatLabel="Boosted V #rightarrow J";
      if (ailabel!="") strCatLabel = strCatLabel + ", " + ailabel + " analysis";
      if (isPostFit) strCatLabel = strCatLabel + " (postfit)";
      TPaveText pt_cat(0.20, 0.83, 0.40, 0.90, "brNDC");
      pt_cat.SetBorderSize(0);
      pt_cat.SetFillStyle(0);
      pt_cat.SetTextAlign(12);
      pt_cat.SetTextFont(42);
      pt_cat.SetTextSize(0.03);
      text = pt_cat.AddText(0.02, 0.45, strCatLabel);
      text->SetTextSize(0.044);

      cout << "\t- Preparing the cut label..." << endl;
      TString strCutLabel;
      if ((int) idim!=massDim && massDim>=0 && cutvals.at(idim).at(massDim)>0.){
        TString strUnit = "GeV";
        if (is_2l2nu) strCutLabel += "m_{T}^{ZZ}#geq";
        else if (is_3l1nu) strCutLabel += "m_{T}^{WZ}#geq";
        else strCutLabel += "m_{4l}#geq";
        strCutLabel += Form("%.0f", cutvals.at(idim).at(massDim));
        if (strUnit!="") strCutLabel = strCutLabel + " " + strUnit;
      }
      if (is_2l2nu && catname.Contains("Nj_geq_2")){
        if (strCutLabel!="") strCutLabel += ", ";
        if (catname.Contains("pTmiss_lt_200")) strCutLabel += "p_{T}^{miss}<200 GeV";
        else if (catname.Contains("pTmiss_ge_200")) strCutLabel += "p_{T}^{miss}#geq200 GeV";
      }
      if ((int) idim!=kdDim && kdDim>=0 && cutvals.at(idim).at(kdDim)>0.){
        TString strUnit = "";
        if (strCutLabel!="") strCutLabel += ", ";
        if (is_2l2nu){
          strCutLabel += (catname.Contains("Nj_geq_2") ? "D_{2jet}^{VBF}#geq" : "p_{T}^{miss}#geq");
          if (!catname.Contains("Nj_geq_2")) strUnit = "GeV";
        }
        else strCutLabel += "D_{bkg}#geq";
        if (cutvals.at(idim).at(kdDim)>1.) strCutLabel += Form("%.0f", cutvals.at(idim).at(kdDim));
        else strCutLabel += Form("%.1f", cutvals.at(idim).at(kdDim));
        if (strUnit!="") strCutLabel = strCutLabel + " " + strUnit;
      }
      TPaveText pt_cut(0.20, 0.75, 0.40, 0.83, "brNDC");
      pt_cut.SetBorderSize(0);
      pt_cut.SetFillStyle(0);
      pt_cut.SetTextAlign(12);
      pt_cut.SetTextFont(42);
      pt_cut.SetTextSize(0.03);
      text = pt_cut.AddText(0.02, 0.45, strCutLabel);
      text->SetTextSize(0.044);

      float ymax=-1, ymin=-1;
      float xmin=-1, xmax=-1;
      TGraphAsymmErrors* tgdata=nullptr;
      for (unsigned int ip=0; ip<proc_order.size(); ip++){
        TString const& procname = proc_order.at(ip);
        TString const& proclabel = proc_label.at(ip);
        cout << "Adjusting process " << procname << " at index " << ip << endl;
        TH1F*& prochist = procdist[procname].at(idim);
        if (!prochist) cout << procname << " histogram is null!" << endl;
        else cout << procname << " histogram is present." << endl;
        if (onORoffshell && idim==0 && is_4l){ xmin=220; xmax=1000; }
        else{
          xmin=prochist->GetXaxis()->GetBinLowEdge(1);
          xmax=prochist->GetXaxis()->GetBinUpEdge(prochist->GetNbinsX());
        }
        cout << "\t- Process " << procname << " color: " << proc_color[ip] << endl;
        if (!prochist) cout << "ERROR: PDF for " << procname << " missing!" << endl;
        if (procname.Contains("ALT")){
          if (procname.Contains("GGsm")) prochist->SetLineStyle(2);
          else if (procname.Contains("fai1")) prochist->SetLineStyle(7);
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColor(proc_color[ip]);
        }
        else if (proclabel=="Total SM"){
          prochist->SetMarkerColor(kRed);
          prochist->SetLineColor(kRed);
        }
        else if (proclabel=="gg#rightarrow4l SM s+b+i"){
          prochist->SetMarkerColor(kRed);
          prochist->SetLineColor(kRed);
          prochist->SetFillColor(kRed-7);
          prochist->SetFillStyle(3345);
        }
        else if (proclabel=="VBF+VH SM"){
          prochist->SetMarkerColor(kBlue);
          prochist->SetLineColor(kBlue);
        }
        else if (proclabel=="EW SM s+b+i"){
          prochist->SetMarkerColor(kBlue);
          prochist->SetLineColor(kBlue);
          prochist->SetFillColor(kBlue-7);
          //prochist->SetFillStyle(3354);
          prochist->SetFillStyle(3344);
        }
        else if (proc_code.at(ip)==0){
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColorAlpha(kBlack, 0.5);
          prochist->SetFillColor(proc_color[ip]);
          prochist->SetFillStyle(1001);
        }
        else if (proc_code.at(ip)==-99){ // Data
          prochist->SetMarkerColor(proc_color[ip]);
        }
        else{
          prochist->SetMarkerColor(proc_color[ip]);
          //prochist->SetLineColorAlpha(proc_color[ip], 0.9);
          //prochist->SetFillColor(0);
          prochist->SetLineColorAlpha(kRed+1, 0.5);
          prochist->SetFillColor(proc_color[ip]);
          prochist->SetFillStyle(1001);
        }
        prochist->SetLineWidth(1);

        cout << "\t- Adding overflow content" << endl;
        int binXlow = prochist->GetXaxis()->FindBin(xmin);
        int binXhigh = prochist->GetXaxis()->FindBin(xmax);
        if (prochist->GetXaxis()->GetBinLowEdge(binXhigh)==xmax) binXhigh--;
        for (int ix=binXhigh+1; ix<=prochist->GetNbinsX(); ix++){
          prochist->SetBinContent(binXhigh, prochist->GetBinContent(ix)+prochist->GetBinContent(binXhigh));
          prochist->SetBinError(binXhigh, sqrt(pow(prochist->GetBinError(ix), 2)+pow(prochist->GetBinContent(binXhigh), 2)));
        }
        for (int ix=1; ix<binXlow; ix++){
          prochist->SetBinContent(binXlow, prochist->GetBinContent(ix)+prochist->GetBinContent(binXlow));
          prochist->SetBinError(binXlow, sqrt(pow(prochist->GetBinError(ix), 2)+pow(prochist->GetBinContent(binXlow), 2)));
        }
        prochist->GetXaxis()->SetRangeUser(xmin, xmax);
        prochist->GetXaxis()->SetNdivisions(505);
        prochist->GetXaxis()->SetLabelFont(42);
        prochist->GetXaxis()->SetLabelOffset(0.007);
        prochist->GetXaxis()->SetLabelSize(0.04);
        prochist->GetXaxis()->SetTitleSize(0.06);
        prochist->GetXaxis()->SetTitleOffset(0.9);
        prochist->GetXaxis()->SetTitleFont(42);
        prochist->GetYaxis()->SetNdivisions(505);
        prochist->GetYaxis()->SetLabelFont(42);
        prochist->GetYaxis()->SetLabelOffset(0.007);
        prochist->GetYaxis()->SetLabelSize(0.04);
        prochist->GetYaxis()->SetTitleSize(0.06);
        prochist->GetYaxis()->SetTitleOffset(1.1);
        prochist->GetYaxis()->SetTitleFont(42);
        prochist->GetYaxis()->SetTitle("Events / bin");
        prochist->GetYaxis()->CenterTitle();
        prochist->GetXaxis()->CenterTitle();

        if (procname!="data"){
          for (int ix=binXlow; ix<=binXhigh; ix++){
            float bc = prochist->GetBinContent(ix);
            if (bc!=0.){
              ymax = std::max(bc, ymax);
            }
          }
        }
        else{
          cout << "\t- Obtaining " << procname << " graph" << endl;
          tgdata=getDataGraph(prochist);
          if (tgdata){
            cout << "\t\t- Np = " << tgdata->GetN() << endl;
            for (int ipoint=0; ipoint<tgdata->GetN(); ipoint++){
              float bc = tgdata->GetY()[ipoint]+tgdata->GetEYhigh()[ipoint];
              if (bc!=0.){
                ymax = std::max(bc, ymax);
              }
            }
            cout << "\t\t- Success!" << endl;
          }
          else cout << "-t-t- Failure!" << endl;
        }
      }
      ymin=0;

      float ymaxfactor = 1.25;
      float yminfactor = 1;

      vector<TH1F*> intermediateHistList;
      canvas.cd();
      pads.front()->cd();
      bool drawfirst=true;
      for (unsigned int ip=proc_order.size(); ip>0; ip--){
        TString const& procname = proc_order.at(ip-1);
        TString const& proclabel = proc_label.at(ip-1);
        TH1F*& prochist = procdist[procname].at(idim);
        if (procname=="data") continue;
        if (!proc_order.at(ip-1).Contains("ALT")) intermediateHistList.push_back((TH1F*) prochist->Clone(Form("%s_copy", prochist->GetName())));
        prochist->GetXaxis()->SetTitle("");
        prochist->GetXaxis()->SetLabelSize(0);
        prochist->GetXaxis()->SetTitleSize(0);
        cout << "\t- Drawing " << procname << endl;
        prochist->GetYaxis()->SetRangeUser(ymin*yminfactor*0.8, ymax*ymaxfactor);
        prochist->Draw((drawfirst ? "hist" : "histsame"));
        drawfirst=false;
      }
      // Re-draw ALT
      for (unsigned int ip=proc_order.size(); ip>0; ip--){
        if (!proc_order.at(ip-1).Contains("ALT")) continue;
        cout << "\t- Drawing " << proc_order.at(ip-1) << endl;
        procdist[proc_order.at(ip-1)].at(idim)->Draw("histsame");
      }
      if (tgdata && tgdata->GetN()>0 && !isBlind) tgdata->Draw("e1psame");
      pt.Draw();
      pt_cat.Draw();
      if (strCutLabel!="") pt_cut.Draw();

      canvas.cd();
      pads.back()->cd();
      drawfirst=false;
      for (int ix=1; ix<=intermediateHistList.front()->GetNbinsX(); ix++){
        double bc_sum = intermediateHistList.front()->GetBinContent(ix);
        for (auto& hh:intermediateHistList) hh->SetBinContent(ix, hh->GetBinContent(ix)/bc_sum);
      }
      for (auto& hh:intermediateHistList){
        hh->GetYaxis()->SetRangeUser(0, 1);
        hh->GetXaxis()->SetLabelSize(0.04*800./(1100.*(2./11. + 0.13 * 8./11. + 0.5/11.)));
        hh->GetXaxis()->SetTitleSize(0.06*800./(1100.*(2./11. + 0.13 * 8./11. + 0.5/11.)));
        hh->GetYaxis()->SetLabelSize(0.04*800./(1100.*(2./11. + 0.13 * 8./11. + 0.5/11.)));
        hh->GetYaxis()->SetTitleSize(0.06*800./(1100.*(2./11. + 0.13 * 8./11. + 0.5/11.)));
        hh->GetYaxis()->SetTitle("Ratio");
        hh->GetYaxis()->SetTitleOffset(1.1*(2./11. + 0.13 * 8./11. + 0.5/11.) / (1.-(2./11. + 0.13 * 8./11. + 0.5/11.)));
        hh->Draw((drawfirst ? "hist" : "histsame"));
        drawfirst=false;
      }

      for (auto& pad:pads){
        pad->RedrawAxis();
        pad->Modified();
        pad->Update();
      }

      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(TString(canvas.GetName()) + (isBlind ? "_blind" : "") + ".pdf");
      for (TH1F*& htmp:intermediateHistList) delete htmp;
      for (auto& pad:pads) pad->Close();
      canvas.Close();
      curdir->cd();

      delete tgdata;
      for (auto it=procdist.begin(); it!=procdist.end(); it++){
        cout << "\t- Deleting histogram " << it->second.at(idim)->GetName() << endl;
        delete it->second.at(idim);
      }
    }

  }
}


template<> void divideBinWidth<TH1F>(TH1F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    histo->SetBinContent(binx, histo->GetBinContent(binx)/binwidthX);
    histo->SetBinError(binx, histo->GetBinError(binx)/binwidthX);
  }
}
template<> void divideBinWidth<TH2F>(TH2F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      float binwidth=binwidthX*binwidthY;
      histo->SetBinContent(binx, biny, histo->GetBinContent(binx, biny)/binwidth);
      histo->SetBinError(binx, biny, histo->GetBinError(binx, biny)/binwidth);
    }
  }
}
template<> void divideBinWidth<TH3F>(TH3F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  TAxis const* zaxis = histo->GetZaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      for (int binz=1; binz<=histo->GetNbinsZ(); binz++){
        float binwidthZ = zaxis->GetBinWidth(binz);
        float binwidth=binwidthX*binwidthY*binwidthZ;
        histo->SetBinContent(binx, biny, binz, histo->GetBinContent(binx, biny, binz)/binwidth);
        histo->SetBinError(binx, biny, binz, histo->GetBinError(binx, biny, binz)/binwidth);
      }
    }
  }
}

template<> void multiplyBinWidth<TH1F>(TH1F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    histo->SetBinContent(binx, histo->GetBinContent(binx)*binwidthX);
    histo->SetBinError(binx, histo->GetBinError(binx)*binwidthX);
  }
}
template<> void multiplyBinWidth<TH2F>(TH2F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      float binwidth=binwidthX*binwidthY;
      histo->SetBinContent(binx, biny, histo->GetBinContent(binx, biny)*binwidth);
      histo->SetBinError(binx, biny, histo->GetBinError(binx, biny)*binwidth);
    }
  }
}
template<> void multiplyBinWidth<TH3F>(TH3F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  TAxis const* zaxis = histo->GetZaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      for (int binz=1; binz<=histo->GetNbinsZ(); binz++){
        float binwidthZ = zaxis->GetBinWidth(binz);
        float binwidth=binwidthX*binwidthY*binwidthZ;
        histo->SetBinContent(binx, biny, binz, histo->GetBinContent(binx, biny, binz)*binwidth);
        histo->SetBinError(binx, biny, binz, histo->GetBinError(binx, biny, binz)*binwidth);
      }
    }
  }
}

template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };

      res=histo->IntegralAndError(xb[0], xb[1], reserror, "width");

      double integralinside, integralerrorinside;
      integralinside=histo->IntegralAndError(xb[0], xb[1], integralerrorinside, "");

      double integraloutside, integralerroroutside;
      integraloutside=histo->IntegralAndError(ix, jx, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, iy, jy, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };
      int yb[2]={ std::max(1, std::min(histo->GetNbinsY(), iy)), std::max(1, std::min(histo->GetNbinsY(), jy)) };

      res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], reserror, "width");

      double integralinside, integralerrorinside;
      integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], integralerrorinside, "");

      double integraloutside, integralerroroutside;
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };
      int yb[2]={ std::max(1, std::min(histo->GetNbinsY(), iy)), std::max(1, std::min(histo->GetNbinsY(), jy)) };
      int zb[2]={ std::max(1, std::min(histo->GetNbinsZ(), iz)), std::max(1, std::min(histo->GetNbinsZ(), jz)) };

      res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], reserror, "width");

      double integralinside, integralerrorinside;
      integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], integralerrorinside, "");

      double integraloutside, integralerroroutside;
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}

TH1F* getHistogramSlice(TH1F const* histo, TString newname){
  if (!histo) return nullptr;
  if (newname=="") newname=Form("Slice_%s", histo->GetName());

  const TAxis* xaxis=histo->GetXaxis();
  vector<float> bins;
  for (int i=1; i<=xaxis->GetNbins()+1; i++) bins.push_back(xaxis->GetBinLowEdge(i));
  TH1F* res = new TH1F(newname, "", bins.size()-1, bins.data());

  for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
    double integral=0, integralerror=0;
    integral = getHistogramIntegralAndError(histo, ii, ii, false, &integralerror);
    res->SetBinContent(ii, integral);
    res->SetBinError(ii, integralerror);
  }

  return res;
}
TH1F* getHistogramSlice(TH2F const* histo, unsigned char XDirection, int iy, int jy, TString newname){
  if (!histo || XDirection>=2) return nullptr;
  if (newname=="") newname=Form("Slice_%s_%i_%i_%s", (XDirection==0 ? "X" : "Y"), iy, jy, histo->GetName());

  const TAxis* xaxis=histo->GetXaxis();
  const TAxis* yaxis=histo->GetYaxis();
  vector<float> bins;
  if (XDirection==0){
    for (int i=1; i<=xaxis->GetNbins()+1; i++) bins.push_back(xaxis->GetBinLowEdge(i));
    iy = std::max(0, iy); jy = std::min(yaxis->GetNbins()+1, jy);
  }
  else{
    for (int i=1; i<=yaxis->GetNbins()+1; i++) bins.push_back(yaxis->GetBinLowEdge(i));
    iy = std::max(0, iy); jy = std::min(xaxis->GetNbins()+1, jy);
  }
  if (iy>jy) cerr << "getHistogramSlice: iy>jy!" << endl;
  TH1F* res = new TH1F(newname, "", bins.size()-1, bins.data());

  if (XDirection==0){
    for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
      double integral=0, integralerror=0;
      integral = getHistogramIntegralAndError(histo, ii, ii, iy, jy, false, &integralerror);
      res->SetBinContent(ii, integral);
      res->SetBinError(ii, integralerror);
    }
  }
  else{
    for (int ii=0; ii<=yaxis->GetNbins()+1; ii++){
      double integral=0, integralerror=0;
      integral = getHistogramIntegralAndError(histo, iy, jy, ii, ii, false, &integralerror);
      res->SetBinContent(ii, integral);
      res->SetBinError(ii, integralerror);
    }
  }

  return res;
}
TH1F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, int iy, int jy, int iz, int jz, TString newname){
  if (!histo || XDirection>=3) return nullptr;
  if (newname=="") newname=Form("Slice_%s_%i_%i_%i_%i_%s", (XDirection==0 ? "X" : (XDirection==1 ? "Y" : "Z")), iy, jy, iz, jz, histo->GetName());

  const TAxis* xaxis;
  const TAxis* yaxis;
  const TAxis* zaxis;
  vector<float> bins;
  if (XDirection==0){
    xaxis=histo->GetXaxis();
    yaxis=histo->GetYaxis();
    zaxis=histo->GetZaxis();
  }
  else if (XDirection==1){
    xaxis=histo->GetYaxis();
    yaxis=histo->GetZaxis();
    zaxis=histo->GetXaxis();
  }
  else{
    xaxis=histo->GetZaxis();
    yaxis=histo->GetXaxis();
    zaxis=histo->GetYaxis();
  }
  for (int i=1; i<=xaxis->GetNbins()+1; i++) bins.push_back(xaxis->GetBinLowEdge(i));
  iy = std::max(0, iy); jy = std::min(yaxis->GetNbins()+1, jy);
  iz = std::max(0, iz); jz = std::min(zaxis->GetNbins()+1, jz);
  if (iy>jy) cerr << "getHistogramSlice: iy>jy!" << endl;
  if (iz>jz) cerr << "getHistogramSlice: iz>jz!" << endl;
  TH1F* res = new TH1F(newname, "", bins.size()-1, bins.data());

  for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
    double integral=0, integralerror=0;
    int IX, JX, IY, JY, IZ, JZ;
    if (XDirection==0){
      IX=ii; JX=ii;
      IY=iy; JY=jy;
      IZ=iz; JZ=jz;
    }
    else if (XDirection==1){
      IX=iz; JX=jz;
      IY=ii; JY=ii;
      IZ=iy; JZ=jy;
    }
    else{
      IX=iy; JX=jy;
      IY=iz; JY=jz;
      IZ=ii; JZ=ii;
    }
    integral = getHistogramIntegralAndError(histo, IX, JX, IY, JY, IZ, JZ, false, &integralerror);
    res->SetBinContent(ii, integral);
    res->SetBinError(ii, integralerror);
  }

  return res;
}
TH2F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, unsigned char YDirection, int iz, int jz, TString newname){
  if (!histo || XDirection==YDirection || XDirection>=3 || YDirection>=3) return nullptr;
  if (newname=="") newname=Form("Slice_%s%s_%i_%i_%s", (XDirection==0 ? "X" : (XDirection==1 ? "Y" : "Z")), (YDirection==0 ? "X" : (YDirection==1 ? "Y" : "Z")), iz, jz, histo->GetName());

  unsigned char ZDirection=3-XDirection-YDirection; // 0+1+2=3
  const TAxis* xaxis;
  const TAxis* yaxis;
  const TAxis* zaxis;
  vector<float> xbins, ybins;
  if (XDirection==0) xaxis=histo->GetXaxis();
  else if (XDirection==1) xaxis=histo->GetYaxis();
  else xaxis=histo->GetZaxis();
  if (YDirection==0) yaxis=histo->GetXaxis();
  else if (YDirection==1) yaxis=histo->GetYaxis();
  else yaxis=histo->GetZaxis();
  if (ZDirection==0) zaxis=histo->GetXaxis();
  else if (ZDirection==1) zaxis=histo->GetYaxis();
  else zaxis=histo->GetZaxis();

  for (int i=1; i<=xaxis->GetNbins()+1; i++) xbins.push_back(xaxis->GetBinLowEdge(i));
  for (int i=1; i<=yaxis->GetNbins()+1; i++) ybins.push_back(yaxis->GetBinLowEdge(i));
  iz = std::max(0, iz); std::min(zaxis->GetNbins()+1, jz);
  TH2F* res = new TH2F(newname, "", xbins.size()-1, xbins.data(), ybins.size()-1, ybins.data());

  for (int ii=0; ii<=xaxis->GetNbins()+1; ii++){
    for (int jj=0; jj<=yaxis->GetNbins()+1; jj++){
      double integral=0, integralerror=0;
      int IX=0, JX=0, IY=0, JY=0, IZ=0, JZ=0;
      if (XDirection==0){ IX=ii; JX=ii; }
      else if (XDirection==1){ IY=ii; JY=ii; }
      else{ IZ=ii; JZ=ii; }
      if (YDirection==0){ IX=jj; JX=jj; }
      else if (YDirection==1){ IY=jj; JY=jj; }
      else{ IZ=jj; JZ=jj; }
      if (ZDirection==0){ IX=iz; JX=jz; }
      else if (ZDirection==1){ IY=iz; JY=jz; }
      else{ IZ=iz; JZ=jz; }
      integral = getHistogramIntegralAndError(histo, IX, JX, IY, JY, IZ, JZ, false, &integralerror);
      res->SetBinContent(ii, jj, integral);
      res->SetBinError(ii, jj, integralerror);
    }
  }

  return res;
}

void splitOption(const string rawoption, string& wish, string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
    while (value.find(delimiter)==0) value=value.substr(1);
  }
  else{
    wish="";
    value=rawoption;
  }
}
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="") splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}

void splitOption(const TString rawoption, TString& wish, TString& value, char delimiter){
  string srawoption = rawoption.Data();
  string swish, svalue;
  splitOption(srawoption, swish, svalue, delimiter);
  wish = swish.data();
  value = svalue.data();
}
void splitOptionRecursive(const TString rawoption, vector<TString>& splitoptions, char delimiter){
  string srawoption = rawoption.Data();
  vector<string> soptions;
  splitOptionRecursive(srawoption, soptions, delimiter);
  splitoptions.clear(); splitoptions.reserve(soptions.size());
  for (auto const& sopt:soptions) splitoptions.push_back(sopt.data());
}
