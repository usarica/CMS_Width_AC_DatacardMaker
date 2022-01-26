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
#include <HiggsAnalysis/CombinedLimit/interface/AsymPow.h>
#include <IvyFramework/IvyDataTools/interface/HelperFunctions.h>
#include <IvyFramework/IvyDataTools/interface/HostHelpersCore.h>


using namespace RooFit;
using namespace std;

struct process_spec{
  RooAbsPdf* pdf;
  RooAbsReal* norm;
  double rate;
  TString name;
  std::vector<RooAbsReal*> rateModifiers;

  std::unordered_map<TString, TH1F*> syst_h1D_map;
  std::unordered_map<TString, TH2F*> syst_h2D_map;
  std::unordered_map<TString, TH3F*> syst_h3D_map;

  process_spec() : pdf(0), norm(0), rate(0){}
  process_spec(TString name_, double rate_) : pdf(nullptr), norm(nullptr), rate(rate_), name(name_){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()){}
  process_spec(const process_spec& other) : pdf(other.pdf), norm(other.norm), rate(other.rate), name(other.name){}
  ~process_spec();

  void setRateModifiers(std::vector<RooAbsReal*> const& vlist){ rateModifiers=vlist; }
  double getExtraNorm();

  bool dependsOn(RooAbsReal* var);

  void setDistribution(TH1F* tpl, TString const& procname);
  void setDistribution(TH2F* tpl, TString const& procname);
  void setDistribution(TH3F* tpl, TString const& procname);
};
process_spec::~process_spec(){
  for (auto& pp:syst_h1D_map) delete pp.second;
  for (auto& pp:syst_h2D_map) delete pp.second;
  for (auto& pp:syst_h3D_map) delete pp.second;
}

double process_spec::getExtraNorm(){
  double res = 1;
  for (auto const& var:rateModifiers) res *= var->getVal();
  return res;
}

bool process_spec::dependsOn(RooAbsReal* var){
  bool res = false;
  if (!res && pdf) res |= pdf->dependsOn(*var);
  if (!res && norm) res |= norm->dependsOn(*var);
  for (auto const& rateModifier:rateModifiers){
    if (!res) res |= rateModifier->dependsOn(*var);
  }
  return res;
}

void process_spec::setDistribution(TH1F* tpl, TString const& procname){
  if (syst_h1D_map.find(procname)!=syst_h1D_map.end()){
    cerr << "syst_h1D_map[" << name << "] already contains " << procname << "." << endl;
    exit(1);
  }
  TH1F* tplnew = new TH1F(*tpl); tplnew->SetName(procname);
  syst_h1D_map[procname] = tplnew;
}
void process_spec::setDistribution(TH2F* tpl, TString const& procname){
  if (syst_h2D_map.find(procname)!=syst_h2D_map.end()){
    cerr << "syst_h2D_map[" << name << "] already contains " << procname << "." << endl;
    exit(1);
  }
  TH2F* tplnew = new TH2F(*tpl); tplnew->SetName(procname);
  syst_h2D_map[procname] = tplnew;
}
void process_spec::setDistribution(TH3F* tpl, TString const& procname){
  if (syst_h3D_map.find(procname)!=syst_h3D_map.end()){
    cerr << "syst_h3D_map[" << name << "] already contains " << procname << "." << endl;
    exit(1);
  }
  TH3F* tplnew = new TH3F(*tpl); tplnew->SetName(procname);
  syst_h3D_map[procname] = tplnew;
}



using namespace HelperFunctions;



unsigned int extractTemplates(
  process_spec& proc, RooDataSet* data,
  unordered_map<TString, TH1F>& procshape_1D, unordered_map<TString, TH2F>& procshape_2D, unordered_map<TString, TH3F>& procshape_3D,
  TString newname="",
  bool forceAddShapeToProcSpec=false
){
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

  bool addShapeToProcSpec = true;
  TString procname=proc.name;
  if (newname!=""){
    procname=newname;
    addShapeToProcSpec = forceAddShapeToProcSpec;
  }
  TString tplname = "T_";
  tplname += procname;

  if (ndims==3){
    TH3F* tpl=(TH3F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    multiplyBinWidth(tpl);
    double normval = proc.rate * proc.getExtraNorm(); if (proc.norm) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    //cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;

    tpl->SetTitle("");
    tpl->Scale(scale);
    //cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), 1, tpl->GetNbinsZ(), false, nullptr) << "?=" << normval << endl;

    if (addShapeToProcSpec) proc.setDistribution(tpl, procname);

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
    //cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_3D[procname]), 1, procshape_3D[procname].GetNbinsX(), 1, procshape_3D[procname].GetNbinsY(), 1, procshape_3D[procname].GetNbinsZ(), false, nullptr) << "?=" << normval << endl;
    delete tpl;
  }
  else if (ndims==2){
    TH2F* tpl=(TH2F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    multiplyBinWidth(tpl);
    double normval = proc.rate * proc.getExtraNorm(); if (proc.norm) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    //cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;

    tpl->SetTitle("");
    tpl->Scale(scale);
    //cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), false, nullptr) << "?=" << normval << endl;

    if (addShapeToProcSpec) proc.setDistribution(tpl, procname);

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
    //cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_2D[procname]), 1, procshape_2D[procname].GetNbinsX(), 1, procshape_2D[procname].GetNbinsY(), false, nullptr) << "?=" << normval << endl;
    delete tpl;
  }
  else if (ndims==1){
    TH1F* tpl=(TH1F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd, zcmd);
    multiplyBinWidth(tpl);
    double normval = proc.rate * proc.getExtraNorm(); if (proc.norm) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    //cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;

    tpl->SetTitle("");
    tpl->Scale(scale);
    //cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), false, nullptr) << "?=" << normval << endl;

    if (addShapeToProcSpec) proc.setDistribution(tpl, procname);

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
    //cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_1D[procname]), 1, procshape_1D[procname].GetNbinsX(), false, nullptr) << "?=" << normval << endl;
    delete tpl;
  }
  //cout << "Work on template " << tplname << " is complete!" << endl;
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
    //cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), 1, tpl->GetNbinsZ(), false, nullptr) << endl;

    bool isAdded=false;
    for (auto it=procshape_3D.begin(); it!=procshape_3D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      //cout << "\t- Creating " << tplname << endl;
      procshape_3D[procname]=TH3F(*tpl);
      procshape_3D[procname].SetName(tplname);
    }
    //cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_3D[procname]), 1, procshape_3D[procname].GetNbinsX(), 1, procshape_3D[procname].GetNbinsY(), 1, procshape_3D[procname].GetNbinsZ(), false, nullptr) << endl;
    delete tpl;
  }
  else if (ndims==2){
    TH2F* tpl=(TH2F*) data->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd);
    //cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), false, nullptr) << endl;

    bool isAdded=false;
    for (auto it=procshape_2D.begin(); it!=procshape_2D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      //cout << "\t- Creating " << tplname << endl;
      procshape_2D[procname]=TH2F(*tpl);
      procshape_2D[procname].SetName(tplname);
    }
    //cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_2D[procname]), 1, procshape_2D[procname].GetNbinsX(), 1, procshape_2D[procname].GetNbinsY(), false, nullptr) << endl;
    delete tpl;
  }
  else if (ndims==1){
    TH1F* tpl=(TH1F*) data->createHistogram(tplname+"_Copy", *(deps.at(0)), xcmd, ycmd);
    //cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), false, nullptr) << endl;

    bool isAdded=false;
    for (auto it=procshape_1D.begin(); it!=procshape_1D.end(); it++){
      if (it->first==procname){
        it->second.Add(tpl, 1.);
        isAdded=true;
        break;
      }
    }
    if (!isAdded){
      //cout << "\t- Creating " << tplname << endl;
      procshape_1D[procname]=TH1F(*tpl);
      procshape_1D[procname].SetName(tplname);
    }
    //cout << procname << " final integral = " << getHistogramIntegralAndError(&(procshape_1D[procname]), 1, procshape_1D[procname].GetNbinsX(), false, nullptr) << endl;
    delete tpl;
  }
  //cout << "Work on data " << tplname << " is complete!" << endl;
  return ndims;
}


bool isTrueNuisance(TString const& nuisname){
  return !(
    nuisname.BeginsWith("RF") || nuisname.BeginsWith("RV") || nuisname.BeginsWith("R_") || nuisname=="R"
    ||
    nuisname.BeginsWith("rf") || nuisname.BeginsWith("rv") || nuisname.BeginsWith("r_") || nuisname=="r"
    ||
    nuisname=="CMS_zz4l_fai1" || nuisname=="GGsm"
    );
}

std::unordered_map<TString, std::pair<RooRealVar*, RooRealVar*>> getZippedNuisanceMap(
  std::vector<RooRealVar*> const& nuisanceVars,
  std::vector<RooRealVar*> const& lnNmodvars
){
  std::unordered_map<TString, std::pair<RooRealVar*, RooRealVar*>> res;
  for (auto const& nuis:nuisanceVars){
    TString nuisname = nuis->GetName();
    res[nuisname] = std::pair<RooRealVar*, RooRealVar*>(nuis, nullptr);
  }
  for (auto const& nuis:lnNmodvars){
    TString nuisname = nuis->GetName();
    auto it = res.find(nuisname);
    if (it!=res.end()) it->second.second = nuis;
    else res[nuisname] = std::pair<RooRealVar*, RooRealVar*>(nullptr, nuis);
  }
  return res;
}

void extractTemplateSystVariations(
  process_spec& proc, RooDataSet* data,
  std::unordered_map<TString, std::pair<RooRealVar*, RooRealVar*>> const& zipped_nuisances,
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > const& postfit_vals_map
){
  std::unordered_map<TString, TH1F> dummymap_1D;
  std::unordered_map<TString, TH2F> dummymap_2D;
  std::unordered_map<TString, TH3F> dummymap_3D;
  TString const& procname = proc.name;
  for (auto const& nuisname_nuispair_pair:zipped_nuisances){
    TString nuisname = nuisname_nuispair_pair.first;
    // These are POIs, not nuisances. Skip them.
    if (!isTrueNuisance(nuisname)) continue;

    RooRealVar* nuis_shape = nuisname_nuispair_pair.second.first;
    RooRealVar* nuis_norm = nuisname_nuispair_pair.second.second;
    if (
      !(
        (nuis_shape && proc.dependsOn(nuis_shape))
        ||
        (nuis_norm && proc.dependsOn(nuis_norm))
        )
      ) continue;

    auto const& postfit_vals = postfit_vals_map.find(nuisname)->second;
    double const& vnom = postfit_vals.first;
    double const& vdn = postfit_vals.second.first;
    double const& vup = postfit_vals.second.second;

    TString strhname;

    strhname = Form("%s_%s_Down", procname.Data(), nuisname.Data());
    if (nuis_shape) nuis_shape->setVal(vdn);
    if (nuis_norm) nuis_norm->setVal(vdn);
    //cout << "\t- Extracting " << strhname << " at " << nuisname << "=" << vdn << "..." << endl;
    extractTemplates(proc, data, dummymap_1D, dummymap_2D, dummymap_3D, strhname, true);

    strhname = Form("%s_%s_Up", procname.Data(), nuisname.Data());
    if (nuis_shape) nuis_shape->setVal(vup);
    if (nuis_norm) nuis_norm->setVal(vup);
    //cout << "\t- Extracting " << strhname << " at " << nuisname << "=" << vup << "..." << endl;
    extractTemplates(proc, data, dummymap_1D, dummymap_2D, dummymap_3D, strhname, true);

    if (nuis_shape) nuis_shape->setVal(vnom);
    if (nuis_norm) nuis_norm->setVal(vnom);
  }
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


TGraphAsymmErrors* getDataGraph(TH1F* hdata, bool errorsOnZero, TString strappend=""){
  TGraphAsymmErrors* tgdata = nullptr;
  if (hdata){
    TH1F* htmp = (TH1F*) hdata->Clone("__data_tmp__");
    for (int bin = 1; bin <= htmp->GetNbinsX(); bin++){
      double bincontent = htmp->GetBinContent(bin);
      if (bincontent>=0.) htmp->SetBinError(bin, std::sqrt(bincontent));
      else htmp->SetBinError(bin, 0);
    }
    HelperFunctions::convertTH1FToTGraphAsymmErrors(htmp, tgdata, errorsOnZero, true, false);
    delete htmp;
    tgdata->SetName(TString("tgdata") + (strappend=="" ? "" : "_") + strappend + (errorsOnZero ? "_withZeros" : ""));
    tgdata->SetMarkerSize(1.2);
    tgdata->SetMarkerStyle(20);
    tgdata->SetMarkerColor(kBlack);
    tgdata->SetLineColor(kBlack);
    tgdata->SetLineWidth(1);
  }
  else cout << "Data histogram is null." << endl;
  return tgdata;
}

void getParameterErrors(RooRealVar const& par, double& errLo, double& errHi){
  double errSym = par.getError();
  double errAsym[2]={ errSym, errSym };
  if (par.hasAsymError()){
    errAsym[0] = std::abs(par.getAsymErrorLo()); // This value is negative.
    errAsym[1] = std::abs(par.getAsymErrorHi());
  }
  errLo = errAsym[0];
  errHi = errAsym[1];
}

void readParameters(
  TString fname,
  std::unordered_map<TString, double>& vals_hypo,
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> >& postfit_vals_map
){
  if (fname=="" || !HostHelpers::FileExists(fname)) return;

  TDirectory* curdir = gDirectory;
  TFile* finput_postfit = TFile::Open(fname, "read");
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
      double err_dn, err_up;
      getParameterErrors(*rvar, err_dn, err_up);
      postfit_vals_map[strvar] = std::pair<double, std::pair<double, double>>(val, std::pair<double, double>(val-err_dn, val+err_up));
      if (vals_hypo.find(strvar)!=vals_hypo.end()) vals_hypo[strvar] = val;
    }
  }
  delete it;
  finput_postfit->Close();
  curdir->cd();
}


using namespace HelperFunctions;


void plotHypoLikelihood(
  bool markPreliminary=false,
  bool useLogY=false,
  TString strFitResultFile_Nominal="",
  TString strFitResultFile_GGsmALT="",
  TString strFitResultFile_NoHSigALT="",
  TString strPostfit=""
){
  TString const strAChypo="SM";
  constexpr bool isPostFit = true;
  constexpr bool runSysts = true;

  // Magic numbers
  constexpr double npixels_stdframe_xy = 800;
  constexpr double relmargin_frame_left = 0.20;
  constexpr double relmargin_frame_right = 0.05;
  constexpr double relmargin_frame_CMS = 0.07;
  constexpr double relmargin_frame_XTitle = 0.15;
  constexpr double relmargin_frame_separation = 0.1;
  constexpr double relsize_frame_ratio = 0.2;
  constexpr double relsize_CMSlogo = 0.98;
  constexpr double relsize_CMSlogo_sqrts = 0.8;
  constexpr double relsize_XYTitle = 0.9;
  constexpr double relsize_XYLabel = 0.8;
  constexpr double offset_xlabel = 0.02;
  constexpr double offset_ylabel = 0.02;
  constexpr double offset_xtitle = 1.09;
  constexpr double offset_ytitle = 1.3;

  constexpr unsigned int npads = 2;
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
  const double npixels_pad_top = int(
    npixels_stdframe_xy*(
      relmargin_frame_CMS
      + 1.
      + relmargin_frame_separation/2.
      ) + 0.5
    );
  const double npixels_pad_bot = int(
    npixels_stdframe_xy*(
      relmargin_frame_separation/2.
      + relsize_frame_ratio
      + relmargin_frame_XTitle
      ) + 0.5
    );
  const double npixels_y = npixels_pad_top + npixels_pad_bot;

  gStyle->SetOptStat(0);

  std::vector<TString> const sqrtsnames{ "13TeV_2016", "13TeV_2017", "13TeV_2018" };
  const unsigned int nsqrts = sqrtsnames.size();

  TDirectory* curdir = gDirectory;

  // Determine alternative model parameters
  constexpr double GHref = 4.07;
  unordered_map<TString, double> val_NominalHypo; val_NominalHypo["RV"]=1; val_NominalHypo["RF"]=1; val_NominalHypo["GGsm"]=1; val_NominalHypo["fai1"]=0;

  constexpr bool includeGGSmALT = false;
  unordered_map<TString, double> val_GGsmALT(val_NominalHypo); val_GGsmALT["GGsm"]=20./GHref;

  constexpr bool includeNoHSigALT = true;
  unordered_map<TString, double> val_NoHSigALT(val_NominalHypo); val_NoHSigALT["GGsm"]=0;

  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_Nominal_map;
  readParameters(strFitResultFile_Nominal, val_NominalHypo, postfit_vals_Nominal_map);
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_GGsmALT_map;
  if (includeGGSmALT) readParameters(strFitResultFile_GGsmALT, val_GGsmALT, postfit_vals_GGsmALT_map);
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_NoHSigALT_map;
  if (includeNoHSigALT) readParameters(strFitResultFile_NoHSigALT, val_NoHSigALT, postfit_vals_NoHSigALT_map);

  cout << "Nominal parameters:" << endl;
  for (auto const& pp:postfit_vals_Nominal_map){
    cout << "\t- " << pp.first << " = " << pp.second.first << endl;
  }
  cout << "GGsmALT parameters:" << endl;
  for (auto const& pp:postfit_vals_GGsmALT_map){
    cout << "\t- " << pp.first << " = " << pp.second.first << endl;
  }
  cout << "NoHSigALT parameters:" << endl;
  for (auto const& pp:postfit_vals_NoHSigALT_map){
    cout << "\t- " << pp.first << " = " << pp.second.first << endl;
  }

  // Get process order/labels
  std::vector<TString> proc_order{ "total_BestFit_ZZTo4L", "total_BestFit"/*, "total_GGsmALT"*/, "total_NoHSigALT" };
  std::vector<TString> proc_label{
    "Best fit, 4l component",
    "Best fit, 2l2#nu component",
    //Form("%s#Gamma_{H}=%s MeV", (aihypo=="" ? "" : "f_{ai}=0, "), getFractionString(val_GGsmALT["GGsm"]*GHref).Data()),
    "No off-shell"
  };
  std::vector<int> proc_color{ int(TColor::GetColor("#ff66ff")), int(TColor::GetColor("#66ff99"))/*, int(kCyan+2)*/, int(kOrange-3) };
  std::vector<int> proc_code{ -2, -2/*, -2*/, -2 };
  proc_order.push_back("data_ZZTo4L");
  proc_color.push_back((int) kBlack);
  proc_code.push_back(-99);
  proc_label.push_back("Observed 4l");
  proc_order.push_back("data");
  proc_color.push_back((int) kBlack);
  proc_code.push_back(-99);
  proc_label.push_back("Observed 2l2#nu+4l");

  constexpr double binwidth_adj = 0.001;
  std::vector<double> const binning_HD{
    0.5, 0.502, 0.504, 0.506, 0.508,
    0.51, 0.515, 0.525,
    0.54, 0.565, 0.6
  };
  const int nbins_HD = binning_HD.size()-1;
  unordered_map<TString, TH1F> procshape_flat_translated;
  unordered_map<TString, std::pair<TH1F, TH1F>> syst_totalshape_flat_translated_map;

  std::vector<TString> const strinputdirs{ "Offshell_2L2Nu", "Offshell_4L" };
  for (auto const& cinputdir:strinputdirs){
    TString cinputdir_lower = cinputdir; cinputdir_lower.ToLower();
    bool const is_2l2nu = cinputdir_lower.Contains("2l2nu");
    bool const is_4l = !is_2l2nu;

    std::vector<TString> catnames;
    if (is_4l) catnames = std::vector<TString>{ "Untagged", "JJVBFTagged", "HadVHTagged" };
    else if (is_2l2nu) catnames = std::vector<TString>{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2_pTmiss_lt_200", "Nj_geq_2_pTmiss_ge_200" };
    else return;
    const unsigned int ncats = catnames.size();

    std::vector<TString> channames;
    if (is_2l2nu) channames = std::vector<TString>{ "2e2nu", "2mu2nu" };
    else /*if (is_4l)*/ channames = std::vector<TString>{ "4mu", "4e", "2e2mu" };
    const unsigned int nchans = channames.size();

    for (unsigned int icat=0; icat<ncats; icat++){
      TString const& catname = catnames.at(icat);
      cout << "Acquiring shapes for category " << catname << endl;

      for (unsigned int ic=0; ic<nchans; ic++){
        TString const& channame = channames.at(ic);
        cout << "\t- Channel:  " << channame << endl;
        for (unsigned int is=0; is<nsqrts; is++){
          TString const& sqrtsname = sqrtsnames.at(is);
          cout << "\t- Data period:  " << sqrtsname << endl;
          TString cinput_main = cinputdir + "/SM/" + sqrtsname + "/hto" + channame + "_" + catname;
          TString cinput_file = cinput_main + ".input.root";
          TString cinput_dc = cinput_main + ".txt";
          if (!HostHelpers::FileExists(cinput_file)){ cout << "File " << cinput_file << " does not exist. Skipping this channel..." << endl; continue; }
          else if (!HostHelpers::FileExists(cinput_dc)){ cout << "File " << cinput_dc << " does not exist. Skipping this channel..." << endl; continue; }

          vector<TString> procname;
          vector<double> procrate;
          vector<RooRealVar*> lnNmodvars;
          vector<RooConstVar*> kappavars;
          unordered_map<TString, std::vector<RooAbsReal*>> procname_rateModifier_map;

          unsigned int ndims=0;
          unordered_map<TString, TH1F> procshape_1D;
          unordered_map<TString, TH2F> procshape_2D;
          unordered_map<TString, TH3F> procshape_3D;

          unordered_map<TString, std::pair<TH1F, TH1F>> syst_totalshape_map_1D;
          unordered_map<TString, std::pair<TH2F, TH2F>> syst_totalshape_map_2D;
          unordered_map<TString, std::pair<TH3F, TH3F>> syst_totalshape_map_3D;

          unordered_map<TString, TH1F> procshape_flat;
          unordered_map<TString, std::pair<TH1F, TH1F>> syst_totalshape_flat_map;

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
              splitOptionRecursive(strline, systdist, ' ', false);
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
            std::vector<double>::iterator it_rate = procrate.begin();
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

          std::vector<RooRealVar*> nuisanceVars;

          // Set postfit values first if they are present.
          for (auto const& pp:postfit_vals_Nominal_map){
            RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
            double const& vnom = pp.second.first;
            double const& vlow = pp.second.second.first;
            double const& vhigh = pp.second.second.second;
            if (tmpvar){
              nuisanceVars.push_back(tmpvar);
              //cout << "\t- Setting " << pp.first << " = " << vnom << " [" << vlow << ", " << vhigh << "] in the workspace." << endl;
              tmpvar->setVal(vnom);
              tmpvar->setAsymError(vlow-vnom, vhigh-vnom);
            }

            for (auto const& lnNmodvar:lnNmodvars){
              if (TString(lnNmodvar->GetName()) == pp.first){
                //cout << "\t- Setting " << pp.first << " = " << vnom << " [" << vlow << ", " << vhigh << "] in the list of lnN nuisances." << endl;
                lnNmodvar->setVal(vnom);
                lnNmodvar->setAsymError(vlow-vnom, vhigh-vnom);
              }
            }
          }
          std::unordered_map<TString, std::pair<RooRealVar*, RooRealVar*>> zipped_nuisances = getZippedNuisanceMap(nuisanceVars, lnNmodvars);

          unordered_map<TString, RooRealVar*> controlVars;
          std::vector<TString> varsToCheck{
            "R", "RF", "RF_13TeV", "RV", "RV_13TeV", "R_13TeV", "CMS_zz4l_fai1", "GGsm", "kbkg_gg", "kbkg_VBF"
          };
          for (auto const& strvar:varsToCheck){
            RooRealVar* tmpvar = ws->var(strvar);
            controlVars[strvar] = tmpvar;
            if (tmpvar){
              if (strvar!="CMS_zz4l_fai1"){
                tmpvar->setVal(1);
                if (strvar=="GGsm") tmpvar->removeRange();
              }
              else tmpvar->setVal(0);
            }
            else cerr << strvar << " could not be found!" << endl;
          }
          RooDataSet* data = (RooDataSet*) ws->data("data_obs");

          // FIXME: HACK lumi 2016 here
          {
            RooRealVar* var_lumi_13TeV_2016 = dynamic_cast<RooRealVar*>(ws->var("LUMI_13TeV_2016"));
            if (var_lumi_13TeV_2016) var_lumi_13TeV_2016->setVal(36.326450);
          }

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

            if (includeGGSmALT){
              setControlVariableValue(controlVars, "GGsm", val_GGsmALT["GGsm"]);
              setControlVariableValue(controlVars, "RV", val_GGsmALT["RV"]);
              setControlVariableValue(controlVars, "RF", val_GGsmALT["RF"]);
              for (auto const& pp:postfit_vals_GGsmALT_map){
                RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
                double const& vnom = pp.second.first;
                double const& vlow = pp.second.second.first;
                double const& vhigh = pp.second.second.second;
                if (tmpvar){
                  tmpvar->setVal(vnom);
                  tmpvar->setAsymError(vlow-vnom, vhigh-vnom);
                }
                for (auto const& lnNmodvar:lnNmodvars){
                  if (TString(lnNmodvar->GetName()) == pp.first){
                    lnNmodvar->setVal(vnom);
                    lnNmodvar->setAsymError(vlow-vnom, vhigh-vnom);
                  }
                }
              }
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_GGsmALT");
              setControlVariableValue(controlVars, "GGsm", val_NominalHypo["GGsm"]);
              setControlVariableValue(controlVars, "RV", val_NominalHypo["RV"]);
              setControlVariableValue(controlVars, "RF", val_NominalHypo["RF"]);
              for (auto const& pp:postfit_vals_Nominal_map){
                RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
                double const& vnom = pp.second.first;
                double const& vlow = pp.second.second.first;
                double const& vhigh = pp.second.second.second;
                if (tmpvar){
                  tmpvar->setVal(vnom);
                  tmpvar->setAsymError(vlow-vnom, vhigh-vnom);
                }
                for (auto const& lnNmodvar:lnNmodvars){
                  if (TString(lnNmodvar->GetName()) == pp.first){
                    lnNmodvar->setVal(vnom);
                    lnNmodvar->setAsymError(vlow-vnom, vhigh-vnom);
                  }
                }
              }
            }
            if (includeNoHSigALT){
              setControlVariableValue(controlVars, "GGsm", val_NoHSigALT["GGsm"]);
              setControlVariableValue(controlVars, "RV", val_NoHSigALT["RV"]);
              setControlVariableValue(controlVars, "RF", val_NoHSigALT["RF"]);
              for (auto const& pp:postfit_vals_NoHSigALT_map){
                RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
                double const& vnom = pp.second.first;
                double const& vlow = pp.second.second.first;
                double const& vhigh = pp.second.second.second;
                if (tmpvar){
                  tmpvar->setVal(vnom);
                  tmpvar->setAsymError(vlow-vnom, vhigh-vnom);
                }
                for (auto const& lnNmodvar:lnNmodvars){
                  if (TString(lnNmodvar->GetName()) == pp.first){
                    lnNmodvar->setVal(vnom);
                    lnNmodvar->setAsymError(vlow-vnom, vhigh-vnom);
                  }
                }
              }
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_NoHSigALT");
              setControlVariableValue(controlVars, "GGsm", val_NominalHypo["GGsm"]);
              setControlVariableValue(controlVars, "RV", val_NominalHypo["RV"]);
              setControlVariableValue(controlVars, "RF", val_NominalHypo["RF"]);
              for (auto const& pp:postfit_vals_Nominal_map){
                RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
                double const& vnom = pp.second.first;
                double const& vlow = pp.second.second.first;
                double const& vhigh = pp.second.second.second;
                if (tmpvar){
                  tmpvar->setVal(vnom);
                  tmpvar->setAsymError(vlow-vnom, vhigh-vnom);
                }
                for (auto const& lnNmodvar:lnNmodvars){
                  if (TString(lnNmodvar->GetName()) == pp.first){
                    lnNmodvar->setVal(vnom);
                    lnNmodvar->setAsymError(vlow-vnom, vhigh-vnom);
                  }
                }
              }
            }

            if (is_4l) ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_BestFit_ZZTo4L");
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_BestFit");

            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D);
            if (runSysts) extractTemplateSystVariations(procSpecs[pname], data, zipped_nuisances, postfit_vals_Nominal_map);
          }
          extractDataTemplates(procSpecs[procname.front()], data, procshape_1D, procshape_2D, procshape_3D, "data");
          if (is_4l) extractDataTemplates(procSpecs[procname.front()], data, procshape_1D, procshape_2D, procshape_3D, "data_ZZTo4L");

          {
            for (auto const& pp:procSpecs){
              auto const& proc = pp.second;
              auto const& pname = proc.name;
              switch (ndims){
              case 1:
              {
                auto const& hist_nom = proc.syst_h1D_map.find(pname)->second;
                for (auto const& hhpp:proc.syst_h1D_map){
                  TString systname = hhpp.first;
                  if (!systname.Contains("Down") && !systname.Contains("Up")) continue;
                  auto const& syst_hist = hhpp.second;
                  bool const isDn = systname.Contains("Down");
                  replaceString<TString, TString const>(systname, "_Down", "");
                  replaceString<TString, TString const>(systname, "_Up", "");
                  replaceString<TString, TString const>(systname, Form("%s_", pname.Data()), "");
                  auto it_syst_totalshape_map = syst_totalshape_map_1D.find(systname);
                  if (it_syst_totalshape_map==syst_totalshape_map_1D.end()){
                    TH1F* htmp_dn = (TH1F*) syst_hist->Clone(systname + "_Down"); htmp_dn->Reset("ICESM");
                    TH1F* htmp_up = (TH1F*) syst_hist->Clone(systname + "_Up"); htmp_up->Reset("ICESM");
                    syst_totalshape_map_1D[systname] = std::pair<TH1F, TH1F>(*htmp_dn, *htmp_up);
                    delete htmp_dn;
                    delete htmp_up;
                    it_syst_totalshape_map = syst_totalshape_map_1D.find(systname);
                  }
                  auto& target_hist = (isDn ? it_syst_totalshape_map->second.first : it_syst_totalshape_map->second.second);
                  target_hist.Add(syst_hist, 1.);
                  target_hist.Add(hist_nom, -1.);
                }
                break;
              }
              case 2:
              {
                auto const& hist_nom = proc.syst_h2D_map.find(pname)->second;
                for (auto const& hhpp:proc.syst_h2D_map){
                  TString systname = hhpp.first;
                  if (!systname.Contains("Down") && !systname.Contains("Up")) continue;
                  auto const& syst_hist = hhpp.second;
                  bool const isDn = systname.Contains("Down");
                  replaceString<TString, TString const>(systname, "_Down", "");
                  replaceString<TString, TString const>(systname, "_Up", "");
                  replaceString<TString, TString const>(systname, Form("%s_", pname.Data()), "");
                  auto it_syst_totalshape_map = syst_totalshape_map_2D.find(systname);
                  if (it_syst_totalshape_map==syst_totalshape_map_2D.end()){
                    TH2F* htmp_dn = (TH2F*) syst_hist->Clone(systname + "_Down"); htmp_dn->Reset("ICESM");
                    TH2F* htmp_up = (TH2F*) syst_hist->Clone(systname + "_Up"); htmp_up->Reset("ICESM");
                    syst_totalshape_map_2D[systname] = std::pair<TH2F, TH2F>(*htmp_dn, *htmp_up);
                    delete htmp_dn;
                    delete htmp_up;
                    it_syst_totalshape_map = syst_totalshape_map_2D.find(systname);
                  }
                  auto& target_hist = (isDn ? it_syst_totalshape_map->second.first : it_syst_totalshape_map->second.second);
                  target_hist.Add(syst_hist, 1.);
                  target_hist.Add(hist_nom, -1.);
                }
                break;
              }
              case 3:
              {
                auto const& hist_nom = proc.syst_h3D_map.find(pname)->second;
                for (auto const& hhpp:proc.syst_h3D_map){
                  TString systname = hhpp.first;
                  if (!systname.Contains("Down") && !systname.Contains("Up")) continue;
                  auto const& syst_hist = hhpp.second;
                  bool const isDn = systname.Contains("Down");
                  replaceString<TString, TString const>(systname, "_Down", "");
                  replaceString<TString, TString const>(systname, "_Up", "");
                  replaceString<TString, TString const>(systname, Form("%s_", pname.Data()), "");
                  auto it_syst_totalshape_map = syst_totalshape_map_3D.find(systname);
                  if (it_syst_totalshape_map==syst_totalshape_map_3D.end()){
                    TH3F* htmp_dn = (TH3F*) syst_hist->Clone(systname + "_Down"); htmp_dn->Reset("ICESM");
                    TH3F* htmp_up = (TH3F*) syst_hist->Clone(systname + "_Up"); htmp_up->Reset("ICESM");
                    syst_totalshape_map_3D[systname] = std::pair<TH3F, TH3F>(*htmp_dn, *htmp_up);
                    delete htmp_dn;
                    delete htmp_up;
                    it_syst_totalshape_map = syst_totalshape_map_3D.find(systname);
                  }
                  auto& target_hist = (isDn ? it_syst_totalshape_map->second.first : it_syst_totalshape_map->second.second);
                  target_hist.Add(syst_hist, 1.);
                  target_hist.Add(hist_nom, -1.);
                }
                break;
              }
              default:
                break;
              }
            }
          }

          int nbins=0;
          for (auto const& pp:procshape_1D){
            TH1F* hflat = flattenHistogram(&(pp.second), false);
            if (nbins==0) nbins = hflat->GetNbinsX();
            procshape_flat[pp.first] = *hflat;
            delete hflat;
          }
          for (auto const& pp:procshape_2D){
            TH1F* hflat = flattenHistogram(&(pp.second), false);
            if (nbins==0) nbins = hflat->GetNbinsX();
            procshape_flat[pp.first] = *hflat;
            delete hflat;
          }
          for (auto const& pp:procshape_3D){
            TH1F* hflat = flattenHistogram(&(pp.second), false);
            if (nbins==0) nbins = hflat->GetNbinsX();
            procshape_flat[pp.first] = *hflat;
            delete hflat;
          }
          for (auto const& pp:syst_totalshape_map_1D){
            TH1F* hflat_dn = flattenHistogram(&(pp.second.first), false);
            TH1F* hflat_up = flattenHistogram(&(pp.second.second), false);
            syst_totalshape_flat_map[pp.first] = std::pair<TH1F, TH1F>(*hflat_dn, *hflat_up);
            delete hflat_dn;
            delete hflat_up;
          }
          for (auto const& pp:syst_totalshape_map_2D){
            TH1F* hflat_dn = flattenHistogram(&(pp.second.first), false);
            TH1F* hflat_up = flattenHistogram(&(pp.second.second), false);
            syst_totalshape_flat_map[pp.first] = std::pair<TH1F, TH1F>(*hflat_dn, *hflat_up);
            delete hflat_dn;
            delete hflat_up;
          }
          for (auto const& pp:syst_totalshape_map_3D){
            TH1F* hflat_dn = flattenHistogram(&(pp.second.first), false);
            TH1F* hflat_up = flattenHistogram(&(pp.second.second), false);
            syst_totalshape_flat_map[pp.first] = std::pair<TH1F, TH1F>(*hflat_dn, *hflat_up);
            delete hflat_dn;
            delete hflat_up;
          }

          TH1F hHTestVals("HTestVals", "", nbins, 0, nbins);
          TH1F const& hBestFit = procshape_flat.find("total_BestFit")->second;
          TH1F const& hBkgOnly = procshape_flat.find("total_NoHSigALT")->second;
          for (int ix=1; ix<=nbins; ix++){
            double bc_bf = hBestFit.GetBinContent(ix);
            double bc_bo = hBkgOnly.GetBinContent(ix);
            double bc_hd = -1;
            if ((bc_bo+bc_bf)>0.) bc_hd = bc_bo/(bc_bf + bc_bo);
            hHTestVals.SetBinContent(ix, bc_hd);
          }

          for (auto const& pp:procshape_flat){
            auto const& pname = pp.first;

            cout << pname << " integral before = " << getHistogramIntegralAndError(&(pp.second), 1, nbins, false, nullptr) << " ?= ";
            double integral_after = 0;
            double integral_translated = 0;

            if (procshape_flat_translated.find(pname)==procshape_flat_translated.end()) procshape_flat_translated[pname] = TH1F(Form("%s_translated", pp.second.GetName()), "", nbins_HD, binning_HD.data());
            integral_translated = getHistogramIntegralAndError(&(procshape_flat_translated[pname]), 1, nbins_HD, false, nullptr);
            for (int ix=1; ix<=nbins; ix++){
              double bc_add = pp.second.GetBinContent(ix); integral_after += bc_add;
              double xc = hHTestVals.GetBinContent(ix);
              int jx = procshape_flat_translated[pname].GetXaxis()->FindBin(xc);
              jx = std::max(1, std::min(nbins_HD, jx));
              double bc_old = procshape_flat_translated[pname].GetBinContent(jx);
              procshape_flat_translated[pname].SetBinContent(jx, bc_add + bc_old);
            }
            integral_translated = getHistogramIntegralAndError(&(procshape_flat_translated[pname]), 1, nbins_HD, false, nullptr) - integral_translated;

            cout << integral_after << " ?= " << integral_translated << endl;
          }
          for (auto const& pp:syst_totalshape_flat_map){
            auto const& pname = pp.first;
            if (syst_totalshape_flat_translated_map.find(pname)==syst_totalshape_flat_translated_map.end()){
              syst_totalshape_flat_translated_map[pname] = std::pair<TH1F, TH1F>(
                TH1F(Form("%s_translated", pp.second.first.GetName()), "", nbins_HD, binning_HD.data()),
                TH1F(Form("%s_translated", pp.second.second.GetName()), "", nbins_HD, binning_HD.data())
                );
            }
            for (int ix=1; ix<=nbins; ix++){
              double bc_add = pp.second.first.GetBinContent(ix);
              double xc = hHTestVals.GetBinContent(ix);
              int jx = syst_totalshape_flat_translated_map[pname].first.GetXaxis()->FindBin(xc);
              jx = std::max(1, std::min(nbins_HD, jx));
              double bc_old = syst_totalshape_flat_translated_map[pname].first.GetBinContent(jx);
              syst_totalshape_flat_translated_map[pname].first.SetBinContent(jx, bc_add + bc_old);
            }
            for (int ix=1; ix<=nbins; ix++){
              double bc_add = pp.second.second.GetBinContent(ix);
              double xc = hHTestVals.GetBinContent(ix);
              int jx = syst_totalshape_flat_translated_map[pname].second.GetXaxis()->FindBin(xc);
              jx = std::max(1, std::min(nbins_HD, jx));
              double bc_old = syst_totalshape_flat_translated_map[pname].second.GetBinContent(jx);
              syst_totalshape_flat_translated_map[pname].second.SetBinContent(jx, bc_add + bc_old);
            }
          }

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
    }
  }

  TGraphAsymmErrors* tgdata = nullptr;
  TGraphAsymmErrors* tgdata_4l = nullptr;
  for (auto& pp:procshape_flat_translated){
    auto const& procname = pp.first;
    TString procname_lower = procname; procname_lower.ToLower();
    TH1F* tpl = &(pp.second);
    TGraphAsymmErrors* tgtmp = nullptr;
    if (procname.BeginsWith("data")){
      tgtmp = getDataGraph(tpl, true, procname);
      for (int ix=0; ix<tgtmp->GetN(); ix++){
        double bw = tpl->GetXaxis()->GetBinWidth(tpl->GetXaxis()->FindBin(tgtmp->GetX()[ix]))/binwidth_adj;
        tgtmp->GetY()[ix] /= bw;
        tgtmp->GetEYhigh()[ix] /= bw;
        tgtmp->GetEYlow()[ix] /= bw;
      }
      if (procname_lower.Contains("4l")) tgdata_4l = tgtmp;
      else tgdata = tgtmp;
    }
  }
  if (tgdata && tgdata_4l){
    for (int ix=0; ix<tgdata->GetN(); ix++){
      double bc_total_low = tgdata->GetY()[ix]-std::abs(tgdata->GetEYlow()[ix]);
      double bc_4l_high = tgdata_4l->GetY()[ix]+tgdata_4l->GetEYhigh()[ix];
      if (bc_4l_high>=bc_total_low){
        tgdata->GetX()[ix] += binwidth_adj*0.5;
        tgdata_4l->GetX()[ix] -= binwidth_adj*0.5;
      }
    }
  }
  // Account for bin widths
  for (auto& pp:procshape_flat_translated){
    TH1F* tpl = &(pp.second);
    double integral_before = getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), false, nullptr);
    double count_lastbin = tpl->GetBinContent(tpl->GetNbinsX());
    divideBinWidth(tpl);
    pp.second.Scale(binwidth_adj);
    double count_lastbin_after = tpl->GetBinContent(tpl->GetNbinsX());
    double integral_after = getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), true, nullptr)/binwidth_adj;
    cout << pp.first << " final integral = " << integral_before << " ?= " << integral_after << ". Last bin count = " << count_lastbin << " ?!= " << count_lastbin_after << "." << endl;
  }
  for (auto& pp:syst_totalshape_flat_translated_map){
    divideBinWidth(&(pp.second.first));
    divideBinWidth(&(pp.second.second));
    pp.second.first.Scale(binwidth_adj);
    pp.second.second.Scale(binwidth_adj);
  }

  curdir->cd();

  std::pair<TH1F*, TH1F*> systband;
  TH1F* hsyst_dn = new TH1F("allprocs_syst_dn", "", nbins_HD, binning_HD.data()); systband.first = hsyst_dn;
  TH1F* hsyst_up = new TH1F("allprocs_syst_up", "", nbins_HD, binning_HD.data()); systband.second = hsyst_up;
  for (auto& pp:syst_totalshape_flat_translated_map){
    TH1F* hist_dn = &(pp.second.first);
    TH1F* hist_up = &(pp.second.second);

    for (int ix=1; ix<=nbins_HD; ix++){
      double bc_add_dn = hist_dn->GetBinContent(ix);
      double bc_dn = hsyst_dn->GetBinContent(ix);
      double bc_add_up = hist_up->GetBinContent(ix);
      double bc_up = hsyst_up->GetBinContent(ix);
      if (bc_add_dn>bc_add_up) std::swap(bc_add_dn, bc_add_up);
      bc_add_dn = std::min(bc_add_dn, 0.);
      bc_add_up = std::max(bc_add_up, 0.);
      bc_dn = -sqrt(std::pow(bc_dn, 2) + std::pow(bc_add_dn, 2));
      bc_up = sqrt(std::pow(bc_up, 2) + std::pow(bc_add_up, 2));
      hsyst_dn->SetBinContent(ix, bc_dn);
      hsyst_up->SetBinContent(ix, bc_up);
    }
  }

  // Make a canvas for the legend
  {
    TH1F* hsum_nominal = &(procshape_flat_translated.find("total_BestFit")->second);
    std::vector<TH1F*> hlist_ALT;
    for (auto& pp:procshape_flat_translated){
      if (pp.first.Contains("ALT")) hlist_ALT.push_back(&(pp.second));
    }
    TGraphAsymmErrors* tg_systband = nullptr;
    {
      std::vector<double> xx, xl, xh, yy, yl, yh;
      for (int ix=1; ix<=hsum_nominal->GetNbinsX(); ix++){
        xx.push_back(hsum_nominal->GetXaxis()->GetBinCenter(ix));
        xl.push_back(-hsum_nominal->GetXaxis()->GetBinLowEdge(ix) + xx.back());
        xh.push_back(hsum_nominal->GetXaxis()->GetBinUpEdge(ix) - xx.back());
        yy.push_back(hsum_nominal->GetBinContent(ix));
        yl.push_back(-systband.first->GetBinContent(ix));
        yh.push_back(systband.second->GetBinContent(ix));
      }
      tg_systband = new TGraphAsymmErrors(xx.size(), xx.data(), yy.data(), xl.data(), xh.data(), yl.data(), yh.data());
      tg_systband->SetName("tg_systband");
      tg_systband->SetTitle("");
      tg_systband->SetLineColor(kBlack);
      tg_systband->SetMarkerColor(kBlack);
      tg_systband->SetFillColor(kBlack);
      tg_systband->SetFillStyle(3345);
    }

    // Draw
    cout << "Creating the canvas..." << endl;
    TCanvas canvas(
      TString("c_HypoLL") + TString((markPreliminary ? "_Preliminary" : "")) + "_SM" + (!useLogY ? "" : "_LogY"),
      "", 8, 30, npixels_x, npixels_y
    );
    canvas.cd();

    int nleft=0;
    int nright=0;
    if (tgdata) nleft++;
    if (tgdata_4l) nleft++;
    if (procshape_flat_translated.find("total_BestFit")!=procshape_flat_translated.end()) nright++;
    if (procshape_flat_translated.find("total_BestFit_ZZTo4L")!=procshape_flat_translated.end()) nright++;
    if (procshape_flat_translated.find("total_NoHSigALT")!=procshape_flat_translated.end()) nright++;


    double leg_xmin=(relmargin_frame_left + 0.03)/(1. + relmargin_frame_left + relmargin_frame_right);
    double leg_xmax=(relmargin_frame_left + 0.4)/(1. + relmargin_frame_left + relmargin_frame_right);
    double leg_ymax=(npixels_y - npixels_stdframe_xy*(0.05 + relmargin_frame_CMS))/npixels_y;
    double leg_ymin=leg_ymax; leg_ymin -= npixels_XYTitle*nleft/npixels_y*1.25;
    TLegend legend_left(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
    legend_left.SetBorderSize(0);
    legend_left.SetTextFont(43);
    legend_left.SetTextSize(npixels_XYTitle);
    legend_left.SetLineColor(1);
    legend_left.SetLineStyle(1);
    legend_left.SetLineWidth(1);
    legend_left.SetFillColor(0);
    legend_left.SetFillStyle(0);

    leg_xmax = (0.97 + relmargin_frame_left)/(1. + relmargin_frame_left + relmargin_frame_right);
    leg_xmin = (0.13 + relmargin_frame_left)/(1. + relmargin_frame_left + relmargin_frame_right);
    leg_ymax = leg_ymin - npixels_XYTitle/npixels_y*0.5;
    leg_ymin=leg_ymax; leg_ymin -= npixels_XYTitle*nright/npixels_y*1.25;
    TLegend legend_right(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
    legend_right.SetBorderSize(0);
    legend_right.SetTextFont(43);
    legend_right.SetTextSize(npixels_XYTitle);
    legend_right.SetLineColor(1);
    legend_right.SetLineStyle(1);
    legend_right.SetLineWidth(1);
    legend_right.SetFillColor(0);
    legend_right.SetFillStyle(0);

    cout << "\t- Creating the pads..." << endl;
    TPad* pad_main = nullptr;
    TPad* pad_ratio = nullptr;
    std::vector<TPad*> pads;
    canvas.cd();
    pads.push_back(
      new TPad(
        Form("%s_top", canvas.GetName()), "",
        0, (1.-npixels_pad_top/npixels_y), 1, 1
      )
    );
    pad_main = pads.back();
    canvas.cd();
    pads.push_back(
      new TPad(
        Form("%s_bot", canvas.GetName()), "",
        0, 0, 1, npixels_pad_bot/npixels_y
      )
    );
    pad_ratio = pads.back();
    for (auto& pad:pads){
      pad->cd();
      pad->SetFillColor(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(2);
      pad->SetTickx(1);
      pad->SetTicky(1);
      pad->SetLeftMargin(relmargin_frame_left/(1.+relmargin_frame_left+relmargin_frame_right));
      pad->SetRightMargin(relmargin_frame_right/(1.+relmargin_frame_left+relmargin_frame_right));
      if (pad==pad_main){
        pad->SetTopMargin(npixels_stdframe_xy*relmargin_frame_CMS/npixels_pad_top);
        pad->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_separation/2./npixels_pad_top);
      }
      else{
        pad->SetTopMargin(npixels_stdframe_xy*relmargin_frame_separation/2./npixels_pad_bot);
        pad->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_XTitle/npixels_pad_bot);
      }
      pad->SetFrameFillStyle(0);
      pad->SetFrameBorderMode(0);
      pad->SetFrameFillStyle(0);
      pad->SetFrameBorderMode(0);
      if (useLogY && pad==pad_main) pad->SetLogy(true);
      canvas.cd();
    }
    // Draw in reverse order
    for (auto rit=pads.rbegin(); rit!=pads.rend(); rit++) (*rit)->Draw();

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
    cout << "Size of the CMS logo: " << npixels_CMSlogo << endl;
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
    TString cErgTev = Form("#leq138 fb^{-1} (%i TeV)", theSqrts);
    text = pt.AddText(0.999, 0.45, cErgTev);
    text->SetTextFont(43);
    text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
    text->SetTextAlign(32);

    TPaveText ptc(
      (0.43 + relmargin_frame_left)/(1. + relmargin_frame_left + relmargin_frame_right),
      leg_ymin - npixels_XYTitle/npixels_y*0.5 - npixels_XYTitle/npixels_y*1.25,
      (0.97 + relmargin_frame_left)/(1. + relmargin_frame_left + relmargin_frame_right),
      leg_ymin - npixels_XYTitle/npixels_y*0.5,
      "brNDC"
    );
    ptc.SetBorderSize(0);
    ptc.SetFillStyle(0);
    ptc.SetTextAlign(12);
    ptc.SetTextFont(43);
    ptc.SetTextSize(npixels_XYTitle);
    ptc.AddText(0.001, 0.5, Form("Best fit #Gamma_{H}: %s MeV", castValueToString(val_NominalHypo["GGsm"]*GHref, 1).data()));

    double ymax=-1, ymin=-1;
    if (useLogY) ymin=9e9;
    for (unsigned int ip=0; ip<proc_order.size(); ip++){
      TString const& procname = proc_order.at(ip);
      TString procname_lower = procname; procname_lower.ToLower();
      TString const& proclabel = proc_label.at(ip);
      cout << "Adjusting process " << procname << " at index " << ip << endl;
      TH1F* prochist = &(procshape_flat_translated.find(procname)->second);
      cout << "\t- Process " << procname << " color: " << proc_color[ip] << endl;
      prochist->SetLineWidth(1);
      if (procname.Contains("ALT")){
        if (procname.Contains("GGsm")) prochist->SetLineStyle(2);
        else if (procname.Contains("NoHSig")) prochist->SetLineStyle(4);
        prochist->SetMarkerColor(proc_color[ip]);
        prochist->SetLineColor(proc_color[ip]);
        prochist->SetLineWidth(4);
      }
      else if (proc_code.at(ip)==-99){ // Data
        prochist->SetMarkerColor(proc_color[ip]);
      }
      else{
        prochist->SetMarkerColor(proc_color[ip]);
        prochist->SetLineColorAlpha(proc_color[ip]+1, 0.5);
        prochist->SetFillColor(proc_color[ip]);
        prochist->SetFillStyle(1001);
      }

      prochist->GetXaxis()->SetNdivisions(505);
      prochist->GetXaxis()->SetLabelFont(43);
      prochist->GetXaxis()->SetLabelOffset(offset_xlabel);
      prochist->GetXaxis()->SetLabelSize(npixels_XYLabel);
      prochist->GetXaxis()->SetTitleFont(42);
      prochist->GetXaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_top);
      prochist->GetXaxis()->SetTitleOffset(offset_xtitle);
      prochist->GetYaxis()->SetLabelFont(43);
      prochist->GetYaxis()->SetLabelOffset(offset_ylabel);
      prochist->GetYaxis()->SetLabelSize(npixels_XYLabel);
      prochist->GetYaxis()->SetTitleFont(42);
      prochist->GetYaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_top);
      prochist->GetYaxis()->SetTitleOffset(offset_ytitle);
      prochist->GetYaxis()->SetTitle(Form("Events / %.3f", binwidth_adj));
      prochist->GetYaxis()->CenterTitle();
      prochist->GetXaxis()->CenterTitle();

      if (!procname.Contains("data")){
        for (int ix=1; ix<=prochist->GetNbinsX(); ix++){
          double bc = prochist->GetBinContent(ix);
          if (bc!=0.){
            ymax = std::max(bc, ymax);
            ymin = std::min(bc, ymin);
          }
        }
      }
      else{
        cout << "\t- Obtaining " << procname << " graph" << endl;
        TGraphAsymmErrors* tgdata_tmp = (procname_lower.Contains("4l") ? tgdata_4l : tgdata);
        if (tgdata_tmp){
          for (int ipoint=0; ipoint<tgdata_tmp->GetN(); ipoint++){
            double bc = tgdata_tmp->GetY()[ipoint]+tgdata_tmp->GetEYhigh()[ipoint];
            double bc_low = tgdata_tmp->GetY()[ipoint]-std::abs(tgdata_tmp->GetEYlow()[ipoint]);
            if (bc!=0.){
              ymax = std::max(bc, ymax);
            }
            else{
              ymin = std::min(bc, ymin);
            }
            if (bc_low!=0.){
              ymin = std::min(bc_low, ymin);
            }
          }
        }
      }
    }
    if (tg_systband){
      for (int ipoint=0; ipoint<tg_systband->GetN(); ipoint++){
        double bc = tg_systband->GetY()[ipoint]+tg_systband->GetEYhigh()[ipoint];
        double bc_low = tg_systband->GetY()[ipoint]-std::abs(tg_systband->GetEYlow()[ipoint]);
        ymax = std::max(bc, ymax);
        ymin = std::min(bc_low, ymin);
      }
    }
    if (tgdata_4l){
      tgdata_4l->SetMarkerStyle(24);
      tgdata_4l->SetLineStyle(7);
    }

    if (!useLogY) ymin=0;

    double ymaxfactor = (!useLogY ? 1.25 : 2.);
    double yminfactor = 0.8;

    // Plot the ratio panel
    TGraphAsymmErrors* tg_systband_unit = nullptr;
    TGraphAsymmErrors* tgdata_unit = nullptr;
    TH1F* hdummy_ratio = nullptr;
    std::vector<TH1F*> hlist_ALT_dummy;
    {
      if (tgdata->GetN()!=tg_systband->GetN()) cerr << "Number of bins for data with zeros and syst. band are not the same." << endl;
      else{
        hdummy_ratio = (TH1F*) hsum_nominal->Clone("hframe_ratio"); hdummy_ratio->Reset("ICESM"); hdummy_ratio->SetLineColor(0); hdummy_ratio->SetMarkerColor(0); hdummy_ratio->SetLineWidth(1);
        tg_systband_unit = (TGraphAsymmErrors*) tg_systband->Clone(Form("%s_unit", tg_systband->GetName()));
        tgdata_unit = (TGraphAsymmErrors*) tgdata->Clone(Form("%s_unit", tgdata->GetName()));
        for (auto const& htmp:hlist_ALT) hlist_ALT_dummy.push_back((TH1F*) htmp->Clone(Form("%s_unit", htmp->GetName())));

        double ymin_ratio = 9e9;
        double ymax_ratio = -9e9;
        for (int ix=0; ix<tgdata->GetN(); ix++){
          double& val_systband = tg_systband_unit->GetY()[ix];
          double& val_systband_errdn = tg_systband_unit->GetEYlow()[ix];
          double& val_systband_errup = tg_systband_unit->GetEYhigh()[ix];
          double& val_data = tgdata_unit->GetY()[ix];
          double& val_data_errdn = tgdata_unit->GetEYlow()[ix];
          double& val_data_errup = tgdata_unit->GetEYhigh()[ix];

          for (auto& hh:hlist_ALT_dummy){
            hh->SetBinContent(ix+1, (val_systband!=0. ? hh->GetBinContent(ix+1)/val_systband : 0.));
            hh->SetBinError(ix+1, 0);
          }

          if (val_systband!=0.){
            val_data_errdn /= val_systband;
            val_data_errup /= val_systband;
            val_data /= val_systband;
            val_systband_errdn /= val_systband;
            val_systband_errup /= val_systband;
          }
          else{
            val_data_errdn = 0.;
            val_data_errup = 0.;
            val_data = 0.;
            val_systband_errdn = 0.;
            val_systband_errup = 0.;
          }
          val_systband = 1.;

          double maxthr_data = val_data + (val_data>0. ? std::abs(val_data_errup) : 0.);
          if (maxthr_data>5. && val_data<5.) maxthr_data = val_data;
          else if (maxthr_data>10. && val_data>=5. && val_data<10.) maxthr_data = val_data;
          ymin_ratio = std::min(ymin_ratio, val_systband - std::abs(val_systband_errdn));
          ymin_ratio = std::min(ymin_ratio, val_data - std::abs(val_data_errdn));
          ymax_ratio = std::max(ymax_ratio, val_systband + std::abs(val_systband_errup));
          ymax_ratio = std::max(ymax_ratio, maxthr_data);

          for (auto& hh:hlist_ALT_dummy){
            double val_ALT = hh->GetBinContent(ix+1);
            ymin_ratio = std::min(ymin_ratio, val_ALT);
            ymax_ratio = std::max(ymax_ratio, val_ALT);
          }
        }
        pad_ratio->cd();

        if (2.-ymax_ratio>0.){
          ymax_ratio -= 1.;
          ymin_ratio = 1.-ymin_ratio;

          ymin_ratio = ymax_ratio = std::max(ymin_ratio, ymax_ratio);
          ymin_ratio = 1. - ymin_ratio;
          ymax_ratio += 1.;
        }

        {
          double dy_ratio = ymax_ratio - ymin_ratio;
          ymax_ratio += dy_ratio/20.;
          if (ymin_ratio>dy_ratio/20.) ymin_ratio -= dy_ratio/20.;
          else ymin_ratio = 0.;

          ymax_ratio = static_cast<double>(static_cast<int>(ymax_ratio*10.)+1)/10.;
          ymin_ratio = static_cast<double>(static_cast<int>(ymin_ratio*10.))/10.;
        }

        hdummy_ratio->SetMinimum(ymin_ratio);
        hdummy_ratio->SetMaximum(ymax_ratio);
        hdummy_ratio->GetXaxis()->SetTitleFont(42);
        hdummy_ratio->GetXaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_bot);
        hdummy_ratio->GetXaxis()->SetTitle("N_{no off-shell} / (N_{#lower[-0.25]{no off-shell}} + N_{best fit})");
        hdummy_ratio->GetXaxis()->SetTitleOffset(offset_xtitle);
        hdummy_ratio->GetYaxis()->SetRangeUser(ymin_ratio, ymax_ratio);
        hdummy_ratio->GetYaxis()->SetTitleFont(42);
        hdummy_ratio->GetYaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_bot);
        hdummy_ratio->GetYaxis()->SetTitle("#times#frac{1}{Best fit}");
        hdummy_ratio->GetYaxis()->SetTitleOffset(offset_ytitle*(npixels_pad_bot/npixels_pad_top));
        hdummy_ratio->GetYaxis()->SetNdivisions(3, 5, 2, true);
        hdummy_ratio->GetYaxis()->SetMaxDigits(2);

        hdummy_ratio->Draw("hist");
        for (auto& hh:hlist_ALT_dummy) hh->Draw("histsame");
        tg_systband_unit->Draw("2same");
        tgdata_unit->Draw("0psame");
      }
    }

    canvas.cd();
    pad_main->cd();
    bool drawfirst=true;
    for (unsigned int ip=proc_order.size(); ip>0; ip--){
      TString const& procname = proc_order.at(ip-1);
      TString procname_lower = procname; procname_lower.ToLower();
      TString const& proclabel = proc_label.at(ip-1);
      TH1F* prochist = &(procshape_flat_translated.find(procname)->second);

      if (!procname.Contains("data")){
        if (procname.Contains("ALT")) legend_right.AddEntry(prochist, proclabel, "l");
        else legend_right.AddEntry(prochist, proclabel, "f");
      }
      else{
        if (procname_lower.Contains("4l")) legend_left.AddEntry(tgdata_4l, proclabel, "e1p");
        else legend_left.AddEntry(tgdata, proclabel, "e1p");
      }
      if (procname.Contains("data")) continue;
      prochist->GetXaxis()->SetTitle("");
      prochist->GetXaxis()->SetLabelSize(0);
      prochist->GetXaxis()->SetTitleSize(0);
      cout << "\t- Drawing " << procname << endl;
      prochist->GetYaxis()->SetRangeUser(ymin*yminfactor, ymax*ymaxfactor);
      prochist->Draw((drawfirst ? "hist" : "histsame"));
      drawfirst=false;
    }
    // Draw error bands
    if (tg_systband) tg_systband->Draw("2same");
    // Re-draw ALT
    for (unsigned int ip=proc_order.size(); ip>0; ip--){
      if (!proc_order.at(ip-1).Contains("ALT")) continue;
      cout << "\t- Drawing " << proc_order.at(ip-1) << endl;
      procshape_flat_translated[proc_order.at(ip-1)].Draw("histsame");
    }
    if (tgdata_4l && tgdata_4l->GetN()>0) tgdata_4l->Draw("0psame");
    if (tgdata && tgdata->GetN()>0) tgdata->Draw("0psame");

    canvas.cd();
    legend_left.Draw();
    legend_right.Draw();
    pt.Draw();
    ptc.Draw();

    for (auto& pad:pads){
      pad->RedrawAxis();
      pad->Modified();
      pad->Update();
    }

    canvas.Modified();
    canvas.Update();
    canvas.SaveAs(TString(canvas.GetName()) + ".pdf");
    for (auto& pad:pads) pad->Close();
    canvas.Close();
    curdir->cd();

    for (auto& hh:hlist_ALT_dummy) delete hh;
    delete hdummy_ratio;
    delete tg_systband_unit;
    delete tg_systband;
    delete tgdata_unit;
    delete tgdata_4l;
    delete tgdata;

    delete systband.first;
    delete systband.second;
  }
}
