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
#include "TLine.h"
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


void replaceLastBinBoundary(TH1F* histo){
  std::vector<double> boundaries;
  for (int ix=1; ix<=histo->GetNbinsX()+1; ix++) boundaries.push_back(histo->GetBinLowEdge(ix));
  unsigned int nb = boundaries.size();
  if (boundaries.back()>=1e4) boundaries.back() = 2.*boundaries.at(nb-2) - boundaries.at(nb-3);

  TString hname = histo->GetName();
  TH1F* res = new TH1F(Form("%s_replaceLastBinBoundary", hname.Data()), histo->GetTitle(), nb-1, boundaries.data()); res->Sumw2();
  res->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
  res->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
  res->SetLineColor(histo->GetLineColor());
  res->SetMarkerColor(histo->GetMarkerColor());
  res->SetFillColor(histo->GetFillColor());
  res->SetLineStyle(histo->GetLineStyle());
  res->SetMarkerStyle(histo->GetMarkerStyle());
  res->SetFillStyle(histo->GetFillStyle());
  res->SetLineWidth(histo->GetLineWidth());
  for (unsigned int ix=1; ix<=nb-1; ix++){
    res->SetBinContent(ix, histo->GetBinContent(ix));
    res->SetBinError(ix, histo->GetBinError(ix));
  }
  *histo = *res;
  delete res;

  histo->SetName(hname);
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
    cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;

    tpl->SetTitle("");
    tpl->Scale(scale);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), 1, tpl->GetNbinsZ(), false, nullptr) << "?=" << normval << endl;

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

    tpl->SetTitle("");
    tpl->Scale(scale);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), 1, tpl->GetNbinsY(), false, nullptr) << "?=" << normval << endl;

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

    tpl->SetTitle("");
    tpl->Scale(scale);
    cout << procname << " contribution final integral = " << getHistogramIntegralAndError(tpl, 1, tpl->GetNbinsX(), false, nullptr) << "?=" << normval << endl;

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
    cout << "\t- Extracting " << strhname << " at " << nuisname << "=" << vdn << "..." << endl;
    extractTemplates(proc, data, dummymap_1D, dummymap_2D, dummymap_3D, strhname, true);

    strhname = Form("%s_%s_Up", procname.Data(), nuisname.Data());
    if (nuis_shape) nuis_shape->setVal(vup);
    if (nuis_norm) nuis_norm->setVal(vup);
    cout << "\t- Extracting " << strhname << " at " << nuisname << "=" << vup << "..." << endl;
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


TGraphAsymmErrors* getDataGraph(TH1F* hdata, bool errorsOnZero){
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
    tgdata->SetName(TString("tgdata") + (errorsOnZero ? "_withZeros" : ""));
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


void plotStacked(
  TString cinputdir,
  int onORoffshell=0, bool isEnriched=true, unsigned char markPreliminary_Supplementary=0,
  bool useLogY=false,
  bool isBlind=false, bool addObsRatio=false,
  TString strFitResultFile_SM="",
  TString strFitResultFile_GGsmALT="",
  TString strFitResultFile_NoOffshellALT="",
  TString strFitResultFile_fai1ALT="",
  TString strFitResultFile_BestFitALT="",
  TString strPostfit=""
){
  // Magic numbers
  constexpr double npixels_stdframe_xy = 800;
  constexpr double relmargin_frame_left = 0.20;
  constexpr double relmargin_frame_right = 0.05;
  constexpr double relmargin_frame_CMS = 0.07;
  constexpr double relmargin_frame_XTitle = 0.15;
  constexpr double relmargin_frame_separation = 0.1;
  constexpr double relsize_frame_ratio = 0.2;
  constexpr double relsize_frame_composition = 0.2;
  constexpr double npixels_pad_xy = 800;
  constexpr double relsize_CMSlogo = 0.98;
  constexpr double relsize_CMSlogo_sqrts = 0.8;
  constexpr double relsize_XYTitle = 0.9;
  constexpr double relsize_XYLabel = 0.8;
  constexpr double offset_xlabel = 0.004;
  constexpr double offset_ylabel = 0.007;
  constexpr double offset_xtitle = 1.09;
  constexpr double offset_ytitle = 1.3;

  unsigned int npads = 2;
  if (addObsRatio) npads++;
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
      + relsize_frame_composition
      + relmargin_frame_XTitle
      ) + 0.5
    );
  const double npixels_pad_ratio = int(
    npixels_stdframe_xy*(
      relmargin_frame_separation + relsize_frame_ratio
      ) + 0.5
    );
  const double npixels_y = npixels_pad_top + npixels_pad_bot + (npads==3 ? npixels_pad_ratio : 0.);

  constexpr double relmargin_top_legend = 0.1;
  constexpr double relmargin_bottom_legend = 0.1;
  constexpr double relmargin_left_legend = 0.05;
  constexpr double relmargin_right_legend = 0.05;
  constexpr double relsize_x_legend = 0.35;
  const double relsize_space_legend = 1. - relmargin_left_legend - relmargin_right_legend - 2.*relsize_x_legend;
  const double legend_size_mult = std::max(80.*(1.+relmargin_frame_left+relmargin_frame_right), npixels_XYTitle*1.1)/80.;
  const double npixels_perentry_y_legend = 80.*legend_size_mult; // This is for 5 entries along the vertical.
  const double npixels_x_legend = 280.*npixels_perentry_y_legend/80.;
  const double npixels_x_canvas_legend = npixels_x_legend/relsize_x_legend;

  const bool markPreliminary = markPreliminary_Supplementary==1;
  const bool markSupplementary = markPreliminary_Supplementary==2;

  gStyle->SetOptStat(0);

  TString cinputdir_lower = cinputdir; cinputdir_lower.ToLower();
  bool const is_2l2nu = cinputdir_lower.Contains("2l2nu");
  bool const is_3l1nu = cinputdir_lower.Contains("3l1nu");
  bool const is_4l = !is_2l2nu && !is_3l1nu;

  TString strChannelPlotLabel;
  if (is_4l) strChannelPlotLabel = "ZZTo4L";
  else if (is_2l2nu) strChannelPlotLabel = "ZZTo2L2Nu";
  else if (is_3l1nu) strChannelPlotLabel = "ZWTo3L1Nu";
  else{
    cerr << "\t- Cannot identify channel plot label." << endl;
    return;
  }

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

  bool const isPostFit = (strFitResultFile_SM!="");
  bool const hasBestFit = (isPostFit && strFitResultFile_BestFitALT!="");

  if (addObsRatio){
    if (!isPostFit){ cerr << "Ratio plots must include a fit result." << endl; return; }
    if (isBlind){ cerr << "Ratio plots must use unblinded results." << endl; return; }
  }

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
  unordered_map<TString, double> val_default; val_default["RV"]=1; val_default["RF"]=1; val_default["GGsm"]=1; val_default["fai1"]=0;
  unordered_map<TString, double> val_fai1ALT;
  unordered_map<TString, double> val_GGsmALT;
  unordered_map<TString, double> val_NoOffshellALT;
  unordered_map<TString, double> val_BestFitALT;
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
      if (cinputdir.Contains("/a3")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=0.05; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/a2")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=0.05; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/L1")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=0.05; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      val_GGsmALT=val_default; val_GGsmALT["GGsm"]=20./GHref;
      val_NoOffshellALT=val_default; val_NoOffshellALT["GGsm"]=0;
      if (hasBestFit) val_BestFitALT=val_default;
    }
    else if (is_3l1nu){
      if (cinputdir.Contains("/a3")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=1; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/a2")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=1; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
      else if (cinputdir.Contains("/L1")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=1; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
    }
  }
  else{ // Onshell
    if (cinputdir.Contains("/a3")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=0.5; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
    else if (cinputdir.Contains("/a2")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=-0.5; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
    else if (cinputdir.Contains("/L1")){ val_fai1ALT=val_default; val_fai1ALT["fai1"]=0.5; val_fai1ALT["RV"]=getACMuV(cinputdir, val_fai1ALT["fai1"]); val_fai1ALT["RF"]=getACMuF(cinputdir, val_fai1ALT["fai1"]); }
  }

  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_SM_map;
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_GGsmALT_map;
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_NoOffshellALT_map;
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_fai1ALT_map;
  std::unordered_map<TString, std::pair<double, std::pair<double, double>> > postfit_vals_BestFitALT_map;
  if (isPostFit && strFitResultFile_SM!=""){
    TFile* finput_postfit = TFile::Open(strFitResultFile_SM, "read");
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
        postfit_vals_SM_map[strvar] = std::pair<double, std::pair<double, double>>(val, std::pair<double, double>(val-err_dn, val+err_up));
      }
    }
    delete it;
    finput_postfit->Close();
    curdir->cd();
  }
  if (isPostFit && strFitResultFile_GGsmALT!=""){
    TFile* finput_postfit = TFile::Open(strFitResultFile_GGsmALT, "read");
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
        postfit_vals_GGsmALT_map[strvar] = std::pair<double, std::pair<double, double>>(val, std::pair<double, double>(val-err_dn, val+err_up));
      }
    }
    delete it;
    finput_postfit->Close();
    curdir->cd();
  }
  if (isPostFit && strFitResultFile_NoOffshellALT!=""){
    TFile* finput_postfit = TFile::Open(strFitResultFile_NoOffshellALT, "read");
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
        postfit_vals_NoOffshellALT_map[strvar] = std::pair<double, std::pair<double, double>>(val, std::pair<double, double>(val-err_dn, val+err_up));
      }
    }
    delete it;
    finput_postfit->Close();
    curdir->cd();
  }
  if (isPostFit && strFitResultFile_fai1ALT!=""){
    TFile* finput_postfit = TFile::Open(strFitResultFile_fai1ALT, "read");
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
        postfit_vals_fai1ALT_map[strvar] = std::pair<double, std::pair<double, double>>(val, std::pair<double, double>(val-err_dn, val+err_up));
      }
    }
    delete it;
    finput_postfit->Close();
    curdir->cd();
  }
  if (isPostFit && strFitResultFile_BestFitALT!=""){
    TFile* finput_postfit = TFile::Open(strFitResultFile_BestFitALT, "read");
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
        postfit_vals_BestFitALT_map[strvar] = std::pair<double, std::pair<double, double>>(val, std::pair<double, double>(val-err_dn, val+err_up));
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

    unordered_map<TString, std::pair<TH1F, TH1F>> syst_totalshape_map_1D;
    unordered_map<TString, std::pair<TH2F, TH2F>> syst_totalshape_map_2D;
    unordered_map<TString, std::pair<TH3F, TH3F>> syst_totalshape_map_3D;

    // Get process order/labels
    std::vector<TString> proc_order, proc_label;
    std::vector<int> proc_color;
    std::vector<int> proc_code;
    if (is_4l){
      if (onORoffshell==1){ // Offshell
        proc_order=std::vector<TString>{ "Zjets", "qqZZ", "VVZZ_offshell", "ggZZ_offshell", "total_GGsmALT", "total_NoOffshellALT" };
        if (!cinputdir.Contains("/SM")) proc_order.push_back("total_fai1ALT");
        if (hasBestFit) proc_order.push_back("total_BestFitALT");
      }
      else{
        proc_order=std::vector<TString>{ "zjets", "bkg_zz", "VVZZ", "ggH", "total_fai1ALT" };
      }
    }
    else if (is_2l2nu){
      //proc_order=std::vector<TString>{ "tZX", "qqZZ_offshell", "qqWZ_offshell", "NRB_2l2nu", "InstrMET", "VVVV_onshell", "VVVV_offshell", "ggZZ_offshell", "total_NoOffshellALT", "total_GGsmALT" };
      proc_order=std::vector<TString>{ "tZX", "qqWZ_offshell", "qqZZ_offshell", "NRB_2l2nu", "InstrMET", "VVVV_combined", "ggZZ_offshell", "total_GGsmALT", "total_NoOffshellALT" };
      if (!cinputdir.Contains("/SM")) proc_order.push_back("total_fai1ALT");
      if (hasBestFit) proc_order.push_back("total_BestFitALT");
    }
    else if (is_3l1nu){
      //proc_order=std::vector<TString>{ "tWX", "tZX", "ttbar_2l2nu", "qqZG", "DY_2l", "qqZZ", "qqWZ", "VVVV_onshell", "VVVV_offshell" };
      proc_order=std::vector<TString>{ "tVX", "ttbar_2l2nu", "qqZG", "DY_2l", "qqZZ", "qqWZ", "VVVV_combined" };
      //if (!cinputdir.Contains("/SM")) proc_order.push_back("total_fai1ALT");
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
            if (p.Contains("GGsmALT")){
              proc_label.push_back(Form("Total (%s#Gamma_{H}=%s MeV)", (aihypo=="" ? "" : "f_{ai}=0, "), getFractionString(val_GGsmALT["GGsm"]*GHref).Data()));
              proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2);
            }
            else if (p.Contains("NoOffshellALT")){
              proc_label.push_back(Form("Total (%s#Gamma_{H}=0 MeV)", (aihypo=="" ? "" : "f_{ai}=0, ")));
              proc_color.push_back(int(kOrange-3)); proc_code.push_back(-2);
            }
            else if (p.Contains("BestFitALT")){
              proc_label.push_back("Total best fit");
              proc_color.push_back(int(kAzure-6)); proc_code.push_back(-2);
            }
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
        if (p=="qqZZ_offshell"){ proc_color.push_back(int(TColor::GetColor("#99ccff"))); proc_label.push_back("q#bar{q}#rightarrowZZ"); proc_code.push_back(0); }
        else if (p=="qqWZ_offshell"){ proc_color.push_back(int(TColor::GetColor("#00ff00"))); proc_label.push_back("q#bar{q}'#rightarrowWZ"); proc_code.push_back(0); }
        else if (p=="InstrMET"){ proc_color.push_back(int(TColor::GetColor("#669966"))); proc_label.push_back("Instr. p_{T}^{miss}"); proc_code.push_back(0); }
        else if (p=="NRB_2l2nu"){ proc_color.push_back(int(kGray+1)); proc_label.push_back("Nonresonant"); proc_code.push_back(0); }
        else if (p=="tZX"){ proc_color.push_back(int(kOrange-6)); proc_label.push_back("tZ+X"); proc_code.push_back(0); }
        else if (p=="ggZZ_offshell"){
          proc_color.push_back(int(TColor::GetColor("#ffdcdc")));
          proc_label.push_back("gg SM total"); proc_code.push_back(2);
        }
        else if (p=="VVVV_combined"){
          proc_color.push_back(int(TColor::GetColor("#ff9b9b")));
          proc_label.push_back("EW SM total"); proc_code.push_back(1);
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
            if (p.Contains("GGsmALT")){
              proc_label.push_back(Form("Total (%s#Gamma_{H}=%s MeV)", (aihypo=="" ? "" : "f_{ai}=0, "), getFractionString(val_GGsmALT["GGsm"]*GHref).Data()));
              proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2);
            }
            else if (p.Contains("NoOffshellALT")){
              proc_label.push_back(Form("Total (%s#Gamma_{H}=0 MeV)", (aihypo=="" ? "" : "f_{ai}=0, ")));
              proc_color.push_back(int(kOrange-3)); proc_code.push_back(-2);
            }
            else if (p.Contains("BestFitALT")){
              proc_label.push_back("Total best fit");
              proc_color.push_back(int(kAzure-6)); proc_code.push_back(-2);
            }
            else if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total (#bar{f}_{ai}=%s, #Gamma_{H}=#Gamma_{H}^{SM})", getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-2); }
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
        if (p=="qqZZ"){ proc_color.push_back(int(TColor::GetColor("#99ccff"))); proc_label.push_back("q#bar{q}#rightarrowZZ"); proc_code.push_back(0); }
        else if (p=="qqWZ"){ proc_color.push_back(int(TColor::GetColor("#00ff00"))); proc_label.push_back("q#bar{q}'#rightarrowWZ"); proc_code.push_back(0); }
        else if (p=="qqZG"){ proc_color.push_back(int(TColor::GetColor("#f1c232"))); proc_label.push_back("q#bar{q}#rightarrowZ#gamma"); proc_code.push_back(0); }
        else if (p=="DY_2l"){ proc_color.push_back(int(TColor::GetColor("#669966"))); proc_label.push_back("Drell-Yan"); proc_code.push_back(0); }
        else if (p=="ttbar_2l2nu"){ proc_color.push_back(int(kGray+1)); proc_label.push_back("t#bar{t}"); proc_code.push_back(0); }
        else if (p=="tZX"){ proc_color.push_back(int(kOrange-6)); proc_label.push_back("tZ+X"); proc_code.push_back(0); }
        else if (p=="tWX"){ proc_color.push_back(int(TColor::GetColor("#674ea7"))); proc_label.push_back("tW+X"); proc_code.push_back(0); }
        else if (p=="tVX"){ proc_color.push_back(int(kOrange-6)); proc_label.push_back("tV+X"); proc_code.push_back(0); }
        else if (p=="VVVV_combined"){
          proc_color.push_back(int(TColor::GetColor("#ff9b9b")));
          proc_label.push_back("EW SM total"); proc_code.push_back(1);
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
            if (p.Contains("GGsmALT")){
              proc_label.push_back(Form("Total (%s#Gamma_{H}=%s GeV)", (aihypo=="" ? "" : "f_{ai}=0, "), getFractionString(val_GGsmALT["GGsm"]*GHref/1000.).Data()));
              proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2);
            }
            else if (p.Contains("NoOffshellALT")){
              proc_label.push_back(Form("Total (%s#Gamma_{H}=0 MeV)", (aihypo=="" ? "" : "f_{ai}=0, ")));
              proc_color.push_back(int(kOrange-3)); proc_code.push_back(-2);
            }
            else if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total (#bar{f}_{ai}=%s, #Gamma_{H}=#Gamma_{H}^{SM})", getFractionString(val_fai1ALT["fai1"]).Data())); proc_color.push_back(int(kViolet)); proc_code.push_back(-2); }
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
        TString cinput_main = cinputdir + "/" + sqrtsname + "/hto" + channame + "_" + catname;
        TString cinput_file = cinput_main + ".input.root";
        TString cinput_dc = cinput_main + ".txt";
        if (!HostHelpers::FileExists(cinput_file)){ cout << "File " << cinput_file << " does not exist. Skipping this channel..." << endl; continue; }
        else if (!HostHelpers::FileExists(cinput_dc)){ cout << "File " << cinput_dc << " does not exist. Skipping this channel..." << endl; continue; }

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
          std::vector<double>::iterator it_rate=procrate.begin();
          for (TString const& pname:procname){
            cout << "Process " << pname << " is found in the datacard" << endl;
            procSpecs[pname] = process_spec(pname, *it_rate);
            auto const& rateModifiers = procname_rateModifier_map.find(pname)->second;
            procSpecs[pname].setRateModifiers(rateModifiers);
            it_rate++;
          }
        }

        cout << "Opening the workspace file " << cinput_file << "..." << endl;

        // Open the input
        TFile* finput = TFile::Open(cinput_file, "read");
        RooWorkspace* ws = (RooWorkspace*) finput->Get("w");

        std::vector<RooRealVar*> nuisanceVars;

        // Set postfit values first if they are present.
        for (auto const& pp:postfit_vals_SM_map){
          RooRealVar* tmpvar = dynamic_cast<RooRealVar*>(ws->var(pp.first));
          double const& vnom = pp.second.first;
          double const& vlow = pp.second.second.first;
          double const& vhigh = pp.second.second.second;
          if (tmpvar){
            nuisanceVars.push_back(tmpvar);
            cout << "\t- Setting " << pp.first << " = " << vnom << " [" << vlow << ", " << vhigh << "] in the workspace." << endl;
            tmpvar->setVal(vnom);
            tmpvar->setAsymError(vlow-vnom, vhigh-vnom);
          }

          for (auto const& lnNmodvar:lnNmodvars){
            if (TString(lnNmodvar->GetName()) == pp.first){
              cout << "\t- Setting " << pp.first << " = " << vnom << " [" << vlow << ", " << vhigh << "] in the list of lnN nuisances." << endl;
              lnNmodvar->setVal(vnom);
              lnNmodvar->setAsymError(vlow-vnom, vhigh-vnom);
            }
          }
        }
        cout << "Building the zipped map..." << endl;
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

        // Reset all parameters to SM fitted values
        for (auto const& pp:postfit_vals_SM_map){
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
            }
            else if (pname.Contains("ttH") || pname.Contains("bbH")){
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "ggH");
            }
            else if (pname.Contains("VBF")){
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVZZ");
              controlVars["kbkg_VBF"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_vv");
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);
            }
            else if (pname.Contains("ZH") || pname.Contains("WH")){
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVZZ");
              controlVars["kbkg_VBF"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_vv");
              ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);
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
            if (pname.Contains("VVVV_offshell")) extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVVV_combined");
          }
          else if (onORoffshell && (pname.Contains("ggZZ_onshell") || pname.Contains("VVZZ_onshell") || pname.Contains("VVVV_onshell"))){ // On-shell process in off-shell data cards
            if (pname.Contains("VVVV_onshell")) extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "VVVV_combined");
          }

          if (is_3l1nu){
            if (pname.Contains("tZX") || pname.Contains("tWX")) extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "tVX");
          }

          if (!val_GGsmALT.empty()){
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_GGsmALT["fai1"]);
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
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_default["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_default["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_default["RV"]);
            setControlVariableValue(controlVars, "RF", val_default["RF"]);
            for (auto const& pp:postfit_vals_SM_map){
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

          if (!val_NoOffshellALT.empty()){
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_NoOffshellALT["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_NoOffshellALT["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_NoOffshellALT["RV"]);
            setControlVariableValue(controlVars, "RF", val_NoOffshellALT["RF"]);
            for (auto const& pp:postfit_vals_NoOffshellALT_map){
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
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_NoOffshellALT");
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_default["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_default["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_default["RV"]);
            setControlVariableValue(controlVars, "RF", val_default["RF"]);
            for (auto const& pp:postfit_vals_SM_map){
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

          if (!val_BestFitALT.empty()){
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_BestFitALT["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_BestFitALT["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_BestFitALT["RV"]);
            setControlVariableValue(controlVars, "RF", val_BestFitALT["RF"]);
            for (auto const& pp:postfit_vals_BestFitALT_map){
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
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_BestFitALT");
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_default["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_default["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_default["RV"]);
            setControlVariableValue(controlVars, "RF", val_default["RF"]);
            for (auto const& pp:postfit_vals_SM_map){
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

          if (!val_fai1ALT.empty()){
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_fai1ALT["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_fai1ALT["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_fai1ALT["RV"]);
            setControlVariableValue(controlVars, "RF", val_fai1ALT["RF"]);
            for (auto const& pp:postfit_vals_fai1ALT_map){
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
            ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D, "total_fai1ALT");
            setControlVariableValue(controlVars, "CMS_zz4l_fai1", val_default["fai1"]);
            setControlVariableValue(controlVars, "GGsm", val_default["GGsm"]);
            setControlVariableValue(controlVars, "RV", val_default["RV"]);
            setControlVariableValue(controlVars, "RF", val_default["RF"]);
            for (auto const& pp:postfit_vals_SM_map){
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

          ndims = extractTemplates(procSpecs[pname], data, procshape_1D, procshape_2D, procshape_3D);
          if (addObsRatio) extractTemplateSystVariations(procSpecs[pname], data, zipped_nuisances, postfit_vals_SM_map);
        }
        extractDataTemplates(procSpecs[procname.front()], data, procshape_1D, procshape_2D, procshape_3D, "data");

        if (addObsRatio){
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
    std::vector<std::pair<TH1F*, TH1F*>> totalsysts(ndims, std::pair<TH1F*, TH1F*>(nullptr, nullptr)); // Dn/Up variations of each dimension
    for (unsigned int idim=0; idim<ndims; idim++){
      TString dimname = Form("_KD%i", idim+1);
      cout << "Preparing projections on dimension " << idim << "..." << endl;
      bool needSysts = true;
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

            const float valKDCut=(onORoffshell==1 ? (varlabels.at(kdDim).Contains("p_{T}^{miss}") ? 200. : (is_2l2nu ? 0.8 : 0.6)) : 0.5);
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
          replaceLastBinBoundary(htmp);
          htmp->GetXaxis()->SetTitle(varlabels.at(idim));
          cout << "\t- Constructed histogram " << htmp->GetName() << endl;
          cout << "\t\t- Initial histogram integral: " << getHistogramIntegralAndError(hist, 1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1, hist->GetNbinsZ(), false, nullptr) << endl;
          if (procdist.find(it->first)==procdist.end()) procdist[it->first] = std::vector<TH1F*>();
          procdist[it->first].push_back(htmp);
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first].back(), 1, procdist[it->first].back()->GetNbinsX(), false, nullptr) << endl;

          if (addObsRatio && needSysts){
            auto& totalsysthists = totalsysts.at(idim);
            for (auto syst_totalshape:syst_totalshape_map_3D){
              TH1F* htmp_dn=getHistogramSlice(
                &(syst_totalshape.second.first), idim,
                iy, jy, iz, jz,
                syst_totalshape.first + "_" + dimname + "_hist_dn"
              );
              replaceLastBinBoundary(htmp_dn);
              htmp_dn->GetXaxis()->SetTitle(varlabels.at(idim));
              TH1F* htmp_up=getHistogramSlice(
                &(syst_totalshape.second.second), idim,
                iy, jy, iz, jz,
                syst_totalshape.first + "_" + dimname + "_hist_up"
              );
              replaceLastBinBoundary(htmp_up);
              htmp_up->GetXaxis()->SetTitle(varlabels.at(idim));

              if (!totalsysthists.first){
                totalsysthists.first = (TH1F*) htmp_dn->Clone(Form("allsysts_%s_dn", dimname.Data())); totalsysthists.first->Reset("ICESM");
                totalsysthists.second = (TH1F*) htmp_up->Clone(Form("allsysts_%s_up", dimname.Data())); totalsysthists.second->Reset("ICESM");
              }

              for (int ix=1; ix<=htmp->GetNbinsX(); ix++){
                double bdn = htmp_dn->GetBinContent(ix);
                double bup = htmp_up->GetBinContent(ix);
                if (bdn>bup) std::swap(bdn, bup);
                totalsysthists.first->SetBinContent(ix, -std::sqrt(std::pow(totalsysthists.first->GetBinContent(ix), 2) + std::pow(std::min(0., bdn), 2)));
                totalsysthists.second->SetBinContent(ix, std::sqrt(std::pow(totalsysthists.second->GetBinContent(ix), 2) + std::pow(std::max(0., bup), 2)));
              }

              delete htmp_dn;
              delete htmp_up;
            }
            needSysts = false;
          }
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
          replaceLastBinBoundary(htmp);
          htmp->GetXaxis()->SetTitle(varlabels.at(idim));
          cout << "\t- Constructed histogram " << htmp->GetName() << endl;
          cout << "\t\t- Initial histogram integral: " << getHistogramIntegralAndError(hist, 1, hist->GetNbinsX(), 1, hist->GetNbinsY(), false, nullptr) << endl;
          if (procdist.find(it->first)==procdist.end()) procdist[it->first] = std::vector<TH1F*>();
          procdist[it->first].push_back(htmp);
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first].back(), 1, procdist[it->first].back()->GetNbinsX(), false, nullptr) << endl;

          if (addObsRatio && needSysts){
            auto& totalsysthists = totalsysts.at(idim);
            for (auto syst_totalshape:syst_totalshape_map_2D){
              TH1F* htmp_dn=getHistogramSlice(
                &(syst_totalshape.second.first), idim,
                iy, jy,
                syst_totalshape.first + "_" + dimname + "_hist_dn"
              );
              replaceLastBinBoundary(htmp_dn);
              htmp_dn->GetXaxis()->SetTitle(varlabels.at(idim));
              TH1F* htmp_up=getHistogramSlice(
                &(syst_totalshape.second.second), idim,
                iy, jy,
                syst_totalshape.first + "_" + dimname + "_hist_up"
              );
              replaceLastBinBoundary(htmp_up);
              htmp_up->GetXaxis()->SetTitle(varlabels.at(idim));

              if (!totalsysthists.first){
                totalsysthists.first = (TH1F*) htmp_dn->Clone(Form("allsysts_%s_dn", dimname.Data())); totalsysthists.first->Reset("ICESM");
                totalsysthists.second = (TH1F*) htmp_up->Clone(Form("allsysts_%s_up", dimname.Data())); totalsysthists.second->Reset("ICESM");
              }

              for (int ix=1; ix<=htmp->GetNbinsX(); ix++){
                double bdn = htmp_dn->GetBinContent(ix);
                double bup = htmp_up->GetBinContent(ix);
                if (bdn>bup) std::swap(bdn, bup);
                totalsysthists.first->SetBinContent(ix, -std::sqrt(std::pow(totalsysthists.first->GetBinContent(ix), 2) + std::pow(std::min(0., bdn), 2)));
                totalsysthists.second->SetBinContent(ix, std::sqrt(std::pow(totalsysthists.second->GetBinContent(ix), 2) + std::pow(std::max(0., bup), 2)));
              }

              delete htmp_dn;
              delete htmp_up;
            }
            needSysts = false;
          }
        }
      }
      else{
        for (auto it=procshape_1D.begin(); it!=procshape_1D.end(); it++){
          TH1F* hist = &(it->second);
          TH1F* htmp = getHistogramSlice(
            hist, it->first + dimname + "_hist"
          );
          replaceLastBinBoundary(htmp);
          htmp->GetXaxis()->SetTitle(varlabels.at(idim));
          cout << "\t- Constructed histogram " << htmp->GetName() << endl;
          cout << "\t\t- Initial histogram integral: " << getHistogramIntegralAndError(hist, 1, hist->GetNbinsX(), false, nullptr) << endl;
          if (procdist.find(it->first)==procdist.end()) procdist[it->first] = std::vector<TH1F*>();
          procdist[it->first].push_back(htmp);
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first].back(), 1, procdist[it->first].back()->GetNbinsX(), false, nullptr) << endl;

          if (addObsRatio && needSysts){
            auto& totalsysthists = totalsysts.at(idim);
            for (auto syst_totalshape:syst_totalshape_map_1D){
              TH1F* htmp_dn=getHistogramSlice(
                &(syst_totalshape.second.first),
                syst_totalshape.first + "_" + dimname + "_hist_dn"
              );
              replaceLastBinBoundary(htmp_dn);
              htmp_dn->GetXaxis()->SetTitle(varlabels.at(idim));
              TH1F* htmp_up=getHistogramSlice(
                &(syst_totalshape.second.second),
                syst_totalshape.first + "_" + dimname + "_hist_up"
              );
              replaceLastBinBoundary(htmp_up);
              htmp_up->GetXaxis()->SetTitle(varlabels.at(idim));

              if (!totalsysthists.first){
                totalsysthists.first = (TH1F*) htmp_dn->Clone(Form("allsysts_%s_dn", dimname.Data())); totalsysthists.first->Reset("ICESM");
                totalsysthists.second = (TH1F*) htmp_up->Clone(Form("allsysts_%s_up", dimname.Data())); totalsysthists.second->Reset("ICESM");
              }

              for (int ix=1; ix<=htmp->GetNbinsX(); ix++){
                double bdn = htmp_dn->GetBinContent(ix);
                double bup = htmp_up->GetBinContent(ix);
                if (bdn>bup) std::swap(bdn, bup);
                totalsysthists.first->SetBinContent(ix, -std::sqrt(std::pow(totalsysthists.first->GetBinContent(ix), 2) + std::pow(std::min(0., bdn), 2)));
                totalsysthists.second->SetBinContent(ix, std::sqrt(std::pow(totalsysthists.second->GetBinContent(ix), 2) + std::pow(std::max(0., bup), 2)));
              }

              delete htmp_dn;
              delete htmp_up;
            }
            needSysts = false;
          }
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
        if (proccode!=0 || (is_3l1nu && (procname=="qqWZ" || procname=="qqZZ"))) nleft++;
      }
      unsigned int nright = proc_order.size()-nleft;

      TString strLegendCanvasNameCore = TString((onORoffshell ? "Offshell_" : "Onshell_")) + strChannelPlotLabel + "_" + (aihypo=="" ? "SM" : aihypo.Data()) + (markPreliminary ? "_Preliminary" : (markSupplementary ? "_Supplementary" : ""));

      const unsigned int nmax_legend = std::max(nleft, nright);
      const double npixels_y_canvas_legend = npixels_perentry_y_legend*nmax_legend/(1. - relmargin_top_legend - relmargin_bottom_legend);
      TCanvas canvas(
        TString("c_") + strLegendCanvasNameCore + "_legend",
        "", 8, 30, npixels_x_canvas_legend, npixels_y_canvas_legend
      );
      cout << "Legend pixel dimensions: " << npixels_x_canvas_legend << " x " << npixels_y_canvas_legend << endl;
      canvas.cd();
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(relmargin_left_legend);
      canvas.SetRightMargin(relmargin_right_legend);
      canvas.SetTopMargin(relmargin_top_legend);
      canvas.SetBottomMargin(relmargin_bottom_legend);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);

      double leg_xmin=relmargin_left_legend + relsize_space_legend + relsize_x_legend;
      double leg_xmax=1. - relmargin_right_legend;
      double leg_ymax=1. - relmargin_top_legend;
      double leg_ymin=relmargin_bottom_legend; if (nright<nleft) leg_ymin += (1. - relmargin_top_legend - relmargin_bottom_legend)/nmax_legend*(nleft-nright);
      TLegend legend_right(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
      legend_right.SetBorderSize(0);
      legend_right.SetTextFont(43);
      legend_right.SetTextSize(npixels_XYTitle);
      legend_right.SetLineColor(1);
      legend_right.SetLineStyle(1);
      legend_right.SetLineWidth(1);
      legend_right.SetFillColor(0);
      legend_right.SetFillStyle(0);

      leg_xmin=relmargin_left_legend;
      leg_xmax=relmargin_left_legend + relsize_x_legend;
      leg_ymin=relmargin_bottom_legend; if (nright>nleft) leg_ymin += (1. - relmargin_top_legend - relmargin_bottom_legend)/nmax_legend*(nright-nleft);
      TLegend legend_left(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
      legend_left.SetBorderSize(0);
      legend_left.SetTextFont(43);
      legend_left.SetTextSize(npixels_XYTitle);
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

        prochist->SetLineWidth(1);
        if (procname.Contains("ALT")){
          if (procname.Contains("GGsm")) prochist->SetLineStyle(2);
          else if (procname.Contains("NoOffshell")) prochist->SetLineStyle(4);
          else if (procname.Contains("fai1")) prochist->SetLineStyle(7);
          else if (procname.Contains("BestFit")) prochist->SetLineStyle(9);
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColor(proc_color[ip]);
          prochist->SetLineWidth(2);
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

        if (procname=="data") tgdata=getDataGraph(prochist, false);
      }

      for (unsigned int ip=proc_order.size(); ip>0; ip--){
        TString const& procname = proc_order.at(ip-1);
        TString const& proclabel = proc_label.at(ip-1);
        auto const& proccode = proc_code.at(ip-1);
        cout << "Adding process " << procname << " to legend..." << endl;
        TH1F*& prochist = procdist[procname].front();
        if (!prochist) cout << procname << " histogram is null!" << endl;

        TLegend* legend_chosen = nullptr;
        if (proccode!=0 || (is_3l1nu && (procname=="qqWZ" || procname=="qqZZ"))){
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

      auto& systband = totalsysts.at(idim);
      TH1F* hsum_nominal = nullptr;
      std::vector<TH1F*> hlist_ALT;
      for (unsigned int ip=0; ip<proc_order.size(); ip++){
        if (ip>0){
          if (!proc_order.at(ip).Contains("ALT") && proc_order.at(ip)!="data"){
            procdist[proc_order.at(ip)].at(idim)->Add(procdist[proc_order.at(ip-1)].at(idim));
            hsum_nominal = procdist[proc_order.at(ip)].at(idim);
          }
          else if (proc_order.at(ip)!="data"){
            for (unsigned int jp=ip-1; jp>0; jp--){
              if (proc_code[jp]==0){
                //procdist[proc_order.at(ip)].at(idim)->Add(procdist[proc_order.at(jp)].at(idim)); // Add the first bkg process and break
                break;
              }
            }
            hlist_ALT.push_back(procdist[proc_order.at(ip)].at(idim));
          }
        }
        cout << proc_order.at(ip) << " integral: " << getHistogramIntegralAndError(procdist[proc_order.at(ip)].at(idim), 1, procdist[proc_order.at(ip)].at(idim)->GetNbinsX(), false, nullptr) << endl;
      }

      TGraphAsymmErrors* tg_systband = nullptr;
      if (addObsRatio){
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

      // Divide by bin width when plotting off-shell mass
      //if (onORoffshell==1 && (int) idim==massDim){ for (auto proc:proc_order) divideBinWidth(procdist[proc]); }

      // Draw
      cout << "Creating the canvas..." << endl;

      TString strCanvasNameCore = TString((onORoffshell ? "Offshell_" : "Onshell_")) + strChannelPlotLabel + "_" + (aihypo=="" ? "SM" : aihypo.Data()) + (isEnriched ? "_SignalEnriched" : "") + (markPreliminary ? "_Preliminary" : (markSupplementary ? "_Supplementary" : ""));

      TCanvas canvas(
        TString("c_") + strCanvasNameCore + "_" + catname + dimname + (isPostFit ? "_postfit" : ""),
        "", 8, 30, npixels_x, npixels_y
      );
      canvas.cd();

      cout << "\t- Creating the pads..." << endl;
      TPad* pad_main = nullptr;
      TPad* pad_ratio = nullptr;
      TPad* pad_comp = nullptr;
      std::vector<TPad*> pads;
      canvas.cd();
      pads.push_back(
        new TPad(
          Form("%s_top", canvas.GetName()), "",
          0, (1.-npixels_pad_top/npixels_y), 1, 1
        )
      );
      pad_main = pads.back();
      if (addObsRatio){
        canvas.cd();
        pads.push_back(
          new TPad(
            Form("%s_rat", canvas.GetName()), "",
            0, npixels_pad_bot/npixels_y, 1, (1.-npixels_pad_top/npixels_y)
          )
        );
        pad_ratio = pads.back();
      }
      canvas.cd();
      pads.push_back(
        new TPad(
          Form("%s_comp", canvas.GetName()), "",
          0, 0, 1, npixels_pad_bot/npixels_y
        )
      );
      pad_comp = pads.back();
      {
        unsigned int ipad=0;
        for (auto& pad:pads){
          pad->cd();
          pad->SetFillColor(0);
          pad->SetBorderMode(0);
          pad->SetBorderSize(2);
          pad->SetFrameFillStyle(0);
          pad->SetFrameBorderMode(0);
          pad->SetTickx(1);
          pad->SetTicky(1);
          pad->SetLeftMargin(relmargin_frame_left/(1.+relmargin_frame_left+relmargin_frame_right));
          pad->SetRightMargin(relmargin_frame_right/(1.+relmargin_frame_left+relmargin_frame_right));
          if (pad==pad_main){
            pad->SetTopMargin(npixels_stdframe_xy*relmargin_frame_CMS/npixels_pad_top);
            pad->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_separation/2./npixels_pad_top);
          }
          else if (pad==pad_comp){
            pad->SetTopMargin(npixels_stdframe_xy*relmargin_frame_separation/2./npixels_pad_bot);
            pad->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_XTitle/npixels_pad_bot);
          }
          else{
            pad->SetTopMargin(npixels_stdframe_xy*relmargin_frame_separation/2./npixels_pad_ratio);
            pad->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_separation/2./npixels_pad_ratio);
          }
          if (varlabels.at(idim).Contains("m_{T}^{ZZ}") || varlabels.at(idim).Contains("m_{T}^{WZ}") || varlabels.at(idim).Contains("p_{T}^{miss}")){
            pad->SetLogx(true);
            if (useLogY && pad==pad_main) pad->SetLogy(true);
          }
          canvas.cd();
          ipad++;
        }
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
        text = pt.AddText(npixels_CMSlogo*2.2/npixels_pad_xy, 0.45, "Preliminary");
        text->SetTextFont(53);
        text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
        text->SetTextAlign(12);
      }
      else if (markSupplementary){
        text = pt.AddText(npixels_CMSlogo*2.2/npixels_pad_xy, 0.45, "Supplementary");
        text->SetTextFont(53);
        text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
        text->SetTextAlign(12);
      }
      int theSqrts=13;
      TString cErgTev = Form("138 fb^{-1} (%i TeV)", theSqrts);
      text = pt.AddText(0.999, 0.45, cErgTev);
      text->SetTextFont(43);
      text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
      text->SetTextAlign(32);

      TString strCatLabel;
      if (catname=="Untagged") strCatLabel="Untagged";
      else if (catname=="JJVBFTagged") strCatLabel="VBF-tagged";
      else if (catname=="HadVHTagged") strCatLabel="VH-tagged";
      else if (catname=="Nj_eq_0") strCatLabel="N_{j}=0";
      else if (catname=="Nj_eq_1") strCatLabel="N_{j}=1";
      else if (catname.Contains("Nj_geq_2")) strCatLabel="N_{j}#geq2";
      else if (catname=="BoostedHadVH") strCatLabel="Boosted V #rightarrow J";
      if (ailabel!="") strCatLabel = strCatLabel + ", " + ailabel + " analysis";
      if (isPostFit){
        if (strPostfit=="") strCatLabel = strCatLabel + " (postfit)";
        else strCatLabel = strCatLabel + Form(" (%s)", strPostfit.Data());
      }
      else if (strPostfit!="") strCatLabel = strCatLabel + Form(" (%s)", strPostfit.Data());

      TPaveText pt_cat(
        npixels_stdframe_xy*relmargin_frame_left/npixels_x,
        (npixels_y - npixels_stdframe_xy*relmargin_frame_CMS-1 - 0.03*npixels_stdframe_xy - npixels_XYTitle*1.25)/npixels_y,
        1.-npixels_stdframe_xy*relmargin_frame_right/npixels_x,
        (npixels_y - npixels_stdframe_xy*relmargin_frame_CMS-1 - 0.03*npixels_stdframe_xy)/npixels_y,
        "brNDC"
      );
      pt_cat.SetBorderSize(0);
      pt_cat.SetFillStyle(0);
      pt_cat.SetTextAlign(12);
      pt_cat.SetTextFont(43);
      pt_cat.SetTextSize(npixels_XYTitle);
      text = pt_cat.AddText(0.05, 0.45, strCatLabel);

      cout << "\t- Preparing the cut label..." << endl;
      TString strCutLabel;
      if ((int) idim!=massDim && massDim>=0 && cutvals.at(idim).at(massDim)>0.){
        TString strUnit = "GeV";
        if (is_2l2nu) strCutLabel += "m_{T}^{ZZ}>";
        else if (is_3l1nu) strCutLabel += "m_{T}^{WZ}>";
        else strCutLabel += "m_{4l}>";
        strCutLabel += Form("%.0f", cutvals.at(idim).at(massDim));
        if (strUnit!="") strCutLabel = strCutLabel + " " + strUnit;
      }
      if (is_2l2nu && catname.Contains("Nj_geq_2")){
        if (strCutLabel!="") strCutLabel += ", ";
        if (catname.Contains("pTmiss_lt_200")) strCutLabel += "p_{T}^{miss}<200 GeV";
        else if (catname.Contains("pTmiss_ge_200")) strCutLabel += "p_{T}^{miss}>200 GeV";
      }
      if ((int) idim!=kdDim && kdDim>=0 && cutvals.at(idim).at(kdDim)>0.){
        TString strUnit = "";
        if (strCutLabel!="") strCutLabel += ", ";
        if (is_2l2nu){
          strCutLabel += (catname.Contains("Nj_geq_2") ? "D_{2jet}^{VBF}>" : "p_{T}^{miss}>");
          if (!catname.Contains("Nj_geq_2")) strUnit = "GeV";
        }
        else strCutLabel += "D_{bkg}>";
        if (cutvals.at(idim).at(kdDim)>1.) strCutLabel += Form("%.0f", cutvals.at(idim).at(kdDim));
        else strCutLabel += Form("%.1f", cutvals.at(idim).at(kdDim));
        if (strUnit!="") strCutLabel = strCutLabel + " " + strUnit;
      }
      TPaveText pt_cut(
        npixels_stdframe_xy*relmargin_frame_left/npixels_x,
        (npixels_y - npixels_stdframe_xy*relmargin_frame_CMS-1 - 0.03*npixels_stdframe_xy - npixels_XYTitle*2.*1.25)/npixels_y,
        1.-npixels_stdframe_xy*relmargin_frame_right/npixels_x,
        (npixels_y - npixels_stdframe_xy*relmargin_frame_CMS-1 - 0.03*npixels_stdframe_xy - npixels_XYTitle*1.25)/npixels_y,
        "brNDC"
      );
      pt_cut.SetBorderSize(0);
      pt_cut.SetFillStyle(0);
      pt_cut.SetTextAlign(12);
      pt_cut.SetTextFont(43);
      pt_cut.SetTextSize(npixels_XYTitle);
      text = pt_cut.AddText(0.05, 0.45, strCutLabel);

      float ymax=-1, ymin=-1;
      if (useLogY) ymin=9e9;
      float xmin=-1, xmax=-1;
      TGraphAsymmErrors* tgdata=nullptr;
      TGraphAsymmErrors* tgdata_withZeros=nullptr;
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
        prochist->SetLineWidth(1);
        if (procname.Contains("ALT")){
          if (procname.Contains("GGsm")) prochist->SetLineStyle(2);
          else if (procname.Contains("NoOffshell")) prochist->SetLineStyle(4);
          else if (procname.Contains("fai1")) prochist->SetLineStyle(7);
          else if (procname.Contains("BestFit")) prochist->SetLineStyle(9);
          prochist->SetMarkerColor(proc_color[ip]);
          prochist->SetLineColor(proc_color[ip]);
          prochist->SetLineWidth(2);
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
        prochist->GetXaxis()->SetLabelFont(43);
        prochist->GetXaxis()->SetLabelOffset(offset_xlabel);
        prochist->GetXaxis()->SetLabelSize(npixels_XYLabel);
        prochist->GetXaxis()->SetTitleFont(42);
        prochist->GetXaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_top);
        prochist->GetXaxis()->SetTitleOffset(offset_xtitle);
        prochist->GetYaxis()->SetNdivisions(505);
        prochist->GetYaxis()->SetLabelFont(43);
        prochist->GetYaxis()->SetLabelOffset(offset_ylabel);
        prochist->GetYaxis()->SetLabelSize(npixels_XYLabel);
        prochist->GetYaxis()->SetTitleFont(42);
        prochist->GetYaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_top);
        prochist->GetYaxis()->SetTitleOffset(offset_ytitle);
        prochist->GetYaxis()->SetTitle("Events / bin");
        prochist->GetYaxis()->CenterTitle();
        prochist->GetXaxis()->CenterTitle();

        if (procname!="data"){
          for (int ix=binXlow; ix<=binXhigh; ix++){
            float bc = prochist->GetBinContent(ix);
            if (bc!=0.){
              ymax = std::max(bc, ymax);
              ymin = std::min(bc, ymin);
            }
          }
        }
        else{
          cout << "\t- Obtaining " << procname << " graph" << endl;
          tgdata=getDataGraph(prochist, false);
          tgdata_withZeros = getDataGraph(prochist, true);
          if (tgdata){
            cout << "\t\t- Np = " << tgdata->GetN() << endl;
            for (int ipoint=0; ipoint<tgdata->GetN(); ipoint++){
              float bc = tgdata->GetY()[ipoint]+tgdata->GetEYhigh()[ipoint];
              float bc_low = tgdata->GetY()[ipoint]-std::abs(tgdata->GetEYlow()[ipoint]);
              if (bc!=0.){
                ymax = std::max(bc, ymax);
              }
              if (bc_low!=0.){
                ymin = std::min(bc_low, ymin);
              }
            }
            cout << "\t\t- Success!" << endl;
          }
          else cout << "-t-t- Failure!" << endl;
        }
      }
      if (tg_systband){
        for (int ipoint=0; ipoint<tg_systband->GetN(); ipoint++){
          float bc = tg_systband->GetY()[ipoint]+tg_systband->GetEYhigh()[ipoint];
          float bc_low = tg_systband->GetY()[ipoint]-std::abs(tg_systband->GetEYlow()[ipoint]);
          ymax = std::max(bc, ymax);
          ymin = std::min(bc_low, ymin);
        }
      }

      if (!useLogY) ymin=0;

      float ymaxfactor = (!useLogY ? 1.25 : 250.);
      float yminfactor = 1;

      vector<TH1F*> intermediateHistList;
      canvas.cd();
      pad_main->cd();
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
      // Draw error bands
      if (tg_systband) tg_systband->Draw("2same");
      // Re-draw ALT
      for (unsigned int ip=proc_order.size(); ip>0; ip--){
        if (!proc_order.at(ip-1).Contains("ALT")) continue;
        cout << "\t- Drawing " << proc_order.at(ip-1) << endl;
        procdist[proc_order.at(ip-1)].at(idim)->Draw("histsame");
      }
      if (tgdata && tgdata->GetN()>0 && !isBlind) tgdata->Draw("e1psame");

      canvas.cd();
      pt.Draw();
      pt_cat.Draw();
      if (strCutLabel!="") pt_cut.Draw();

      pad_comp->cd();
      drawfirst=false;
      for (int ix=1; ix<=intermediateHistList.front()->GetNbinsX(); ix++){
        double bc_sum = intermediateHistList.front()->GetBinContent(ix);
        for (auto& hh:intermediateHistList) hh->SetBinContent(ix, hh->GetBinContent(ix)/bc_sum);
      }
      for (auto& hh:intermediateHistList){
        hh->GetYaxis()->SetRangeUser(0, 1);
        hh->GetYaxis()->SetNdivisions(504);
        hh->GetYaxis()->SetTitle("Comp.");
        //hh->GetYaxis()->SetTitleOffset(1.1*(2./11. + 0.13 * 8./11. + 0.5/11.) / (1.-(2./11. + 0.13 * 8./11. + 0.5/11.)));

        hh->GetXaxis()->SetTitleFont(42);
        hh->GetXaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_bot);
        hh->GetXaxis()->SetTitleOffset(offset_xtitle);
        hh->GetYaxis()->SetTitleFont(42);
        hh->GetYaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_bot);
        hh->GetYaxis()->SetTitleOffset(offset_ytitle*(npixels_pad_bot/npixels_pad_top));

        hh->Draw((drawfirst ? "hist" : "histsame"));
        drawfirst=false;
      }


      canvas.cd();

      TGraphAsymmErrors* tg_systband_unit = nullptr;
      TGraphAsymmErrors* tgdata_withZeros_unit = nullptr;
      TH1F* hdummy_ratio = nullptr;
      std::vector<TH1F*> hlist_ALT_dummy;
      if (addObsRatio){
        if (tgdata_withZeros->GetN()!=tg_systband->GetN()) cerr << "Number of bins for data with zeros and syst. band are not the same." << endl;
        else{
          hdummy_ratio = (TH1F*) hsum_nominal->Clone("hframe_ratio"); hdummy_ratio->Reset("ICESM"); hdummy_ratio->SetLineColor(0); hdummy_ratio->SetMarkerColor(0); hdummy_ratio->SetLineWidth(1);
          tg_systband_unit = (TGraphAsymmErrors*) tg_systband->Clone(Form("%s_unit", tg_systband->GetName()));
          tgdata_withZeros_unit = (TGraphAsymmErrors*) tgdata_withZeros->Clone(Form("%s_unit", tgdata_withZeros->GetName()));
          for (auto const& htmp:hlist_ALT) hlist_ALT_dummy.push_back((TH1F*) htmp->Clone(Form("%s_unit", htmp->GetName())));

          double ymin_ratio = 9e9;
          double ymax_ratio = -9e9;
          for (int ix=0; ix<tgdata_withZeros->GetN(); ix++){
            double& val_systband = tg_systband_unit->GetY()[ix];
            double& val_systband_errdn = tg_systband_unit->GetEYlow()[ix];
            double& val_systband_errup = tg_systband_unit->GetEYhigh()[ix];
            double& val_data = tgdata_withZeros_unit->GetY()[ix];
            double& val_data_errdn = tgdata_withZeros_unit->GetEYlow()[ix];
            double& val_data_errup = tgdata_withZeros_unit->GetEYhigh()[ix];

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
          hdummy_ratio->GetYaxis()->SetRangeUser(ymin_ratio, ymax_ratio);
          hdummy_ratio->GetXaxis()->SetLabelSize(0);
          hdummy_ratio->GetXaxis()->SetTitleSize(0);
          hdummy_ratio->GetYaxis()->SetTitleFont(42);
          hdummy_ratio->GetYaxis()->SetTitleSize(npixels_XYTitle/npixels_pad_ratio);
          hdummy_ratio->GetYaxis()->SetTitle("#times#frac{1}{SM tot.}");
          hdummy_ratio->GetYaxis()->SetTitleOffset(offset_ytitle*(npixels_pad_ratio/npixels_pad_top));
          hdummy_ratio->GetYaxis()->SetNdivisions(3, 5, 2, true);
          hdummy_ratio->GetYaxis()->SetMaxDigits(2);

          hdummy_ratio->Draw("hist");
          for (auto& hh:hlist_ALT_dummy) hh->Draw("histsame");
          tg_systband_unit->Draw("2same");
          tgdata_withZeros_unit->Draw("0psame");
        }
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

      for (auto& hh:hlist_ALT_dummy) delete hh;
      delete hdummy_ratio;
      delete tg_systband_unit;
      delete tg_systband;
      delete tgdata_withZeros_unit;
      delete tgdata_withZeros;
      delete tgdata;
      for (auto it=procdist.begin(); it!=procdist.end(); it++){
        cout << "\t- Deleting histogram " << it->second.at(idim)->GetName() << endl;
        delete it->second.at(idim);
      }
      delete systband.first;
      delete systband.second;
    }

  }
}
