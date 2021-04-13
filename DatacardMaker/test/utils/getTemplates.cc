#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <unistd.h>
#include "TIterator.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include "RooCmdArg.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooRealIntegral.h"
#include "RooProdPdf.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"

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
  unordered_map<string, RooAbsPdf*> pdf_shape;

  RooAbsReal* norm;
  double rate;
  TString name;
  vector<TH2F*> templates2D;
  vector<TH3F*> templates3D;
  unordered_map<string, pair<string, string>> systematics;

  process_spec() : pdf(0), norm(0), rate(0){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_, vector<TH2F*> templates_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()), templates2D(templates_){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_, vector<TH3F*> templates_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()), templates3D(templates_){}
  process_spec(const process_spec& other, TString newname="") : pdf(other.pdf), pdf_shape(other.pdf_shape), norm(other.norm), rate(other.rate), name(newname=="" ? other.name: newname), templates2D(other.templates2D), templates3D(other.templates3D), systematics(other.systematics){}

  void setSystematic(string systname, string systtype, string systline){ systematics[systname] = pair<string, string>(systtype, systline); }
  void setShapePdf(string systname, RooAbsPdf* pdf_){ if (pdf_) pdf_shape[systname] = pdf_; else cout << name << "_" << systname << " pdf does not exist!" << endl; }
  void writeTemplates(TFile* fin){ for (auto& tpl:templates2D) fin->WriteTObject(tpl); for (auto& tpl:templates3D) fin->WriteTObject(tpl); }
  void deleteTemplates(){ for (auto& tpl:templates2D) delete tpl; for (auto& tpl:templates3D) delete tpl; }

  template<typename T> void assignTemplates(std::vector<T*> templates);

};
template<> void process_spec::assignTemplates<TH2F>(std::vector<TH2F*> templates){ templates2D=templates; }
template<> void process_spec::assignTemplates<TH3F>(std::vector<TH3F*> templates){ templates3D=templates; }


void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);
template<typename TH_t> void extractTemplates(process_spec& proc, RooDataSet* data, string shapename);
template void extractTemplates<TH2F>(process_spec& proc, RooDataSet* data, string shapename);
template void extractTemplates<TH3F>(process_spec& proc, RooDataSet* data, string shapename);

double ArcCot(double arg){
  double res = -TMath::ATan(arg);
  if (arg>=0.) res += TMath::Pi()/2.;
  else res -= TMath::Pi()/2.;
  return res;
}
float getBreitWignerIntegral(
  float mass, float width,
  float massmin, float massmax
){
  double MG = mass*width;
  double msq = pow(mass, 2);
  double msqmax = pow(massmax, 2);
  double msqmin = pow(massmin, 2);
  double int_high = -ArcCot(MG/(msq-msqmax)) / MG;
  double int_low = -ArcCot(MG/(msq-msqmin)) / MG;
  double res = int_high - int_low;
  return res;
}
float getBreitWignerIntegralRatio(
  float mass_old, float width_old,
  float mass_new, float width_new,
  float massmin, float massmax
){
  return getBreitWignerIntegral(mass_new, width_new, massmin, massmax)/getBreitWignerIntegral(mass_old, width_old, massmin, massmax);
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
TH1F* getHistogramSlice(TH2F const* histo, unsigned char XDirection, int iy, int jy, TString newname="");
TH1F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, int iy, int jy, int iz, int jz, TString newname=""); // "y" and "z" are cylical, so if Xdirection==1 (Y), "y"=Z and "z"=X
TH2F* getHistogramSlice(TH3F const* histo, unsigned char XDirection, unsigned char YDirection, int iz, int jz, TString newname="");
template <typename T> void wipeOverUnderFlows(T* hwipe, bool rescale=false);
template<> void wipeOverUnderFlows<TH1F>(TH1F* hwipe, bool rescale);
template<> void wipeOverUnderFlows<TH2F>(TH2F* hwipe, bool rescale);
template<> void wipeOverUnderFlows<TH3F>(TH3F* hwipe, bool rescale);
template <typename T> void conditionalizeHistogram(T* histo, unsigned int iaxis, std::vector<std::pair<T*, float>> const* conditionalsReference=nullptr, bool useWidth=true, bool useEffErr=false);
template<> void conditionalizeHistogram<TH2F>(TH2F* histo, unsigned int iaxis, std::vector<std::pair<TH2F*, float>> const* conditionalsReference, bool useWidth, bool useEffErr);
template<> void conditionalizeHistogram<TH3F>(TH3F* histo, unsigned int iaxis, std::vector<std::pair<TH3F*, float>> const* conditionalsReference, bool useWidth, bool useEffErr);

template<typename T> bool checkVarNanInf(T const& val);
template<typename T> bool checkVarNanInf(T const& val){
  return !(std::isnan(val) || std::isinf(val));
}
template<typename T> bool checkNanInf(std::vector<T> const& vars){
  for (T const& v:vars){ if (!checkVarNanInf<T>(v)) return false; }
  return true;
}

template <typename T> void multiplyHistograms(T const* h1, T const* h2, T*& hAssign, bool useEffErr);
template <typename T> void multiplyHistograms(T const* h1, TH1F const* h2, unsigned int matchDimension, T*& hAssign, bool useEffErr);
template<> void multiplyHistograms<TH1F>(TH1F const* h1, TH1F const* h2, TH1F*& hAssign, bool useEffErr);
template<> void multiplyHistograms<TH2F>(TH2F const* h1, TH2F const* h2, TH2F*& hAssign, bool useEffErr);
template<> void multiplyHistograms<TH3F>(TH3F const* h1, TH3F const* h2, TH3F*& hAssign, bool useEffErr);
template<> void multiplyHistograms<TH2F>(TH2F const* h1, TH1F const* h2, unsigned int matchDimension, TH2F*& hAssign, bool useEffErr);
template<> void multiplyHistograms<TH3F>(TH3F const* h1, TH1F const* h2, unsigned int matchDimension, TH3F*& hAssign, bool useEffErr);

template<typename T> void rescaleOffshellTemplates(std::vector<T*>& tpls, TString const& strSqrts, TString const& strPeriod){
  if (tpls.empty()) return;
  vector<float> binLowEdgeList;
  for (int ix=1; ix<=tpls.back()->GetNbinsX()+1; ix++) binLowEdgeList.push_back(tpls.back()->GetXaxis()->GetBinLowEdge(ix));
  if (strSqrts=="7TeV" || strSqrts=="8TeV"){
    TH1F scalingHist("tmp_scale", "", binLowEdgeList.size()-1, binLowEdgeList.data());
    TH1F scalingHist_sqrt("tmp_scale_sqrt", "", binLowEdgeList.size()-1, binLowEdgeList.data());
    for (unsigned int ix=0; ix<binLowEdgeList.size()-1; ix++){
      float const& xlow = binLowEdgeList.at(ix);
      float const& xhigh = binLowEdgeList.at(ix+1);
      float const bwscale = getBreitWignerIntegralRatio(125.6, 0.00415, 125, 0.00407, xlow, xhigh);
      //cout << "BW scale factor for bin [" << xlow << ", " << xhigh << "]: " << bwscale << endl;
      scalingHist.SetBinContent(ix+1, bwscale);
      scalingHist_sqrt.SetBinContent(ix+1, sqrt(bwscale));
    }
    const float fLQscale = pow(125.6/125., 2);
    cout << "rescaleOffshellTemplates: WARNING! Scaling for fLQ contributions. Single power scale is " << fLQscale << endl;
    for (auto& tpl:tpls){
      TString tplname = tpl->GetName();
      if (tplname.Contains("Int")) multiplyHistograms(tpl, &scalingHist_sqrt, 0, tpl, false);
      else if (tplname.Contains("Sig")) multiplyHistograms(tpl, &scalingHist, 0, tpl, false);

      // Scale fLQ contributions
      if (tplname.Contains("ai1_1")) tpl->Scale(fLQscale);
      else if (tplname.Contains("ai1_2")) tpl->Scale(pow(fLQscale, 2));
      else if (tplname.Contains("ai1_3")) tpl->Scale(pow(fLQscale, 3));
      else if (tplname.Contains("ai1_4")) tpl->Scale(pow(fLQscale, 4));
    }
  }
}

void getSqrtsPeriod(TString const& cinput, TString& strSqrtsPeriod, TString& strSqrts, TString& strPeriod){
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  TString strCWD(cwd);

  TString cinput_test = strCWD + '/' + cinput;
  if (cinput_test.Contains("7TeV")){
    strSqrts = "7TeV";
    strPeriod = "2011";
  }
  else if (cinput_test.Contains("8TeV")){
    strSqrts = "8TeV";
    strPeriod = "2012";
  }
  else if (cinput_test.Contains("13TeV_2015") || cinput_test.Contains("13_15TeV")){
    strSqrts = "13TeV";
    strPeriod = "2015";
  }
  else if (cinput_test.Contains("13TeV_2016") || cinput_test.Contains("13_16TeV")){
    strSqrts = "13TeV";
    strPeriod = "2016";
  }
  else if (cinput_test.Contains("13TeV_2017")){
    strSqrts = "13TeV";
    strPeriod = "2017";
  }
  else if (cinput_test.Contains("13TeV_2018")){
    strSqrts = "13TeV";
    strPeriod = "2018";
  }
  else{
    cerr << "Need a valid strSqrtsPeriod!" << endl;
    cerr << "\t- Input test strings: " << cinput_test << endl;
    assert(0);
  }
  strSqrtsPeriod = strSqrts + '_' + strPeriod;
}
TString getSystRename(TString const& systname, TString const& systLine, TString const& strSqrts, TString const& strPeriod, TString const& strCategory, TString const& strChannel){
  TString res=systname;
  if (res.Contains("lumi")) res = "lumiUnc";
  else if (res == "pdf_qq") res = "pdf_qqbar";
  else if (res == "QCDscale_ggVV_bonly") res = "kbkg_gg";
  else if (res == "CMS_btag_comb") res = Form("CMS_btag_comb_%s_%s", strSqrts.Data(), strPeriod.Data());
  else if (res == "CMS_eff_e") res = Form("CMS_eff_stat_e_%s_%s", strSqrts.Data(), strPeriod.Data());
  else if (res == "CMS_eff_mu" || res == "CMS_eff_m") res = Form("CMS_eff_stat_mu_%s_%s", strSqrts.Data(), strPeriod.Data());
  else if (res == "CMS_zz4mu_zjets") res = "CMS_hzz4l_zz4mu_zjets";
  else if (res == "CMS_zz4e_zjets") res = "CMS_hzz4l_zz4e_zjets";
  else if (res == "CMS_zz2e2mu_zjets") res = "CMS_hzz4l_zz2e2mu_zjets";
  else if (res == "Res4mu") res = "CMS_hzz4l_zz4mu_res";
  else if (res == "Res4e") res = "CMS_hzz4l_zz4e_res";
  else if (res == "Res2e2mu") res = "CMS_hzz4l_zz2e2mu_res";
  else if (res == "CMS_zz4l_smd_zjets_bkg_4mu") res = "CMS_hzz4l_zz4mu_shape_zjets";
  else if (res == "CMS_zz4l_smd_zjets_bkg_4e") res = "CMS_hzz4l_zz4e_shape_zjets";
  else if (res == "CMS_zz4l_smd_zjets_bkg_2e2mu") res = "CMS_hzz4l_zz2e2mu_shape_zjets";
  else if (res == "EWKcorr_VV") res = "EWcorr_VV";
  else if (res == "CMS_zz4l_ZXshape_syst"){
    if (strChannel=="4mu") res = "CMS_hzz4l_zz4mu_shape_zjets";
    else if (strChannel=="4e") res = "CMS_hzz4l_zz4e_shape_zjets";
    else if (strChannel=="2e2mu") res = "CMS_hzz4l_zz2e2mu_shape_zjets";
  }
  return res;
}
float getProcessRescale(TString const& procname, TString const& strSqrts, TString const& strPeriod, TString const& strChannel){
  float res=1;
  if (strSqrts=="13TeV" && strPeriod=="2015"){
    float ggH[2] ={ 1.218e-2, 1.234e-2 };
    float ttH[2] ={ 3.933e-4, 3.902e-4 };
    float bbH[2] ={ 0, 1.347e-4 };
    float VBFH[2] ={ 1.044e-3, 1.035e-3 };
    float ZH[2] ={ 6.677e-4, 5.643e-4 };
    float WH[2] ={ 3.788E-4, 3.800e-4 };
    if (procname.Contains("ggH")) res = (ggH[1]+ttH[1]+bbH[1])/(ggH[0]+ttH[0]+bbH[0]);
    else if (procname.Contains("gg")) res = ggH[1]/ggH[0];
    else if (procname.Contains("VBF")) res = VBFH[1]/VBFH[0];
    else if (procname.Contains("ZH")) res = ZH[1]/ZH[0];
    else if (procname.Contains("WH")) res = WH[1]/WH[0];
    else if (procname.Contains("ttH")) res = ttH[1]/ttH[0];
    else if (procname.Contains("VVH")) res = (VBFH[1]+ZH[1]+WH[1])/(VBFH[0]+ZH[0]+WH[0]);
  }
  else if (strSqrts=="8TeV" && strPeriod=="2012"){
    if (strChannel=="2e2mu"){
      if (procname=="ggH" || procname=="ttH") res = 1.063000995;
      else if (procname=="qqH" || procname=="ZH" || procname=="WH") res = 0.999383935;
      else if (procname=="bkg2d_zjets") res = 1.037345976;
    }
    else if (strChannel=="4e"){
      if (procname=="ggH" || procname=="ttH") res = 1.056010454;
      else if (procname=="qqH" || procname=="ZH" || procname=="WH") res = 1.042473623;
      else if (procname=="bkg2d_zjets") res = 1.841115359;
    }
    else if (strChannel=="4mu"){
      if (procname=="ggH" || procname=="ttH") res = 1.047930429;
      else if (procname=="qqH" || procname=="ZH" || procname=="WH") res = 1.049536928;
    }
  }
  else if (strSqrts=="7TeV" && strPeriod=="2011"){
    if (strChannel=="2e2mu"){
      if (procname=="ggH" || procname=="ttH") res = 1.061918564;
      else if (procname=="qqH" || procname=="ZH" || procname=="WH") res = 1.00309823;
      else if (procname=="bkg2d_ggzz") res = 1.024011343;
      else if (procname=="bkg2d_zjets") res = 1.037375111;
    }
    else if (strChannel=="4e"){
      if (procname=="ggH" || procname=="ttH") res = 1.050730425;
      else if (procname=="qqH" || procname=="ZH" || procname=="WH") res = 1.071638182;
      else if (procname=="bkg2d_ggzz") res = 1.021341508;
      else if (procname=="bkg2d_zjets") res = 1.840982222;
    }
    else if (strChannel=="4mu"){
      if (procname=="ggH" || procname=="ttH") res = 1.044852093;
      else if (procname=="qqH" || procname=="ZH" || procname=="WH") res = 1.071797032;
      else if (procname=="bkg2d_qqzz") res = 0.999944344;
      else if (procname=="bkg2d_ggzz") res = 1.026622352;
    }
  }
  return res;
}

bool checkProcessIsBkg(TString procname){
  procname.ToLower();
  return (procname.Contains("bkg") || procname.Contains("qqzz") || procname.Contains("zjets"));
}

RooAbsPdf* searchMassPdf(RooDataSet* data, RooAbsPdf* rawpdf){
  RooProdPdf* prodpdf = dynamic_cast<RooProdPdf*>(rawpdf);
  if (!prodpdf) return nullptr;
  RooArgList pdfList = prodpdf->pdfList();
  TIterator* pdfIter = pdfList.createIterator();
  RooAbsPdf* res=nullptr;
  RooAbsPdf* pdf;
  while ((pdf = (RooAbsPdf*) pdfIter->Next())){
    cout << "searchMassPdf: Checking if pdf " << pdf->GetName() << " might be the mass pdf..." << endl;
    RooArgSet* depList = pdf->getDependents(data);
    RooProdPdf* subprodpdf = dynamic_cast<RooProdPdf*>(pdf);
    if (subprodpdf) res=searchMassPdf(data, subprodpdf);
    if (!res && depList->getSize()==1){
      bool namecheck=false;
      TIterator* varIter = depList->createIterator();
      RooAbsReal* var;
      while ((var = (RooAbsReal*) varIter->Next())){
        TString varname_lower=var->GetName(); varname_lower.ToLower();
        namecheck |= varname_lower.Contains("mass");
      }
      delete varIter;
      if (namecheck){
        cout << "searchMassPdf: " << pdf->GetName() << " is the mass pdf!" << endl;
        res=pdf;
        break;
      }
    }
  }
  delete pdfIter;
  return res;
}

void renameDataObservables(RooWorkspace* ws, RooDataSet* data){
  RooArgSet* args = (RooArgSet*) data->get();
  if (args->getSize()==0){ cerr << "renameDataObservables: Number of observables is 0!" << endl; exit(1); }
  TIterator* coefIter = args->createIterator();
  RooAbsArg* coef;
  unsigned int ikd=0;
  while ((coef = (RooAbsArg*) coefIter->Next())){
    RooAbsReal* rvar = dynamic_cast<RooAbsReal*>(coef);
    TString rname = rvar->GetName();
    TString newname=rname;
    if (rname.Contains("mass") || rname.Contains("Mass")) newname="mass";
    else{
      newname=Form("KD%i", ikd+1);
      ikd++;
    }
    cout << "renameDataObservables: Renaming " << rname << " to " << newname << endl;
    RooRealVar* wvar = ws->var(rname);
    rvar->SetName(newname);
    rvar->SetTitle(newname);
    wvar->SetName(newname);
    wvar->SetTitle(newname);
  }
  delete coefIter;
}
void getDataTree(TString cinput, TString coutput_main){
  TString strSqrtsPeriod, strSqrts, strPeriod;
  getSqrtsPeriod(cinput, strSqrtsPeriod, strSqrts, strPeriod);

  string strinput = cinput.Data();
  vector<string> splitinput;
  splitOptionRecursive(strinput, splitinput, '/');
  strinput = splitinput.at(splitinput.size()-1); splitinput.clear();
  splitOptionRecursive(strinput, splitinput, '/');

  const unsigned int nchans=3;
  string channames[nchans]={ "4mu", "4e", "2e2mu" };
  string channame;
  for (unsigned int ic=0; ic<nchans; ic++){
    if (strinput.find(channames[ic])!=string::npos) channame=channames[ic];
  }

  const unsigned int ncats=4;
  string catnames[ncats]={ "Inclusive", "JJVBFTagged", "HadVHTagged", "Untagged" };
  string catnamesALT[ncats]={ "Inclusive", "VBFtagged", "VHHadrtagged", "Untagged" };
  string catname = catnames[0];
  {
    int icat=-1;
    for (unsigned int ic=0; ic<ncats; ic++){
      if (strinput.find(catnames[ic])!=string::npos){ catname=catnames[ic]; icat=ic; }
    }
    if (icat<0){
      for (unsigned int ic=0; ic<ncats; ic++){
        if (strinput.find(catnamesALT[ic])!=string::npos){ catname=catnames[ic]; icat=ic; }
      }
    }
  }

  TFile* finput = TFile::Open(cinput+".input.root", "read");
  RooWorkspace* ws = (RooWorkspace*)finput->Get("w");
  RooDataSet* data = (RooDataSet*)ws->data("data_obs");
  int nevents=data->sumEntries();
  int nvars=((const RooArgSet*)data->get())->getSize();
  float* KD = new float[nvars];
  float mass;

  TString coutput_data = "Decompilation/Data/" + coutput_main;
  gSystem->Exec("mkdir -p " + coutput_data);
  TString coutput_root;
  TFile* foutput;
  coutput_root = Form("%s/hzz%s_%s_%s.root", coutput_data.Data(), channame.c_str(), catname.c_str(), strSqrtsPeriod.Data());
  foutput = TFile::Open(coutput_root, "recreate");

  TTree* t = new TTree("data_obs", "");
  //for (int iv=0; iv<nvars; iv++) t->Branch(Form("KD%i", iv+1), KD+iv);

  unsigned int nKDs=0;
  for (int ev=0; ev<nevents; ev++){
    RooArgSet* args = (RooArgSet*)data->get(ev);
    TIterator* coefIter = args->createIterator();
    RooAbsArg* coef;
    unsigned int ik=0;
    unsigned int ikd=0;
    while ((coef = (RooAbsArg*)coefIter->Next())){
      RooAbsReal* rvar = dynamic_cast<RooAbsReal*>(coef);
      TString rname = rvar->GetName();
      if (rname.Contains("mass") || rname.Contains("Mass")){
        mass=rvar->getVal();
        if (t->GetBranchStatus("mass")==0) t->Branch("mass", &mass);
      }
      else{
        KD[ik]=rvar->getVal();
        TString KDname = Form("KD%i", ikd+1);
        if (t->GetBranchStatus(KDname)==0) t->Branch(KDname, KD+ik);
        ikd++;
      }
      ik++;
    }
    delete coefIter;
    t->Fill();
  }

  foutput->WriteTObject(t);
  delete t;
  delete[] KD;
  foutput->Close();
  finput->Close();
}
void getTemplates(TString cinput, TString coutput_main, double lumiScale=1, bool copy_ggH_to_VVH=false, bool rescale_xsec=false, bool rescaleOffshellComponents=false, bool hasExtMassShapes=false){
  TString strSqrtsPeriod, strSqrts, strPeriod;
  getSqrtsPeriod(cinput, strSqrtsPeriod, strSqrts, strPeriod);

  getDataTree(cinput, coutput_main);

  string strinput = cinput.Data();
  vector<string> splitinput;
  splitOptionRecursive(strinput, splitinput, '/');
  strinput = splitinput.back(); splitinput.clear();
  splitOptionRecursive(strinput, splitinput, '/');

  TString coutput_templates = "Decompilation/Templates/" + coutput_main;
  TString coutput_inputs = "Decompilation/Inputs/" + coutput_main;
  TString coutput_extshapes = "Decompilation/ExternalShapes/" + coutput_main;

  gSystem->Exec("mkdir -p " + coutput_templates);
  gSystem->Exec("mkdir -p " + coutput_inputs);
  if (hasExtMassShapes) gSystem->Exec("mkdir -p " + coutput_extshapes);

  const unsigned int nchans=3;
  string channames[nchans]={ "4mu", "4e", "2e2mu" };
  string channame;
  for (unsigned int ic=0; ic<nchans; ic++){
    if (strinput.find(channames[ic])!=string::npos) channame=channames[ic];
  }

  const unsigned int ncats=4;
  string catnames[ncats]={ "Inclusive", "JJVBFTagged", "HadVHTagged", "Untagged" };
  string catnamesALT[ncats]={ "Inclusive", "VBFtagged", "VHHadrtagged", "Untagged" };
  string catname = catnames[0];
  {
    int icat=-1;
    for (unsigned int ic=0; ic<ncats; ic++){
      if (strinput.find(catnames[ic])!=string::npos){ catname=catnames[ic]; icat=ic; }
    }
    if (icat<0){
      for (unsigned int ic=0; ic<ncats; ic++){
        if (strinput.find(catnamesALT[ic])!=string::npos){ catname=catnames[ic]; icat=ic; }
      }
    }
  }
  TString strChannel = channame.c_str();
  TString strCategory = catname.c_str();

  //TString coutput = splitinput.at(splitinput.size()-1);
  //ofstream tout((coutput+".input.txt").Data());
  TString coutput_txt = Form("%s/inputs_%s_%s%s", coutput_inputs.Data(), channame.c_str(), catname.c_str(), ".txt");
  ofstream tout(coutput_txt.Data());

  TFile* finput = TFile::Open(cinput+".input.root", "read");
  RooWorkspace* ws = (RooWorkspace*)finput->Get("w");
  TString varsToCheck[6]={
    "R", "RF", "RF_"+strSqrts, "RV", "RV_"+strSqrts, "R_"+strSqrts
  };
  for (unsigned int v=0; v<6; v++){
    if (ws->var(varsToCheck[v])) ((RooRealVar*)ws->var(varsToCheck[v]))->setVal(1);
  }
  if (lumiScale<0.){
    RooAbsReal* lumivar = (RooAbsReal*) ws->factory(Form("LUMI_%s", strSqrtsPeriod.Data())); // Use factory to get RooConstVars.
    if (lumivar) lumiScale = lumivar->getVal();
    else{
      cout << "No lumi variable is found. Setting lumi scale to 1." << endl;
      lumiScale = 1;
    }
  }
  RooRealVar* MH = ws->var("MH");
  if (MH){ // Set MH to 125
    MH->setRange(0, 13000);
    MH->setVal(125);
  }
  ws->Print("v");

  ifstream tin((cinput+".txt").Data());
  char line[512];

  RooDataSet* data = (RooDataSet*) ws->data("data_obs");

  // Get process names
  while (string(line).find("process")==string::npos) tin.getline(line, 512);
  char* chars_array = strtok(line, " ");
  chars_array = strtok(NULL, " ");
  vector<TString> procname;
  while (chars_array){
    procname.push_back(TString(chars_array));
    chars_array = strtok(NULL, " ");
  }

  // Get process rates
  while (string(line).find("rate")==string::npos || string(line).find('#')==0) tin.getline(line, 512);
  chars_array = strtok(line, " ");
  chars_array = strtok(NULL, " ");
  vector<double> procrate;
  while (chars_array){
    procrate.push_back(atof(chars_array));
    chars_array = strtok(NULL, " ");
  }

  // Get process pdfs
  unordered_map<string, process_spec> procSpecs;
  for (unsigned int ip=0; ip<procname.size(); ip++){
    cout << procname.at(ip) << ": " << procrate.at(ip) << endl;
    RooAbsPdf* pdf = ws->pdf(procname.at(ip));
    RooAbsReal* norm = nullptr; norm = (RooAbsReal*)ws->factory(procname.at(ip)+"_norm");
    if (norm) cout << "Warning: " << procname.at(ip) << "_norm is not found." << endl;
    if (!pdf) cerr << procname.at(ip) << " pdf could not be found." << endl;
    else procSpecs[procname.at(ip).Data()]=process_spec(pdf, norm, procrate.at(ip));
  }

  // Get systematics
  unordered_map<string, string> tplSyst;
  unordered_map<string, string> logSyst;
  unordered_map<string, string> paramSyst;
  while (!tin.eof()){
    tin.getline(line, 512);

    string strline = line;
    std::replace(strline.begin(), strline.end(), ',', ' ');
    strline.erase(std::remove(strline.begin(), strline.end(), '['), strline.end());
    strline.erase(std::remove(strline.begin(), strline.end(), ']'), strline.end());
    bool isShape = strline.find("shape1")!=string::npos;
    bool isParam = strline.find("param")!=string::npos;
    bool isLog = strline.find("lnN")!=string::npos;
    if (isShape || isLog || isParam){
      vector<string> systdist;
      splitOptionRecursive(strline, systdist, ' ');
      string systname = systdist.at(0);
      string systtype = systdist.at(1);
      string accumulate="";
      cout << "Processing systematic " << systname << endl;
      if (systname == "kbkg_gg"){
        // Directly add since this is a special systematic
        for (unsigned int ip=2; ip<systdist.size(); ip++){
          if (ip>2) accumulate += ":";
          accumulate += systdist.at(ip);
        }
        paramSyst[systname] = accumulate;
      }
      else if (isShape || isLog){
        for (unsigned int ip=0; ip<procname.size(); ip++){
          string systline = systdist.at(ip+2);
          if (systline.find("-")==string::npos && systline!=""){
            std::replace(systline.begin(), systline.end(), '/', ':');
            procSpecs[procname.at(ip).Data()].setSystematic(systname, systtype, systline);
            if (isShape){
              RooAbsPdf* up_ = ws->pdf(procname.at(ip)+"_"+systname.c_str()+"Up");
              RooAbsPdf* dn_ = ws->pdf(procname.at(ip)+"_"+systname.c_str()+"Down");
              procSpecs[procname.at(ip).Data()].setShapePdf(systname+"Up", up_);
              procSpecs[procname.at(ip).Data()].setShapePdf(systname+"Down", dn_);
            }
            if (isShape) accumulate += string(procname.at(ip).Data()) + ":0:1 ";
            else accumulate += string(procname.at(ip).Data()) + ":" + systline + " ";
            if (copy_ggH_to_VVH && procname.at(ip)=="ggH"){
              if (isShape) accumulate += string("VVH") + ":0:1 ";
              else accumulate += string("VVH") + ":" + systline + " ";
            }
          }
        }
        if (isLog) logSyst[systname] = accumulate;
        else if (isShape) tplSyst[systname] = accumulate;
      }
      else{
        RooAbsReal* systvar = (RooAbsReal*) ws->var(systname.c_str());
        double defaultVal = systvar->getVal();
        if (systvar->hasClients()){
          for (unsigned int ip=0; ip<procname.size(); ip++){
            RooAbsPdf* pdf = procSpecs[procname.at(ip).Data()].pdf;
            RooAbsReal* norm = procSpecs[procname.at(ip).Data()].norm;
            if (
              (pdf && pdf->dependsOn(*systvar))
              ||
              (norm && norm->dependsOn(*systvar))
              ){
              cout << "Pdf or norm for " << procname.at(ip) << " depends on systematic " << systname << endl;
              string systline="";
              for (unsigned int ip=2; ip<systdist.size(); ip++){
                if (ip>2) systline += ":";
                systline += systdist.at(ip);
              }
              procSpecs[procname.at(ip).Data()].setSystematic(systname, systtype, systline);
              accumulate += string(procname.at(ip).Data()) + ":" + systline + " ";
              if (copy_ggH_to_VVH && procname.at(ip)=="ggH") accumulate += string("VVH") + ":" + systline + " ";
            }
          }
        }
        paramSyst[systname] = accumulate;
      }

    }
  }
  cout << "Processes: ";
  for (auto const& pname:procname) cout << pname << " ";
  cout << endl;
  cout << "Process specification keys: ";
  for (auto it=procSpecs.begin(); it!=procSpecs.end(); it++) cout << it->first << " ";
  cout << endl;
  if (copy_ggH_to_VVH){
    for (auto it=procSpecs.begin(); it!=procSpecs.end(); it++){
      if (TString(it->first)=="ggH"){
        cout << "Copying ggH process into VVH" << endl;
        procSpecs["VVH"]=process_spec(it->second, "VVH");
        procname.emplace_back("VVH");
        cout << "Copy successful" << endl;
      }
    }
  }

  // Write input file
  if (strSqrts!=""){
    TString strSqrtsBare=strSqrts;
    int ipos = strSqrtsBare.Index("TeV");
    if(ipos>=0) strSqrtsBare.Resize(ipos);
    tout << "sqrts " << strSqrtsBare << endl;
  }
  if (strPeriod!="") tout << "period " << strPeriod << endl;
  tout << "decay " << channame << endl;
  tout << "lumi " << lumiScale << endl;
  tout << "category " << catname << endl;

  // Write channels
  for (unsigned int ip=0; ip<procname.size(); ip++){
    unsigned int proccode = (checkProcessIsBkg(procSpecs[procname.at(ip).Data()].name) ? 0 : (procSpecs[procname.at(ip).Data()].name.Contains("offshell") ? 2 : 1));
    tout << "channel " << procname.at(ip) << " 1 -1 " << proccode << endl;
  }

  // Write systematics
  /*
  for (unsigned int ip=0; ip<procname.size(); ip++){
    for (auto syst = procSpecs[procname.at(ip).Data()].systematics.begin(); syst != procSpecs[procname.at(ip).Data()].systematics.end(); ++syst){
      if (syst->second.first=="shape1"){
        tout << "systematic " << syst->first << " template ";
        tout << procname.at(ip) << ":0:1" << endl;
      }
    }
  }
  */
  for (auto syst = tplSyst.begin(); syst != tplSyst.end(); ++syst){
    TString systName = syst->first.c_str();
    TString systLine = syst->second.c_str();
    systName = getSystRename(systName, systLine, strSqrts, strPeriod, strCategory, strChannel);
    if (syst->second!="") tout << "systematic " << systName << " template " << syst->second << endl;
  }
  for (auto syst = logSyst.begin(); syst != logSyst.end(); ++syst){
    TString systName = syst->first.c_str();
    TString systLine = syst->second.c_str();
    systName = getSystRename(systName, systLine, strSqrts, strPeriod, strCategory, strChannel);
    if (syst->second!="") tout << "systematic " << systName << " lnN " << syst->second << endl;
  }
  for (auto syst = paramSyst.begin(); syst != paramSyst.end(); ++syst){
    TString systName = syst->first.c_str();
    TString systLine = syst->second.c_str();
    systName = getSystRename(systName, systLine, strSqrts, strPeriod, strCategory, strChannel);
    if (syst->second!=""){
      if (systName == "kbkg_gg") tout << "systematic " << systName << " param " << syst->second << endl;
      else tout << "systematic " << systName << " template " << syst->second << endl;
    }
  }

  // Search for external shapes but do not record them
  unordered_map<TString,RooAbsPdf*> extShapeProcPdfs;
  if (hasExtMassShapes){
    for (unsigned int ip=0; ip<procname.size(); ip++){
      cout << "Attempting to search for a mass pdf for process " << procname.at(ip) << endl;
      RooAbsPdf* mass_pdf = searchMassPdf(data, procSpecs[procname.at(ip).Data()].pdf);
      if (mass_pdf) extShapeProcPdfs[procname.at(ip)] = mass_pdf;
    }
  }


  for (unsigned int ip=0; ip<procname.size(); ip++){
    cout << "Attempting to extract tpls for process " << procname.at(ip) << endl;

    const unsigned int ndims = procSpecs[procname.at(ip).Data()].pdf->getDependents(data)->getSize();
    cout << "Number of dimensions: " << ndims << endl;

    TString const& theProcName = procname.at(ip);
    bool isBkg = checkProcessIsBkg(theProcName);
    TString theProcNameLower = theProcName; theProcNameLower.ToLower();
    float tplscale=1;
    if (!theProcNameLower.Contains("zjets")) tplscale=1./lumiScale;
    if (rescale_xsec){
      float extrascale = getProcessRescale(theProcName, strSqrts, strPeriod, channame.c_str());
      cout << "Extra scale factor: " << extrascale << endl;
      tplscale *= extrascale;
    }

    TString coutput_root = Form("%s/Hto%s_%s_FinalTemplates_%s_%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procname.at(ip).Data(), "Nominal", ".root");
    TFile* foutput = TFile::Open(coutput_root, "recreate");
    switch (ndims){
    case 2:
      extractTemplates<TH2F>(procSpecs[procname.at(ip).Data()], data, "");
      for (auto& tpl:procSpecs[procname.at(ip).Data()].templates2D) tpl->Scale(tplscale);
      break;
    case 3:
      extractTemplates<TH3F>(procSpecs[procname.at(ip).Data()], data, "");
      for (auto& tpl:procSpecs[procname.at(ip).Data()].templates3D) tpl->Scale(tplscale);
      break;
    }
    if (rescaleOffshellComponents && !isBkg){
      rescaleOffshellTemplates(procSpecs[procname.at(ip).Data()].templates2D, strSqrts, strPeriod);
      rescaleOffshellTemplates(procSpecs[procname.at(ip).Data()].templates3D, strSqrts, strPeriod);
    }
    procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
    if (extShapeProcPdfs.find(procname.at(ip))!=extShapeProcPdfs.cend()){
      vector<TH2F*> extShapeTpls2D;
      vector<TH3F*> extShapeTpls3D;
      for (auto& tpl:procSpecs[procname.at(ip).Data()].templates2D){
        TString tplname = tpl->GetName();
        if (tplname==Form("T_%s", procname.at(ip).Data()) || tplname==Form("T_%s_Bkg", procname.at(ip).Data()) || tplname==Form("T_%s_Sig", procname.at(ip).Data())){
          TH2F* tplcopy=new TH2F(*tpl);
          tplcopy->SetName(Form("%s_condDim0", tplname.Data()));
          conditionalizeHistogram<TH2F>(tplcopy, 0, nullptr, true, false);
          extShapeTpls2D.push_back(tplcopy);
        }
      }
      for (auto& tpl:procSpecs[procname.at(ip).Data()].templates3D){
        TString tplname = tpl->GetName();
        if (tplname==Form("T_%s", procname.at(ip).Data()) || tplname==Form("T_%s_Bkg", procname.at(ip).Data()) || tplname==Form("T_%s_Sig", procname.at(ip).Data())){
          TH3F* tplcopy=new TH3F(*tpl);
          tplcopy->SetName(Form("%s_condDim0", tplname.Data()));
          conditionalizeHistogram<TH3F>(tplcopy, 0, nullptr, true, false);
          extShapeTpls3D.push_back(tplcopy);
        }
      }
      for (auto& tpl:extShapeTpls2D){ foutput->WriteTObject(tpl); delete tpl; }
      for (auto& tpl:extShapeTpls3D){ foutput->WriteTObject(tpl); delete tpl; }
    }
    procSpecs[procname.at(ip).Data()].deleteTemplates();
    foutput->Close();

    for (auto syst = procSpecs[procname.at(ip).Data()].systematics.begin(); syst != procSpecs[procname.at(ip).Data()].systematics.end(); ++syst){
      if (syst->first == "kbkg_gg") continue; // kbkg_gg is special, so skip it here.
      if (syst->second.first=="param"){
        RooRealVar* systvar = (RooRealVar*)ws->var(syst->first.c_str());
        if (!systvar){
          cout << syst->first << " could not be found." << endl;
          continue;
        }
        double systCentralVal = systvar->getVal();
        for (unsigned int is=0; is<2; is++){
          systvar->setVal(double(2*is)-1.+systCentralVal);
          cout << "Setting param systematic " << systvar->GetName() << " to " << systvar->getVal() << endl;
          TString systName = syst->first.c_str();
          TString systLine = (syst->second.first + " " + syst->second.second).c_str();
          systName = getSystRename(systName, systLine, strSqrts, strPeriod, strCategory, strChannel);
          coutput_root = Form("%s/Hto%s_%s_FinalTemplates_%s_%s%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procname.at(ip).Data(), systName.Data(), (is==0 ? "Down" : "Up"), ".root");
          foutput = TFile::Open(coutput_root, "recreate");
          switch (ndims){
          case 2:
            extractTemplates<TH2F>(procSpecs[procname.at(ip).Data()], data, "");
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates2D) tpl->Scale(tplscale);
            break;
          case 3:
            extractTemplates<TH3F>(procSpecs[procname.at(ip).Data()], data, "");
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates3D) tpl->Scale(tplscale);
            break;
          }
          if (rescaleOffshellComponents && !isBkg){
            rescaleOffshellTemplates(procSpecs[procname.at(ip).Data()].templates2D, strSqrts, strPeriod);
            rescaleOffshellTemplates(procSpecs[procname.at(ip).Data()].templates3D, strSqrts, strPeriod);
          }
          procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
          if (extShapeProcPdfs.find(procname.at(ip))!=extShapeProcPdfs.cend()){
            vector<TH2F*> extShapeTpls2D;
            vector<TH3F*> extShapeTpls3D;
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates2D){
              TString tplname = tpl->GetName();
              if (tplname==Form("T_%s", procname.at(ip).Data()) || tplname==Form("T_%s_Bkg", procname.at(ip).Data()) || tplname==Form("T_%s_Sig", procname.at(ip).Data())){
                TH2F* tplcopy=new TH2F(*tpl);
                tplcopy->SetName(Form("%s_condDim0", tplname.Data()));
                conditionalizeHistogram<TH2F>(tplcopy, 0, nullptr, true, false);
                extShapeTpls2D.push_back(tplcopy);
              }
            }
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates3D){
              TString tplname = tpl->GetName();
              if (tplname==Form("T_%s", procname.at(ip).Data()) || tplname==Form("T_%s_Bkg", procname.at(ip).Data()) || tplname==Form("T_%s_Sig", procname.at(ip).Data())){
                TH3F* tplcopy=new TH3F(*tpl);
                tplcopy->SetName(Form("%s_condDim0", tplname.Data()));
                conditionalizeHistogram<TH3F>(tplcopy, 0, nullptr, true, false);
                extShapeTpls3D.push_back(tplcopy);
              }
            }
            for (auto& tpl:extShapeTpls2D){ foutput->WriteTObject(tpl); delete tpl; }
            for (auto& tpl:extShapeTpls3D){ foutput->WriteTObject(tpl); delete tpl; }
          }
          procSpecs[procname.at(ip).Data()].deleteTemplates();
          foutput->Close();
        }
        systvar->setVal(systCentralVal);
      }
      else if (syst->second.first=="shape1"){
        for (unsigned int is=0; is<2; is++){
          string syst_du = syst->first + (is==0 ? "Down" : "Up");
          TString systName = syst->first.c_str();
          TString systLine = (syst->second.first + " " + syst->second.second).c_str();
          systName = getSystRename(systName, systLine, strSqrts, strPeriod, strCategory, strChannel);
          coutput_root = Form("%s/Hto%s_%s_FinalTemplates_%s_%s%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procname.at(ip).Data(), systName.Data(), (is==0 ? "Down" : "Up"), ".root");
          foutput = TFile::Open(coutput_root, "recreate");
          switch (ndims){
          case 2:
            extractTemplates<TH2F>(procSpecs[procname.at(ip).Data()], data, syst_du);
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates2D) tpl->Scale(tplscale);
            break;
          case 3:
            extractTemplates<TH3F>(procSpecs[procname.at(ip).Data()], data, syst_du);
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates3D) tpl->Scale(tplscale);
            break;
          }
          if (rescaleOffshellComponents && !isBkg){
            rescaleOffshellTemplates(procSpecs[procname.at(ip).Data()].templates2D, strSqrts, strPeriod);
            rescaleOffshellTemplates(procSpecs[procname.at(ip).Data()].templates3D, strSqrts, strPeriod);
          }
          procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
          if (extShapeProcPdfs.find(procname.at(ip))!=extShapeProcPdfs.cend()){
            vector<TH2F*> extShapeTpls2D;
            vector<TH3F*> extShapeTpls3D;
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates2D){
              TString tplname = tpl->GetName();
              if (tplname==Form("T_%s", procname.at(ip).Data()) || tplname==Form("T_%s_Bkg", procname.at(ip).Data()) || tplname==Form("T_%s_Sig", procname.at(ip).Data())){
                TH2F* tplcopy=new TH2F(*tpl);
                tplcopy->SetName(Form("%s_condDim0", tplname.Data()));
                conditionalizeHistogram<TH2F>(tplcopy, 0, nullptr, true, false);
                extShapeTpls2D.push_back(tplcopy);
              }
            }
            for (auto& tpl:procSpecs[procname.at(ip).Data()].templates3D){
              TString tplname = tpl->GetName();
              if (tplname==Form("T_%s", procname.at(ip).Data()) || tplname==Form("T_%s_Bkg", procname.at(ip).Data()) || tplname==Form("T_%s_Sig", procname.at(ip).Data())){
                TH3F* tplcopy=new TH3F(*tpl);
                tplcopy->SetName(Form("%s_condDim0", tplname.Data()));
                conditionalizeHistogram<TH3F>(tplcopy, 0, nullptr, true, false);
                extShapeTpls3D.push_back(tplcopy);
              }
            }
            for (auto& tpl:extShapeTpls2D){ foutput->WriteTObject(tpl); delete tpl; }
            for (auto& tpl:extShapeTpls3D){ foutput->WriteTObject(tpl); delete tpl; }
          }
          procSpecs[procname.at(ip).Data()].deleteTemplates();
          foutput->Close();
        }
      }


    }
  }

  if (hasExtMassShapes){
    renameDataObservables(ws, data);
    TString coutput_shapes = Form("%s/Hto%s_%s_FinalMassShape_%s%s", coutput_extshapes.Data(), channame.c_str(), catname.c_str(), "AllProcesses", ".root");
    TFile* foutput_extshapes = TFile::Open(coutput_shapes, "recreate");
    RooWorkspace wws("w", "");
    for (unsigned int ip=0; ip<procname.size(); ip++){
      const TString& theProcName=procname.at(ip);
      cout << "Attempting to extract mass pdf for process " << procname.at(ip) << endl;
      RooAbsPdf* mass_pdf = searchMassPdf(data, procSpecs[theProcName.Data()].pdf);
      if (mass_pdf){
        if (
          !checkProcessIsBkg(theProcName) && (
            theProcName.Contains("gg") || theProcName.Contains("tt") || theProcName.Contains("bb")
            ||
            theProcName.Contains("qqH") || theProcName.Contains("VBF") || theProcName.Contains("ZH") || theProcName.Contains("WH")
            )
          ) mass_pdf->SetName(theProcName+"_Sig_MassShape");
        else mass_pdf->SetName(theProcName+"_MassShape");
        wws.import(*mass_pdf, RecycleConflictNodes());
      }
    }
    foutput_extshapes->WriteTObject(&wws);
    foutput_extshapes->Close();
  }

  tin.close();
  finput->Close();
  tout.close();
}


template<typename TH_t> void extractTemplates_ggLike_fai1(const TString& newname, const vector<TH_t*>& intpl, vector<TH_t*>& outtpl){
  double invA[3][3]={
    { 1, 0, 0 },
  { -1, 2, -1 },
  { 0, 0, 1 }
  };

  if (intpl.size()!=3) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname = Form("%s_Sig", newname.Data());
    if (ot>0) nname += Form("_ai1_%i", ot);
    if (ot==1) nname += "_Re";
    outtpl.push_back((TH_t*) intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}
template<typename TH_t> void extractTemplates_VVLike_fai1(const TString& newname, const vector<TH_t*>& intpl, vector<TH_t*>& outtpl){
  double A[25]={
    1, 0, 0, 0, 0,
    9./16., 3.*sqrt(3.)/16., 3./16., sqrt(3.)/16., 1./16.,
    0.25, 0.25, 0.25, 0.25, 0.25,
    1./16., sqrt(3.)/16., 3./16., 3.*sqrt(3.)/16., 9./16.,
    0, 0, 0, 0, 1
  };
  TMatrixD matA(5, 5, A);
  TMatrixD invA = matA.Invert();

  if (intpl.size()!=5) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname = Form("%s_Sig", newname.Data());
    if (ot>0) nname += Form("_ai1_%i", ot);
    if (ot==1 || ot==3) nname += "_Re";
    else if (ot==2) nname += "_PosDef";
    outtpl.push_back((TH_t*) intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}
template<typename TH_t> void extractTemplates_ggORvv_GGsm(const TString& newname, const vector<TH_t*>& intpl, vector<TH_t*>& outtpl){
  double invA[3][3]={
    { 1, 0, 0 },
  { -1.5, 2, -0.5 },
  { 0.5, -1, 0.5 }
  };

  if (intpl.size()!=3) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname;
    switch (ot){
    case 0:
      nname = Form("%s_Bkg", newname.Data());
      break;
    case 1:
      nname = Form("%s_Int_Re", newname.Data());
      break;
    case 2:
      nname = Form("%s_Sig", newname.Data());
      break;
    default:
      assert(0);
    }
    outtpl.push_back((TH_t*) intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}
template<typename TH_t> void extractTemplates_ggLike_fai1_GGsm(const TString& newname, const vector<TH_t*>& intpl, vector<TH_t*>& outtpl){
  double invA[6][6]={
    { 1, 0, 0, 0, 0, 0 },
  { 0.5, -1, 0.5, 0, 0, 0 },
  { -1.5, 2, -0.5, 0, 0, 0 },
  { 0.5, 0, 0, -1, 0.5, 0 },
  { 5./4., -2, 1./6., -7./4., 1./4., 25./12. },
  { -1.5, 0, 0, 2, -0.5, 0 }
  };

  if (intpl.size()!=6) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname;
    switch (ot){
    case 0:
      nname = Form("%s_Bkg", newname.Data());
      break;
    case 1:
      nname = Form("%s_Sig", newname.Data());
      break;
    case 2:
      nname = Form("%s_Int_Re", newname.Data());
      break;
    case 3:
      nname = Form("%s_Sig_ai1_2", newname.Data());
      break;
    case 4:
      nname = Form("%s_Sig_ai1_1_Re", newname.Data());
      break;
    case 5:
      nname = Form("%s_Int_ai1_1_Re", newname.Data());
      break;
    }
    outtpl.push_back((TH_t*) intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}
template<typename TH_t> void extractTemplates_VVLike_fai1_GGsm(const TString& newname, const vector<TH_t*>& intpl, vector<TH_t*>& outtpl){
  double invA[9][9]={
    { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1./2., -1, 1./2., 0, 0, 0, 0, 0, 0 },
    { -3./2., 2, -1./2., 0, 0, 0, 0, 0, 0 },
    { 1./2., -2, 1./2., 1, 50./7., 0, -50./7., -25./14., 25./14. },
    { 0, -19./60., -25./24., -5./12., -6325./336., 28561./2640., 8075./924., 50./21., -75./56. },
    { 1, 2, 1, 269./144., 16325./448., -142805./6336., -13700./693., -25./14., 25./14. },
    { 0, 53./30., -25./24., -5./2., -8825./336., 28561./2640., 7475./462., 50./21., -75./56. },
    { 0, -25./6., 25./24., 0, 200./21., 0, -75./14., -50./21., 75./56. },
    { -3./2., 2, -1./2., 0, -50./7., 0, 50./7., 25./14., -25./14. }
  };

  if (intpl.size()!=9) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname;
    switch(ot){
    case 0:
      nname="Bkg";
      break;
    case 1:
      nname="Sig";
      break;
    case 2:
      nname="Int_Re";
      break;
    case 3:
      nname="Sig_ai1_4";
      break;
    case 4:
      nname="Sig_ai1_1_Re";
      break;
    case 5:
      nname="Sig_ai1_2_PosDef";
      break;
    case 6:
      nname="Sig_ai1_3_Re";
      break;
    case 7:
      nname="Int_ai1_1_Re";
      break;
    case 8:
      nname="Int_ai1_2_Re";
      break;
    }
    nname = Form("%s_%s", newname.Data(), nname.Data());
    outtpl.push_back((TH_t*) intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}

template void extractTemplates_ggLike_fai1(const TString& newname, const vector<TH2F*>& intpl, vector<TH2F*>& outtpl);
template void extractTemplates_VVLike_fai1(const TString& newname, const vector<TH2F*>& intpl, vector<TH2F*>& outtpl);
template void extractTemplates_ggORvv_GGsm(const TString& newname, const vector<TH2F*>& intpl, vector<TH2F*>& outtpl);
template void extractTemplates_ggLike_fai1_GGsm(const TString& newname, const vector<TH2F*>& intpl, vector<TH2F*>& outtpl);
template void extractTemplates_VVLike_fai1_GGsm(const TString& newname, const vector<TH2F*>& intpl, vector<TH2F*>& outtpl);

template void extractTemplates_ggLike_fai1(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);
template void extractTemplates_VVLike_fai1(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);
template void extractTemplates_ggORvv_GGsm(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);
template void extractTemplates_ggLike_fai1_GGsm(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);
template void extractTemplates_VVLike_fai1_GGsm(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);


template<typename TH_t> void extractTemplates(process_spec& proc, RooDataSet* data, string shapename){
  vector<TH_t*> templates;

  vector<RooRealVar*> deps;
  if (!proc.pdf){
    cerr << "extractTemplates ERROR: PDF IS NULL!" << endl;
    assert(0);
  }
  else cout << "extractTemplates(" << proc.pdf->GetName() << "," << (shapename=="" ? "\"\"" : shapename) << ")" << endl;
  RooArgSet* depList = proc.pdf->getDependents(data);
  depList->Print("v");
  TIterator* coefIter = depList->createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())) deps.push_back((RooRealVar*)coef);
  delete coefIter;
  delete depList;

  vector<RooRealVar*> pars;
  RooRealVar* RV=0;
  RooRealVar* RF=0;
  RooRealVar* fai1=0;
  RooRealVar* GGsm=0;
  RooArgSet* parList;
  if (proc.norm) parList = proc.norm->getParameters(data);
  else parList = proc.pdf->getParameters(data);
  parList->Print("v");
  coefIter = parList->createIterator();
  coef=0;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    pars.push_back((RooRealVar*)coef);
    if (TString(pars.back()->GetName()).Contains("fai1")) fai1 = pars.back();
    else if (TString(pars.back()->GetName()).Contains("GGsm")) GGsm = pars.back();
    else if (TString(pars.back()->GetName())=="RV") RV = pars.back();
    else if (TString(pars.back()->GetName())=="RF") RF = pars.back();
  }
  delete coefIter;
  delete parList;
  cout << "extractTemplates: The following parameters were found as part of pdf normalization components: ";
  if (RV) cout << "RV ";
  if (RF) cout << "RF ";
  if (fai1) cout << "fai1 ";
  if (GGsm) cout << "GGsm ";
  cout << endl;

  RooCmdArg ycmd;
  RooCmdArg zcmd;
  RooCmdArg condObs;
  if (deps.size()>1){
    ycmd = YVar(*(deps.at(1)));
    //condObs = ConditionalObservables(RooArgSet(*(deps.at(0)), *(deps.at(1))));
  }
  if (deps.size()>2){
    zcmd = ZVar(*(deps.at(2)));
    //condObs = ConditionalObservables(RooArgSet(*(deps.at(0)), *(deps.at(1)), *(deps.at(2))));
  }
  TString tplname = "T_";
  tplname += proc.name;

  if (checkProcessIsBkg(proc.name)){
    TH_t* tpl;
    switch (deps.size()){
    case 2:
      tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
      break;
    case 3:
      tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
      break;
    }
    if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
    wipeOverUnderFlows(tpl, false);
    multiplyBinWidth(tpl);
    double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
    double integral = tpl->Integral();
    double scale = normval/integral;
    cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
    if (shapename!=""){
      delete tpl;
      switch (deps.size()){
      case 2:
        tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
        break;
      case 3:
        tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
        break;
      }
    }
    tpl->SetName(tplname);
    tpl->SetTitle("");
    tpl->Scale(scale);
    templates.push_back(tpl);
  }
  else if (proc.name.Contains("gg") || proc.name.Contains("tt") || proc.name.Contains("bb")){
    float RVdefval=0;
    if (RV){
      cout << "\t- Setting RV to 0" << endl;
      RVdefval = RV->getVal();
      RV->setVal(0);
    }
    else{
      cout << "Could not find RV" << endl;
    }
    if (!fai1 && !GGsm){
      TH_t* tpl;
      switch (deps.size()){
      case 2:
        tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
        break;
      case 3:
        tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
        break;
      }
      if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
      wipeOverUnderFlows(tpl, false);
      multiplyBinWidth(tpl);
      double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
      double integral = tpl->Integral();
      double scale = normval/integral;
      cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
      if (shapename!=""){
        delete tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
      }
      tpl->SetName(tplname+"_Sig");
      tpl->SetTitle("");
      tpl->Scale(scale);
      templates.push_back(tpl);
    }
    else if (fai1 && !GGsm){
      vector<TH_t*> intpl;
      for (unsigned int ifv=0; ifv<3; ifv++){
        fai1->setVal(0.5*float(ifv));

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH_t* tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
        if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
        wipeOverUnderFlows(tpl, false);
        multiplyBinWidth(tpl);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          switch (deps.size()){
          case 2:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, condObs);
            break;
          case 3:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd, condObs);
            break;
          }
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      fai1->setVal(0);
      extractTemplates_ggLike_fai1(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
    else if (!fai1 && GGsm){
      vector<TH_t*> intpl;
      for (unsigned int ifv=0; ifv<3; ifv++){
        switch (ifv){
        case 0:
        case 1:
          GGsm->setVal(ifv);
          break;
        case 2:
          GGsm->setVal(4);
          break;
        }

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH_t* tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
        if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
        wipeOverUnderFlows(tpl, false);
        multiplyBinWidth(tpl);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          switch (deps.size()){
          case 2:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, condObs);
            break;
          case 3:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd, condObs);
            break;
          }
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      GGsm->setVal(1);
      extractTemplates_ggORvv_GGsm(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
    else{
      vector<TH_t*> intpl;
      for (unsigned int ifv=0; ifv<6; ifv++){
        switch (ifv){
        case 0:
          GGsm->setVal(0);
          break;
        case 1:
          GGsm->setVal(1);
          break;
        case 2:
          GGsm->setVal(4);
          break;
        case 3:
          GGsm->setVal(1);
          fai1->setVal(1);
          break;
        case 4:
          GGsm->setVal(4);
          fai1->setVal(1);
          break;
        case 5:
          GGsm->setVal(1);
          fai1->setVal(9./25.);
          break;
        }

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH_t* tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
        if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
        wipeOverUnderFlows(tpl, false);
        multiplyBinWidth(tpl);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          switch (deps.size()){
          case 2:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, condObs);
            break;
          case 3:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd, condObs);
            break;
          }
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      GGsm->setVal(1);
      fai1->setVal(0);
      extractTemplates_ggLike_fai1_GGsm(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
    if (RV) RV->setVal(RVdefval);
  }
  else if (proc.name.Contains("qqH") || proc.name.Contains("VBF") || proc.name.Contains("vbf") || proc.name.Contains("WH") || proc.name.Contains("ZH") || proc.name.Contains("VV")){
    float RFdefval=0;
    if (RF){
      cout << "\t- Setting RF to 0" << endl;
      RFdefval = RF->getVal();
      RF->setVal(0);
    }
    else{
      cout << "Could not find RF" << endl;
    }
    if (!fai1 && !GGsm){
      TH_t* tpl;
      switch (deps.size()){
      case 2:
        tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
        break;
      case 3:
        tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
        break;
      }
      if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
      wipeOverUnderFlows(tpl, false);
      multiplyBinWidth(tpl);
      double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
      double integral = tpl->Integral();
      double scale = normval/integral;
      cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
      if (shapename!=""){
        delete tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
      }
      tpl->SetName(tplname+"_Sig");
      tpl->SetTitle("");
      tpl->Scale(scale);
      templates.push_back(tpl);
    }
    else if (fai1 && !GGsm){
      vector<TH_t*> intpl;
      for (unsigned int ifv=0; ifv<5; ifv++){
        fai1->setVal(0.25*float(ifv));

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH_t* tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
        if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
        wipeOverUnderFlows(tpl, false);
        multiplyBinWidth(tpl);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          switch (deps.size()){
          case 2:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, condObs);
            break;
          case 3:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd, condObs);
            break;
          }
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      fai1->setVal(0);
      extractTemplates_VVLike_fai1(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
    else if (!fai1 && GGsm){
      vector<TH_t*> intpl;
      for (unsigned int ifv=0; ifv<3; ifv++){
        switch (ifv){
        case 0:
        case 1:
          GGsm->setVal(ifv);
          break;
        case 2:
          GGsm->setVal(4);
          break;
        }

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH_t* tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
        if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
        wipeOverUnderFlows(tpl, false);
        multiplyBinWidth(tpl);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          switch (deps.size()){
          case 2:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, condObs);
            break;
          case 3:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd, condObs);
            break;
          }
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      GGsm->setVal(1);
      extractTemplates_ggORvv_GGsm(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
    else{
      vector<TH_t*> intpl;
      for (unsigned int ifv=0; ifv<9; ifv++){
        switch (ifv){
        case 0:
          GGsm->setVal(0);
          break;
        case 1:
          GGsm->setVal(1);
          break;
        case 2:
          GGsm->setVal(4);
          break;
        case 3:
          GGsm->setVal(1);
          fai1->setVal(1);
          break;
        case 4:
          GGsm->setVal(1);
          fai1->setVal(9./25.);
          break;
        case 5:
          GGsm->setVal(1);
          fai1->setVal(25./169.);
          break;
        case 6:
          GGsm->setVal(1);
          fai1->setVal(16./25.);
          break;
        case 7:
          GGsm->setVal(4);
          fai1->setVal(9./25.);
          break;
        case 8:
          GGsm->setVal(4);
          fai1->setVal(16./25.);
          break;
        }

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH_t* tpl;
        switch (deps.size()){
        case 2:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, condObs);
          break;
        case 3:
          tpl=(TH_t*) proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd, condObs);
          break;
        }
        if (!tpl){ cerr << "extractTemplates: Template construction failed!" << endl; assert(0); }
        wipeOverUnderFlows(tpl, false);
        multiplyBinWidth(tpl);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          switch (deps.size()){
          case 2:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, condObs);
            break;
          case 3:
            tpl=(TH_t*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd, condObs);
            break;
          }
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      GGsm->setVal(1);
      fai1->setVal(0);
      extractTemplates_VVLike_fai1_GGsm(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
    if (RF) RF->setVal(RFdefval);
  }

  for (auto& tpl:templates) divideBinWidth(tpl);
  proc.assignTemplates(templates);
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
    if (result!=""/* && !checkListVariable(splitoptions, result)*/) splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}
Bool_t checkListVariable(const vector<string>& list, const string& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
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

      double integralinside=0, integralerrorinside=0;
      double integraloutside=0, integralerroroutside=0;
      if (
        xb[0]>=ix && xb[1]<=jx
        ){
        res=histo->IntegralAndError(xb[0], xb[1], reserror, "width");
        integralinside=histo->IntegralAndError(xb[0], xb[1], integralerrorinside, "");
      }
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

      double integralinside=0, integralerrorinside=0;
      double integraloutside=0, integralerroroutside=0;
      if (
        xb[0]>=ix && xb[1]<=jx
        &&
        yb[0]>=iy && yb[1]<=jy
        ){
        res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], reserror, "width");
        integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], integralerrorinside, "");
      }
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

      double integralinside=0, integralerrorinside=0;
      double integraloutside=0, integralerroroutside=0;
      if (
        xb[0]>=ix && xb[1]<=jx
        &&
        yb[0]>=iy && yb[1]<=jy
        &&
        zb[0]>=iz && zb[1]<=jz
        ){
        res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], reserror, "width");
        integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], integralerrorinside, "");
      }
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralerrorinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template double getHistogramIntegralAndError<TH1F>(TH1F const* histo, int ix, int jx, bool useWidth, double* error);
template double getHistogramIntegralAndError<TH2F>(TH2F const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error);
template double getHistogramIntegralAndError<TH3F>(TH3F const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error);


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

template<> void wipeOverUnderFlows<TH1F>(TH1F* hwipe, bool rescale){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    hwipe->SetBinContent(binx, 0);
    hwipe->SetBinError(binx, 0);
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  if (rescale) hwipe->Scale(wipeScale);
}
template<> void wipeOverUnderFlows<TH2F>(TH2F* hwipe, bool rescale){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (
        (binx>=1 && binx<=hwipe->GetNbinsX())
        &&
        (biny>=1 && biny<=hwipe->GetNbinsY())
        ) continue;
      hwipe->SetBinContent(binx, biny, 0);
      hwipe->SetBinError(binx, biny, 0);
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  if (rescale) hwipe->Scale(wipeScale);
}
template<> void wipeOverUnderFlows<TH3F>(TH3F* hwipe, bool rescale){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1, 0, hwipe->GetNbinsZ()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      for (int binz=0; binz<=hwipe->GetNbinsZ()+1; binz++){
        if (
          (binx>=1 && binx<=hwipe->GetNbinsX())
          &&
          (biny>=1 && biny<=hwipe->GetNbinsY())
          &&
          (binz>=1 && binz<=hwipe->GetNbinsZ())
          ) continue;
        hwipe->SetBinContent(binx, biny, binz, 0);
        hwipe->SetBinError(binx, biny, binz, 0);
      }
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  if (rescale) hwipe->Scale(wipeScale);
}

float calculateSimpleProductError(
  float const v1, float const e1, float const p1,
  float const v2, float const e2, float const p2
){
  assert(e1>=0. && e2>=0.);
  float d=0;
  float val=0;
  if (
    ((p1<0. && v1!=0.) || p1>=0.)
    &&
    ((p2<0. && v2!=0.) || p2>=0.)
    ){
    val=pow(v1, p1)*pow(v2, p2);
    if (v1!=0.) d += pow(p1*e1/v1, 2);
    if (v2!=0.) d += pow(p2*e2/v2, 2);
  }
  d=std::abs(val)*sqrt(d);
  return d;
}
float calculateEfficiencyError(
  float const sumW, float const sumWAll,
  float const sumWsq, float const sumWsqAll
){
  float const& sumWp=sumW;
  float const& sumWsqp=sumWsq;
  float const sumWm = sumWAll-sumWp;
  float const sumWsqm = sumWsqAll-sumWsqp;
  float numerator, denominator;
  float ratio=0;
  if (sumWAll!=0.){
    numerator = sqrt(std::max(float(0), float(sumWsqp*pow(sumWm, 2) + sumWsqm*pow(sumWp, 2))));
    denominator = pow(sumWAll, 2);
    ratio = numerator/denominator;
  }
  return ratio;
}
float translateEfficiencyErrorToNumeratorError(
  float const eff, float const sumWAll,
  float const effErr, float const sumWsqAll
){
  float numerator = pow(effErr*pow(sumWAll, 2), 2);
  float sumWp = eff*sumWAll;
  float sumWm = sumWAll-sumWp;
  float res = std::max(float(0), float((numerator - sumWsqAll*pow(sumWp, 2))/(pow(sumWm, 2)-pow(sumWp, 2))));
  return sqrt(res);
}

template<> void conditionalizeHistogram<TH2F>(TH2F* histo, unsigned int iaxis, std::vector<std::pair<TH2F*, float>> const* conditionalsReference, bool useWidth, bool useEffErr){
  const bool forceUseSimpleErr = (conditionalsReference || !useEffErr);
  TAxis* axis[2]={ nullptr };
  switch (iaxis){
  case 0:
    axis[0]=histo->GetXaxis();
    axis[1]=histo->GetYaxis();
    break;
  case 1:
    axis[0]=histo->GetYaxis();
    axis[1]=histo->GetXaxis();
    break;
  default:
    return;
  }
  int nbins[2]; for (unsigned int i=0; i<2; i++) nbins[i]=axis[i]->GetNbins();

  for (int i=0; i<=nbins[0]+1; i++){
    double integral=1;
    double integralerror=0;

    int int_xb[2]={ 0 }, int_yb[2]={ 0 };
    switch (iaxis){
    case 0:
      int_xb[0]=i;
      int_xb[1]=i;
      int_yb[0]=0;
      int_yb[1]=nbins[1]+1;
      break;
    case 1:
      int_yb[0]=i;
      int_yb[1]=i;
      int_xb[0]=0;
      int_xb[1]=nbins[1]+1;
      break;
    }

    if (!conditionalsReference) integral = getHistogramIntegralAndError<TH2F>(histo, int_xb[0], int_xb[1], int_yb[0], int_yb[1], useWidth, &integralerror);
    else{
      for (std::pair<TH2F*, float> const& hh:(*conditionalsReference)){
        double extraintegralerror=0;
        double extraintegral = getHistogramIntegralAndError<TH2F>(hh.first, int_xb[0], int_xb[1], int_yb[0], int_yb[1], useWidth, &extraintegralerror);
        integralerror = calculateSimpleProductError(extraintegral, extraintegralerror, hh.second, integral, integralerror, 1);
        integral *= pow(extraintegral, hh.second);
      }
    }
    for (int j=0; j<=nbins[1]+1; j++){
      int ix=0, iy=0;
      switch (iaxis){
      case 0:
        ix=i;
        iy=j;
        break;
      case 1:
        iy=i;
        ix=j;
        break;
      }

      double width = 1;
      double binerror;
      double bincontent = getHistogramIntegralAndError<TH2F>(histo, ix, ix, iy, iy, useWidth, &binerror);

      if (useWidth && j>=1 && j<=nbins[1]) width *= axis[1]->GetBinWidth(j);

      double hval=0;
      double herr=0;
      if (integral!=0.){
        hval = bincontent/integral;
        if (!forceUseSimpleErr) herr = calculateEfficiencyError(bincontent, integral, pow(binerror, 2), pow(integralerror, 2));
        else herr = calculateSimpleProductError(bincontent, binerror, 1, integral, integralerror, -1);
        hval /= width; herr /= width;
      }

      histo->SetBinContent(ix, iy, hval);
      histo->SetBinError(ix, iy, herr);
    }
  }
}
template<> void conditionalizeHistogram<TH3F>(TH3F* histo, unsigned int iaxis, std::vector<std::pair<TH3F*, float>> const* conditionalsReference, bool useWidth, bool useEffErr){
  const bool forceUseSimpleErr = (conditionalsReference || !useEffErr);
  TAxis* axis[3]={ nullptr };
  switch (iaxis){
  case 0:
    axis[0]=histo->GetXaxis();
    axis[1]=histo->GetYaxis();
    axis[2]=histo->GetZaxis();
    break;
  case 1:
    axis[0]=histo->GetYaxis();
    axis[1]=histo->GetZaxis();
    axis[2]=histo->GetXaxis();
    break;
  case 2:
    axis[0]=histo->GetZaxis();
    axis[1]=histo->GetXaxis();
    axis[2]=histo->GetYaxis();
    break;
  default:
    return;
  }
  int nbins[3]; for (unsigned int i=0; i<3; i++) nbins[i]=axis[i]->GetNbins();

  for (int i=0; i<=nbins[0]+1; i++){
    double integral=1;
    double integralerror=0;

    int int_xb[2]={ 0 }, int_yb[2]={ 0 }, int_zb[2]={ 0 };
    switch (iaxis){
    case 0:
      int_xb[0]=i;
      int_xb[1]=i;
      int_yb[0]=0;
      int_yb[1]=nbins[1]+1;
      int_zb[0]=0;
      int_zb[1]=nbins[2]+1;
      break;
    case 1:
      int_yb[0]=i;
      int_yb[1]=i;
      int_zb[0]=0;
      int_zb[1]=nbins[1]+1;
      int_xb[0]=0;
      int_xb[1]=nbins[2]+1;
      break;
    case 2:
      int_zb[0]=i;
      int_zb[1]=i;
      int_xb[0]=0;
      int_xb[1]=nbins[1]+1;
      int_yb[0]=0;
      int_yb[1]=nbins[2]+1;
      break;
    }

    if (!conditionalsReference) integral = getHistogramIntegralAndError<TH3F>(histo, int_xb[0], int_xb[1], int_yb[0], int_yb[1], int_zb[0], int_zb[1], useWidth, &integralerror);
    else{
      for (std::pair<TH3F*, float> const& hh:(*conditionalsReference)){
        double extraintegralerror=0;
        double extraintegral = getHistogramIntegralAndError<TH3F>(hh.first, int_xb[0], int_xb[1], int_yb[0], int_yb[1], int_zb[0], int_zb[1], useWidth, &extraintegralerror);
        integralerror = calculateSimpleProductError(extraintegral, extraintegralerror, hh.second, integral, integralerror, 1);
        integral *= pow(extraintegral, hh.second);
      }
    }
    for (int j=0; j<=nbins[1]+1; j++){
      for (int k=0; k<=nbins[2]+1; k++){
        int ix=0, iy=0, iz=0;
        switch (iaxis){
        case 0:
          ix=i;
          iy=j;
          iz=k;
          break;
        case 1:
          ix=k;
          iy=i;
          iz=j;
          break;
        case 2:
          ix=j;
          iy=k;
          iz=i;
          break;
        }

        double width = 1;
        double binerror;
        double bincontent = getHistogramIntegralAndError<TH3F>(histo, ix, ix, iy, iy, iz, iz, useWidth, &binerror);

        if (useWidth && j>=1 && j<=nbins[1]) width *= axis[1]->GetBinWidth(j);
        if (useWidth && k>=1 && k<=nbins[2]) width *= axis[2]->GetBinWidth(k);

        double hval=0;
        double herr=0;
        if (integral!=0.){
          hval = bincontent/integral;
          if (!forceUseSimpleErr) herr = calculateEfficiencyError(bincontent, integral, pow(binerror, 2), pow(integralerror, 2));
          else herr = calculateSimpleProductError(bincontent, binerror, 1, integral, integralerror, -1);
          hval /= width; herr /= width;
        }

        histo->SetBinContent(ix, iy, iz, hval);
        histo->SetBinError(ix, iy, iz, herr);
      }
    }
  }
}

template<> void multiplyHistograms<TH1F>(TH1F const* h1, TH1F const* h2, TH1F*& hAssign, bool useEffErr){
  if (h1->GetNbinsX()!=h2->GetNbinsX()) return;
  const int nbinsx = h1->GetNbinsX();
  for (int binx=0; binx<=nbinsx+1; binx++){
    float sumW = h1->GetBinContent(binx);
    float sumWAll = h2->GetBinContent(binx);
    float sumWsq = pow(h1->GetBinError(binx), 2);
    float sumWsqAll = pow(h2->GetBinError(binx), 2);
    float bincontent=sumW*sumWAll;
    float binerror=0;
    if (useEffErr) binerror = translateEfficiencyErrorToNumeratorError(sumW, sumWAll, sumWsq, sumWsqAll);
    else binerror = calculateSimpleProductError(sumW, sqrt(sumWsq), 1, sumWAll, sqrt(sumWsqAll), 1);
    if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
      bincontent=0;
      binerror=0;
    }
    hAssign->SetBinContent(binx, bincontent);
    hAssign->SetBinError(binx, binerror);
  }
}
template<> void multiplyHistograms<TH2F>(TH2F const* h1, TH2F const* h2, TH2F*& hAssign, bool useEffErr){
  if (h1->GetNbinsX()!=h2->GetNbinsX()) return;
  const int nbinsx = h1->GetNbinsX();
  if (h1->GetNbinsY()!=h2->GetNbinsY()) return;
  const int nbinsy = h1->GetNbinsY();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      float sumW = h1->GetBinContent(binx, biny);
      float sumWAll = h2->GetBinContent(binx, biny);
      float sumWsq = pow(h1->GetBinError(binx, biny), 2);
      float sumWsqAll = pow(h2->GetBinError(binx, biny), 2);
      float bincontent=sumW*sumWAll;
      float binerror=0;
      if (useEffErr) binerror = translateEfficiencyErrorToNumeratorError(sumW, sumWAll, sumWsq, sumWsqAll);
      else binerror = calculateSimpleProductError(sumW, sqrt(sumWsq), 1, sumWAll, sqrt(sumWsqAll), 1);
      if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
        bincontent=0;
        binerror=0;
      }
      hAssign->SetBinContent(binx, biny, bincontent);
      hAssign->SetBinError(binx, biny, binerror);
    }
  }
}
template<> void multiplyHistograms<TH3F>(TH3F const* h1, TH3F const* h2, TH3F*& hAssign, bool useEffErr){
  if (h1->GetNbinsX()!=h2->GetNbinsX()) return;
  const int nbinsx = h1->GetNbinsX();
  if (h1->GetNbinsY()!=h2->GetNbinsY()) return;
  const int nbinsy = h1->GetNbinsY();
  if (h1->GetNbinsZ()!=h2->GetNbinsZ()) return;
  const int nbinsz = h1->GetNbinsZ();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      for (int binz=0; binz<=nbinsz+1; binz++){
        float sumW = h1->GetBinContent(binx, biny, binz);
        float sumWAll = h2->GetBinContent(binx, biny, binz);
        float sumWsq = pow(h1->GetBinError(binx, biny, binz), 2);
        float sumWsqAll = pow(h2->GetBinError(binx, biny, binz), 2);
        float bincontent=sumW*sumWAll;
        float binerror=0;
        if (useEffErr) binerror = translateEfficiencyErrorToNumeratorError(sumW, sumWAll, sumWsq, sumWsqAll);
        else binerror = calculateSimpleProductError(sumW, sqrt(sumWsq), 1, sumWAll, sqrt(sumWsqAll), 1);
        if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
          bincontent=0;
          binerror=0;
        }
        hAssign->SetBinContent(binx, biny, binz, bincontent);
        hAssign->SetBinError(binx, biny, binz, binerror);
      }
    }
  }
}
template<> void multiplyHistograms<TH2F>(TH2F const* h1, TH1F const* h2, unsigned int matchDimension, TH2F*& hAssign, bool useEffErr){
  if (matchDimension>1) assert(0);
  if (matchDimension==0 && h1->GetNbinsX()!=h2->GetNbinsX()) return;
  const int nbinsx = h1->GetNbinsX();
  if (matchDimension==1 && h1->GetNbinsY()!=h2->GetNbinsX()) return;
  const int nbinsy = h1->GetNbinsY();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      float sumW = h1->GetBinContent(binx, biny);
      float sumWsq = pow(h1->GetBinError(binx, biny), 2);
      float sumWAll=0, sumWsqAll=0;
      switch (matchDimension){
      case 0:
        sumWAll = h2->GetBinContent(binx);
        sumWsqAll = pow(h2->GetBinError(binx), 2);
        break;
      case 1:
        sumWAll = h2->GetBinContent(biny);
        sumWsqAll = pow(h2->GetBinError(biny), 2);
        break;
      }
      float bincontent=sumW*sumWAll;
      float binerror=0;
      if (useEffErr) binerror = translateEfficiencyErrorToNumeratorError(sumW, sumWAll, sumWsq, sumWsqAll);
      else binerror = calculateSimpleProductError(sumW, sqrt(sumWsq), 1, sumWAll, sqrt(sumWsqAll), 1);
      if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
        bincontent=0;
        binerror=0;
      }
      hAssign->SetBinContent(binx, biny, bincontent);
      hAssign->SetBinError(binx, biny, binerror);
    }
  }
}
template<> void multiplyHistograms<TH3F>(TH3F const* h1, TH1F const* h2, unsigned int matchDimension, TH3F*& hAssign, bool useEffErr){
  if (matchDimension>2) assert(0);
  if (matchDimension==0 && h1->GetNbinsX()!=h2->GetNbinsX()) return;
  const int nbinsx = h1->GetNbinsX();
  if (matchDimension==1 && h1->GetNbinsY()!=h2->GetNbinsX()) return;
  const int nbinsy = h1->GetNbinsY();
  if (matchDimension==2 && h1->GetNbinsZ()!=h2->GetNbinsX()) return;
  const int nbinsz = h1->GetNbinsZ();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      for (int binz=0; binz<=nbinsz+1; binz++){
        float sumW = h1->GetBinContent(binx, biny, binz);
        float sumWsq = pow(h1->GetBinError(binx, biny, binz), 2);
        float sumWAll=0, sumWsqAll=0;
        switch (matchDimension){
        case 0:
          sumWAll = h2->GetBinContent(binx);
          sumWsqAll = pow(h2->GetBinError(binx), 2);
          break;
        case 1:
          sumWAll = h2->GetBinContent(biny);
          sumWsqAll = pow(h2->GetBinError(biny), 2);
          break;
        case 2:
          sumWAll = h2->GetBinContent(binz);
          sumWsqAll = pow(h2->GetBinError(binz), 2);
          break;
        }
        float bincontent=sumW*sumWAll;
        float binerror=0;
        if (useEffErr) binerror = translateEfficiencyErrorToNumeratorError(sumW, sumWAll, sumWsq, sumWsqAll);
        else binerror = calculateSimpleProductError(sumW, sqrt(sumWsq), 1, sumWAll, sqrt(sumWsqAll), 1);
        if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
          bincontent=0;
          binerror=0;
        }
        hAssign->SetBinContent(binx, biny, binz, bincontent);
        hAssign->SetBinError(binx, biny, binz, binerror);
      }
    }
  }
}
