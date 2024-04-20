#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <unistd.h>
#include <regex>
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
  TString name;
  double rate;
  unordered_map<TString, pair<TString, TString>> systematics;
  unordered_map<TString, vector<TH1F*>> syst_templates_map;

  process_spec() : rate(0){}
  process_spec(TString name_, double rate_) : name(name_), rate(rate_){}
  process_spec(const process_spec& other) : name(other.name), rate(other.rate), systematics(other.systematics), syst_templates_map(other.syst_templates_map){}
  virtual ~process_spec(){}

  void setSystematic(TString systname, TString systtype, TString systline);
  void acquireTemplates(TFile* finput);

  void clean();

};



template<typename T, typename U> bool replaceString(T& strinput, U strTakeOut, U strPutIn);
template<> bool replaceString<TString, const TString>(TString& strinput, const TString strTakeOut, const TString strPutIn){
  Ssiz_t ipos=strinput.Index(strTakeOut);
  if (ipos!=-1){ strinput.Replace(ipos, strTakeOut.Length(), strPutIn); return true; }
  else return false;
}
template<> bool replaceString<TString, const char*>(TString& strinput, const char* strTakeOut, const char* strPutIn){
  Ssiz_t ipos=strinput.Index(strTakeOut);
  if (ipos!=-1){ strinput.Replace(ipos, strlen(strTakeOut), strPutIn); return true; }
  else return false;
}
template<> bool replaceString<std::string, const std::string>(std::string& strinput, const std::string strTakeOut, const std::string strPutIn){
  std::string::size_type ipos=strinput.find(strTakeOut);
  if (ipos!=std::string::npos){ strinput.replace(ipos, strTakeOut.length(), strPutIn); return true; }
  else return false;
}
template<> bool replaceString<std::string, const char*>(std::string& strinput, const char* strTakeOut, const char* strPutIn){
  std::string::size_type ipos=strinput.find(strTakeOut);
  if (ipos!=std::string::npos){ strinput.replace(ipos, strlen(strTakeOut), strPutIn); return true; }
  else return false;
}


TString get_process_finalname(TString strproc){
  if (strproc.BeginsWith("offggH") || strproc == "back_ggZZ") return "ggZZ_offshell";
  else if (strproc.BeginsWith("offqqH") || strproc == "back_VVZZ") return "VVZZ_offshell";
  else if (strproc == "qqH_0PM") return "VVZZ_onshell";
  else if (strproc == "ttH_0PM") return "ttH";
  else if (strproc == "back_zjets") return "Zjets";
  else if (strproc == "back_qqZZ") return "qqZZ";
  else{
    cerr << "Process " << strproc << " is not recognized." << endl;
    assert(false);
  }
  return "";
}
TString get_template_corename(TString procname, unsigned int i){
  if (procname=="ttH") return "ttH_0PM";
  else if (procname=="VVZZ_onshell") return "qqH_0PM";
  else if (procname=="Zjets") return "back_zjets";
  else if (procname=="qqZZ") return "back_qqZZ";
  else if (procname=="ggZZ_offshell"){
    TString res = "offggH";
    if (i==0) res += "_0PM";
    else if (i==1) res += "_g11g21";
    else res = "back_ggZZ";
    return res;
  }
  else if (procname=="VVZZ_offshell"){
    TString res = "offqqH";
    if (i==0) res += "_0PM";
    else if (i==1) res += "_g12g22";
    else res = "back_VVZZ";
    return res;
  }
  else{
    cerr << "Process " << procname << " is not recognized." << endl;
    assert(false);
  }
  return "";
}


void process_spec::setSystematic(TString systname, TString systtype, TString systline){
  systematics[systname] = pair<TString, TString>(systtype, systline);
  if (systtype.Index("shape")>=0){
    syst_templates_map[systname+"Down"] = vector<TH1F*>();
    syst_templates_map[systname+"Up"] = vector<TH1F*>();
  }
}
void process_spec::acquireTemplates(TFile* finput){
  constexpr double scale = 1;
  constexpr bool useOnlySM = true;

  bool const is_offshell = name.Contains("offshell");
  unsigned int const ntpls = (is_offshell ? 3 : 1);

  std::vector<TString> tplinputcorenames;
  std::vector<TString> tploutputcorenames;
  for (unsigned int i=0; i<ntpls; i++){
    TString intplname = get_template_corename(name, i);
    if (is_offshell){
      tplinputcorenames.push_back(intplname);

      TString nname;
      if (i==0) nname = Form("T_%s_Sig", name.Data());
      else if (i==1) nname = Form("T_%s_Int_Re", name.Data());
      else nname = Form("T_%s_Bkg", name.Data());
      tploutputcorenames.push_back(nname);
    }
    else{
      tplinputcorenames.push_back(get_template_corename(name, 0));
      TString nname = Form("T_%s", name.Data());
      tploutputcorenames.push_back(nname);
    }
  }

  std::vector<TString> const strappend_Nominal{ "", "_positive", "_negative" };

  syst_templates_map["Nominal"] = vector<TH1F*>();
  for (auto& pp:syst_templates_map){
    auto const& syst = pp.first;
    std::vector<TString> strappend(strappend_Nominal);
    if (syst!="Nominal"){ for (auto& ss:strappend) ss = ss + "_" + syst; }
    for (unsigned int itpl=0; itpl<ntpls; itpl++){
      TH1F* htpl = nullptr;
      TString strouttplname = tploutputcorenames.at(itpl);

      TDirectory* tmpdir = gDirectory;
      finput->cd();
      for (unsigned int isapp=0; isapp<strappend.size(); isapp++){
        if (itpl!=1 && isapp>0) break;
        if (itpl==1 && isapp==0) continue;

        auto const& ss = strappend.at(isapp);
        TString strintplname = tplinputcorenames.at(itpl);
        strintplname += ss;
        TH1F* htmp = dynamic_cast<TH1F*>(finput->Get(strintplname));
        if (!htmp){
          cerr << "Could not acquire template " << strintplname << ". Trying the nominal variation..." << endl;
          strintplname = tplinputcorenames.at(itpl);
          strintplname += strappend_Nominal.at(isapp);
          htmp = dynamic_cast<TH1F*>(finput->Get(strintplname));
          if (!htmp) cerr << "Could not acquire template " << strintplname << " entirely!" << endl;
        }
        if (!htmp){
          ++isapp;
          continue;
        }

        tmpdir->cd();
        htmp = (TH1F*) htmp->Clone(strintplname + "_clone");
        if (strintplname.Contains("negative")) htmp->Scale(-1.);

        if (!htpl){
          htpl = (TH1F*) htmp->Clone(strouttplname);
        }
        else htpl->Add(htmp);

        cout << "[" << name << "][" << syst << "]: Added template " << strintplname << " to " << strouttplname << " (int=" << htmp->Integral("width") << ")." << endl;

        delete htmp;
        finput->cd();
      }

      if (!htpl){
        cerr << "Could not acquire template " << tploutputcorenames.at(itpl) << " from the set " << tplinputcorenames.at(itpl) << endl;
        exit(1);
      }

      pp.second.push_back(htpl);

      tmpdir->cd();
    }
  }
}
void process_spec::clean(){
  for (auto& pp:syst_templates_map){
    for (auto& hh:pp.second) delete hh;
  }
}



void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
void splitOption(const TString rawoption, TString& wish, TString& value, char delimiter);
void splitOptionRecursive(const TString rawoption, vector<TString>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);
Bool_t checkListVariable(const vector<TString>& list, const TString& var);

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

void getSqrtsPeriod(TString& strSqrtsPeriod, TString& strSqrts, TString& strPeriod){
  strSqrts = "13TeV";
  strPeriod = "Run2";
  strSqrtsPeriod = strSqrts + '_' + strPeriod;
}
TString getSystRename(TString const& systname, TString const& systLine, TString const& strSqrts, TString const& strPeriod, TString const& strCategory, TString const& strChannel){
  TString res = systname;
  if (res == "BRhiggs_hzz4l" || res == "hzz_br" || res == "BR_hzz") res = "BRhiggs_hzz";

  else if (res == "zjets_4e") res = "CMS_hzz4l_4e_zjets_13TeV";
  else if (res == "zjets_4mu") res = "CMS_hzz4l_4mu_zjets_13TeV";
  else if (res == "zjets_2e2mu") res = "CMS_hzz4l_2e2mu_zjets_13TeV";
  else if (res == "zjets_2mu2e") res = "CMS_hzz4l_2e2mu_zjets_13TeV";

  else if (res == "CMS_pythia_scale") res = "CMS_scale_pythia";
  else if (res == "CMS_pythia_tune") res = "CMS_tune_pythia";

  else if (res == "CMS16_btag_comb") res = "CMS_btag_comb_13TeV_2016";
  else if (res == "CMS17_btag_comb") res = "CMS_btag_comb_13TeV_2017";
  else if (res == "CMS18_btag_comb") res = "CMS_btag_comb_13TeV_2018";

  else if (res == "CMS16_scale_j") res = "CMS_scale_j_13TeV_2016";
  else if (res == "CMS17_scale_j") res = "CMS_scale_j_13TeV_2017";
  else if (res == "CMS18_scale_j") res = "CMS_scale_j_13TeV_2018";

  else if (res == "CMS16_res_j") res = "CMS_res_j_13TeV_2016";
  else if (res == "CMS17_res_j") res = "CMS_res_j_13TeV_2017";
  else if (res == "CMS18_res_j") res = "CMS_res_j_13TeV_2018";

  else if (res.Contains("QCDscale_muR")) replaceString<TString, TString const>(res, "QCDscale_muR", "QCDscale_ren");
  else if (res.Contains("QCDscale_muF")) replaceString<TString, TString const>(res, "QCDscale_muF", "QCDscale_fac");

  else if (res == "pdf_Higgs_gg") res = "pdf_variation_Higgs_gg";
  else if (res == "pdf_As_Higgs_gg") res = "pdf_asmz_Higgs_gg";
  else if (res == "pdf_Higgs_qq" || res == "pdf_Higgs_qqbar") res = "pdf_variation_Higgs_qqbar";
  else if (res == "pdf_As_Higgs_qq" || res == "pdf_As_Higgs_qqbar") res = "pdf_asmz_Higgs_qqbar";

  else if (res == "pdf_qq" || res == "pdf_qqbar") res = "pdf_variation_qqbar";
  else if (res == "pdf_As_qq" || res == "pdf_As_qqbar") res = "pdf_asmz_qqbar";

  else if (res == "EWcorr_qqZZ") res = "EWcorr_VV";

  else if (res == "CMS_eff_e") res = "CMS_eff_syst_e_13TeV";
  else if (res == "CMS_eff_mu" || res == "CMS_eff_m") res = "CMS_eff_syst_m_13TeV";

  else if (res == "QCDmodl_qqZZ0isr") res = "QCDscale_VV0in";
  else if (res == "QCDmodl_qqZZ1isr") res = "QCDscale_VV1in";
  else if (res == "QCDmodl_qqZZ2isr") res = "QCDscale_VV2in";

  return res;
}

bool checkProcessIsBkg(TString procname){
  procname.ToLower();
  return (procname.Contains("back_qqZZ") || procname.Contains("back_zjets"));
}

void getDataTree(TFile* finput, TString coutput){
  TH1F* data = dynamic_cast<TH1F*>(finput->Get("data_obs"));
  if (!data) return;
  TFile* foutput = TFile::Open(coutput, "recreate");
  foutput->WriteTObject(data);
  foutput->Close();
}

double getBestLumiCurrent(){ return 36.326450 + 41.529152052 + 59.740565209; }

void getCatRename(string& catname){
  if (catname == "VBFtagged") catname = "JJVBFTagged";
  else if (catname == "VHtagged") catname = "HadVHTagged";
}

// flip_asym_lnN: Fixes the bug in lnN systematics to be ordered as "down/up" instead of "up/down" per Combine conventions:
// https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part2/settinguptheanalysis/#a-simple-counting-experiment.
void getTemplates_21019(TString cinput_txt, TString coutput_main, bool flip_asym_lnN = true){
  TString cinput_main;
  {
    vector<TString> splitinput;
    splitOptionRecursive(cinput_txt, splitinput, '/');
    for (unsigned int is=0; is<splitinput.size()-1; ++is) cinput_main += splitinput.at(is) + "/";
    if (cinput_txt.BeginsWith("/")) cinput_main = TString("/") + cinput_main;
    replaceString<TString, TString const>(cinput_txt, cinput_main, "");
  }
  cout << "cinput_main: " << cinput_main << endl;
  cout << "cinput_txt: " << cinput_txt << endl;

  TString cinput = cinput_txt; cinput = cinput(0, cinput.First('.'));
  cout << "cinput: " << cinput << endl;

  TString strSqrtsPeriod, strSqrts, strPeriod;
  getSqrtsPeriod(strSqrtsPeriod, strSqrts, strPeriod);

  string strinput = cinput.Data();
  string channame = "4l";
  string catname = strinput.data();
  getCatRename(catname);

  cout << "Processing channel " << channame << " and category " << catname << "..." << endl;

  TString coutput_templates = "Decompilation/Templates/" + coutput_main;
  TString coutput_inputs = "Decompilation/Inputs/" + coutput_main;
  TString coutput_data = "Decompilation/Data/" + coutput_main;

  gSystem->Exec("mkdir -p " + coutput_templates);
  gSystem->Exec("mkdir -p " + coutput_inputs);
  gSystem->Exec("mkdir -p " + coutput_data);

  TString strChannel = channame.c_str();
  TString strCategory = catname.c_str();

  TString coutput_txt = Form("%s/inputs_%s_%s%s", coutput_inputs.Data(), channame.c_str(), catname.c_str(), ".txt");
  ofstream tout(coutput_txt.Data());

  TFile* finput = TFile::Open(cinput_main + cinput + ".input.root", "read");

  TString coutput_data_root = Form("%s/hto%s_%s_%s.root", coutput_data.Data(), channame.c_str(), catname.c_str(), strSqrtsPeriod.Data());
  getDataTree(finput, coutput_data_root);

  ifstream tin((cinput_main + cinput_txt).Data());
  string line;

  // Get process names
  while (line.find("process")==string::npos) std::getline(tin, line);
  vector<TString> inprocnames;
  vector<TString> procnames;
  {
    std::vector<std::string> splitline;
    splitOptionRecursive(line, splitline, ' ');
    for (auto const& sproc:splitline){
      if (sproc=="" || sproc=="process") continue;
      TString strproc = sproc.c_str();
      inprocnames.push_back(strproc);
      TString proc_fname = get_process_finalname(strproc);
      if (!checkListVariable(procnames, proc_fname)) procnames.push_back(proc_fname);
    }
  }

  // Get process rates
  vector<double> procrates;
  // No need for process rates, templates contain them already.
  procrates.assign(procnames.size(), 1);

  // Get process pdfs
  unordered_map<string, process_spec> procSpecs;
  for (unsigned int ip=0; ip<procnames.size(); ip++){
    auto const& procname = procnames.at(ip);
    auto const& procrate = procrates.at(ip);
    cout << procnames.at(ip) << ": " << procrates.at(ip) << endl;

    procSpecs[procname.Data()] = process_spec(procname, procrate);
  }

  // Get systematics
  std::vector<string> trueParamSyst;
  unordered_map<string, string> tplSyst;
  unordered_map<string, string> logSyst;
  unordered_map<string, string> paramSyst;
  while (!tin.eof()){
    std::getline(tin, line);

    string strline = line;
    std::replace(strline.begin(), strline.end(), ',', ' ');
    strline.erase(std::remove(strline.begin(), strline.end(), '['), strline.end());
    strline.erase(std::remove(strline.begin(), strline.end(), ']'), strline.end());
    bool isMaybeShape = strline.find("shape1?")!=string::npos;
    bool isShape = strline.find("shape1")!=string::npos;
    bool isParam = strline.find("param")!=string::npos;
    bool isLog = strline.find("lnN")!=string::npos || isMaybeShape;
    if (isShape || isLog || isParam){
      vector<string> systdist;
      splitOptionRecursive(strline, systdist, ' ');
      string systname = systdist.at(0);
      string systtype = systdist.at(1);
      string accumulate="";
      cout << "Processing systematic " << systname << " [" << systtype << "]..." << endl;

      if (isMaybeShape){
        bool has_log = false;
        for (unsigned int ip=0; ip<inprocnames.size(); ip++){
          string systline = systdist.at(ip+2);
          if (systline.find("-")==string::npos && systline!=""){
            if (systline.find("/")!=string::npos || systline!="1"){
              has_log = true;
              break;
            }
          }
        }
        if (has_log){
          cout << "\t- Determined that " << systname << " is a lnN systematic." << endl;
          isLog = true;
          isShape = false;
          systtype = "lnN";
        }
        else{
          cout << "\t- Determined that " << systname << " is a shape systematic." << endl;
          isShape = true;
          isLog = false;
          systtype = "shape1";
        }
        isMaybeShape = false;
      }

      if (systname == "kf_ggZZ_back"){
        cout << "Skipping " << systname << "..." << endl;
        continue;
      }

      bool isTrueParam = false;
      isTrueParam |= (systname == "kbkg_gg");

      if (isTrueParam){
        trueParamSyst.push_back(systname);
        // Directly add since this is a special systematic
        for (unsigned int ip=2; ip<systdist.size(); ip++){
          if (ip>2) accumulate += ":";
          accumulate += systdist.at(ip);
        }
        paramSyst[systname] = accumulate;
      }
      else if (isShape || isLog){
        std::vector<TString> found_procnames;
        for (unsigned int ip=0; ip<inprocnames.size(); ip++){
          auto const& strproc = inprocnames.at(ip);
          TString procname = get_process_finalname(strproc);

          string systline = systdist.at(ip+2);
          bool const isNonEmpty = (systline.find("-")==string::npos && systline!="");

          if (isNonEmpty && !checkListVariable(found_procnames, procname)) found_procnames.push_back(procname);
          else continue;

          if (isNonEmpty){
            if (flip_asym_lnN && isLog && systline.find("/")!=string::npos){
              vector<string> splitsystline;
              splitOptionRecursive(systline, splitsystline, '/');
              assert(splitsystline.size()==2);
              systline = splitsystline.at(1) + "/" + splitsystline.at(0);
            }
            std::replace(systline.begin(), systline.end(), '/', ':');
            procSpecs[procname.Data()].setSystematic(systname, systtype, systline);
            if (isShape) accumulate += string(procname.Data()) + ":0:1 ";
            else accumulate += string(procname.Data()) + ":" + systline + " ";
          }
        }
        if (isLog) logSyst[systname] = accumulate;
        else if (isShape) tplSyst[systname] = accumulate;
      }
    }
  }
  cout << "Processes: ";
  for (auto const& pname:procnames) cout << pname << " ";
  cout << endl;
  cout << "Process specification keys: ";
  for (auto it=procSpecs.begin(); it!=procSpecs.end(); it++) cout << it->first << " ";
  cout << endl;


  // Write input file
  if (strSqrts!=""){
    TString strSqrtsBare=strSqrts;
    int ipos = strSqrtsBare.Index("TeV");
    if (ipos>=0) strSqrtsBare.Resize(ipos);
    tout << "sqrts " << strSqrtsBare << endl;
  }
  if (strPeriod!="") tout << "period " << strPeriod << endl;
  tout << "decay " << channame << endl;
  tout << "lumi " << std::setprecision(11) << getBestLumiCurrent() << endl;
  tout << "category " << catname << endl;

  // Write channels
  for (unsigned int ip=0; ip<procnames.size(); ip++){
    auto const pname = procSpecs[procnames.at(ip).Data()].name;
    unsigned int proccode = (checkProcessIsBkg(pname) ? 0 : (pname.Contains("offshell") ? 2 : 1));
    tout << "channel " << procnames.at(ip) << " 1 -1 " << proccode;
    tout << endl;
  }

  // Write systematics
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
      bool const isTrueParamType = (std::find(trueParamSyst.begin(), trueParamSyst.end(), syst->first)!=trueParamSyst.end());
      if (isTrueParamType) tout << "systematic " << systName << " param " << syst->second << endl;
      else tout << "systematic " << systName << " template " << syst->second << endl;
    }
  }
  // Write the kbkg_gg systematic as well if it is not found.
  if (!checkListVariable(trueParamSyst, "kbkg_gg")) tout << "systematic kbkg_gg param 1:0.1:0:2" << endl;


  for (unsigned int ip=0; ip<procnames.size(); ip++){
    cout << "Attempting to extract tpls for process " << procnames.at(ip) << endl;

    procSpecs[procnames.at(ip).Data()].acquireTemplates(finput);
    for (auto const& pp:procSpecs[procnames.at(ip).Data()].syst_templates_map){
      auto syst = pp.first;
      auto const& htpls = pp.second;
      TString systNameCoreOld = syst;
      replaceString<TString, TString const>(systNameCoreOld, "Up", "");
      replaceString<TString, TString const>(systNameCoreOld, "Down", "");
      TString systNameCoreNew = getSystRename(systNameCoreOld, "", strSqrts, strPeriod, strCategory, strChannel);
      replaceString<TString, TString const>(syst, systNameCoreOld, systNameCoreNew);
      TString coutput_root = Form("%s/Hto%s_%s_FinalTemplates_%s_%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procnames.at(ip).Data(), syst.Data(), ".root");
      TFile* foutput = TFile::Open(coutput_root, "recreate");
      cout << "[" << procnames.at(ip) << "][" << syst << "] integrals:";
      for (auto const& hh:htpls){
        cout << " " << hh->Integral("width");
        foutput->WriteTObject(hh);
      }
      cout << endl;
      foutput->Close();
    }

    procSpecs[procnames.at(ip).Data()].clean();
  }

  tin.close();
  finput->Close();
  tout.close();
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
    if (result=="" && suboption.find(delimiter)!=std::string::npos) result=suboption; // This can happen if the string starts with the delimiter.
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

Bool_t checkListVariable(const vector<string>& list, const string& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
}
Bool_t checkListVariable(const vector<TString>& list, const TString& var){
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
