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
#include "RooCmdArg.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
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
  RooAbsReal* norm;
  double rate;
  TString name;

  process_spec() : pdf(0), norm(0), rate(0){}
  process_spec(TString name_, double rate_) : pdf(nullptr), norm(nullptr), rate(rate_), name(name_){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()){}
  process_spec(const process_spec& other) : pdf(other.pdf), norm(other.norm), rate(other.rate), name(other.name){}
};

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


unsigned int extractTemplates(process_spec& proc, RooDataSet* data, unordered_map<TString, TH2F>& procshape_2D, unordered_map<TString, TH3F>& procshape_3D, TString newname=""){
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

  RooCmdArg ycmd;
  RooCmdArg zcmd;
  if (deps.size()>1) ycmd = YVar(*(deps.at(1)));
  if (deps.size()>2) zcmd = ZVar(*(deps.at(2)));

  TString procname=proc.name;
  if (newname!="") procname=newname;
  TString tplname = "T_";
  tplname += procname;

  if (ndims==3){
    TH3F* tpl=(TH3F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), ycmd, zcmd);
    double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
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
    TH2F* tpl=(TH2F*) proc.pdf->createHistogram(tplname+"_Copy", *(deps.at(0)), ycmd, zcmd);
    double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
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
  cout << "Work on template " << tplname << " is complete!" << endl;
  return ndims;
}


void getDistributions(TString cinputdir, int onORoffshell=0, bool isEnriched=true){
  gStyle->SetOptStat(0);

  const unsigned int nsqrts=2;
  TString sqrtsnames[nsqrts]={ "13TeV_2016", "13TeV_2017" };
  const unsigned int nchans=3;
  TString channames[nchans]={ "4mu", "4e", "2e2mu" };
  const unsigned int ncats=3;
  TString catnames[ncats]={ "Untagged", "JJVBFTagged", "HadVHTagged" };

  TDirectory* curdir=gDirectory;

  TString failabel="f_{ai}";
  if (cinputdir.Contains("a3")) failabel="f_{a3}";
  else if (cinputdir.Contains("a2")) failabel="f_{a2}";
  else if (cinputdir.Contains("L1")) failabel="f_{#Lambda1}";

  // Determine alternative model parameters
  constexpr float GHref=4.07;
  unordered_map<TString, float> val_fai1ALT; val_fai1ALT["RV"]=1; val_fai1ALT["RF"]=1; val_fai1ALT["fai1"]=0;
  unordered_map<TString, float> val_GGsmALT; val_GGsmALT["RV"]=1; val_GGsmALT["RF"]=1; val_GGsmALT["GGsm"]=1;
  if (onORoffshell==1){ // Offshell
    if (cinputdir.Contains("a3")){ val_fai1ALT["RV"]=0.22; val_fai1ALT["RF"]=1.06; val_fai1ALT["fai1"]=0.01; }
    else if (cinputdir.Contains("a2")){ val_fai1ALT["RV"]=0.3; val_fai1ALT["RF"]=1.3; val_fai1ALT["fai1"]=-0.015; }
    else if (cinputdir.Contains("L1")){ val_fai1ALT["RV"]=0.1; val_fai1ALT["RF"]=0.85; val_fai1ALT["fai1"]=0.015; }
    val_GGsmALT["RV"]=0.29; val_GGsmALT["RF"]=0.88; val_GGsmALT["GGsm"]=14./GHref;
  }
  else{ // Onshell
    if (cinputdir.Contains("a3")){ val_fai1ALT["RV"]=0.22; val_fai1ALT["RF"]=1.06; val_fai1ALT["fai1"]=0.01; }
    else if (cinputdir.Contains("a2")){ val_fai1ALT["RV"]=0.3; val_fai1ALT["RF"]=1.3; val_fai1ALT["fai1"]=-0.015; }
    else if (cinputdir.Contains("L1")){ val_fai1ALT["RV"]=0.1; val_fai1ALT["RF"]=0.85; val_fai1ALT["fai1"]=0.015; }
    else{ val_GGsmALT["RV"]=0.29; val_GGsmALT["RF"]=0.88; val_GGsmALT["GGsm"]=1; }
  }

  for (unsigned int icat=0; icat<ncats; icat++){
    TString const& catname = catnames[icat];
    unordered_map<TString, TH2F> procshape_2D;
    unordered_map<TString, TH3F> procshape_3D;
    cout << "Acquiring shapes for category " << catname << endl;
    unsigned int ndims=0;

    // Get process order/labels
    std::vector<TString> proc_order, proc_label;
    std::vector<int> proc_color;
    std::vector<int> proc_code;
    if (onORoffshell==1){ // Offshell
      if (catname=="Untagged") proc_order=std::vector<TString>{ "qqZZ", "Zjets", "VVZZ_offshell", "ggZZ_offshell", "total_GGsmALT" };
      else if (catname=="JJVBFTagged") proc_order=std::vector<TString>{ "qqZZ", "Zjets", "ggZZ_offshell", "VVZZ_offshell", "total_GGsmALT" };
      else if (catname=="HadVHTagged") proc_order=std::vector<TString>{ "qqZZ", "Zjets", "ggZZ_offshell", "VVZZ_offshell", "total_GGsmALT" };
      if (!cinputdir.Contains("SM")) proc_order.push_back("total_fai1ALT");
    }
    else{
      if (catname=="Untagged") proc_order=std::vector<TString>{ "bkg_zz", "zjets", "VVZZ", "ggH", "total_fai1ALT" };
      else if (catname=="JJVBFTagged") proc_order=std::vector<TString>{ "bkg_zz", "zjets", "ggH", "VVZZ", "total_fai1ALT" };
      else if (catname=="HadVHTagged") proc_order=std::vector<TString>{ "bkg_zz", "zjets", "ggH", "VVZZ", "total_fai1ALT" };
    }
    for (auto const& p:proc_order){
      if (p=="bkg_qqzz" || p=="qqZZ"){ proc_color.push_back(int(kAzure-2)); proc_label.push_back("q#bar{q}#rightarrow4l bkg."); proc_code.push_back(0); }
      else if (p=="bkg_gg"){ proc_color.push_back(int(kBlue)); proc_label.push_back("gg#rightarrow4l bkg."); proc_code.push_back(0); }
      else if (p=="zjets" || p=="Zjets"){ proc_color.push_back(int(kGreen+2)); proc_label.push_back("Z+jets"); proc_code.push_back(0); }
      else if (p=="bkg_vv"){ proc_color.push_back(int(kPink+9)); proc_label.push_back("EW bkg."); proc_code.push_back(0); }
      else if (p=="bkg_zz"){ proc_color.push_back(int(kAzure-2)); proc_label.push_back("ZZ/Z#gamma*/#gamma*#gamma*#rightarrow4l bkg."); proc_code.push_back(0); }
      else if (p=="ggZZ_offshell" || p=="ggH"){
        proc_color.push_back(int(kOrange-3));
        if (p=="ggH"){ proc_label.push_back("gg#rightarrow4l sig."); proc_code.push_back(1); }
        else if (p=="ggZZ_offshell"){ proc_label.push_back("gg#rightarrow4l SM total"); proc_code.push_back(2); }
      }
      else if (p=="VVZZ" || p=="VVZZ_offshell"){
        proc_color.push_back(int(kViolet));
        if (p=="VVZZ"){ proc_label.push_back("EW sig."); proc_code.push_back(1); }
        else if (p=="VVZZ_offshell"){ proc_label.push_back("EW SM total"); proc_code.push_back(2); }
      }
      else if (p.Contains("ALT")){
        if (onORoffshell==1){ // Offshell
          if (p.Contains("GGsmALT")){ proc_label.push_back(Form("Total (%s=0, #Gamma_{H}=%.0f MeV)", failabel.Data(), val_GGsmALT["GGsm"]*GHref)); proc_color.push_back(int(kCyan+2)); proc_code.push_back(-2); }
          else if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total (%s=%.3f, #Gamma_{H}=#Gamma_{H}^{SM})", failabel.Data(), val_fai1ALT["fai1"])); proc_color.push_back(int(kRed)); proc_code.push_back(-2); }
        }
        else{ // Onshell
          if (p.Contains("fai1ALT")){ proc_label.push_back(Form("Total signal, %s=%.3f", failabel.Data(), val_fai1ALT["fai1"])); proc_color.push_back(int(kRed)); proc_code.push_back(-1); }
        }
      }
      else cerr << p << " is unmatched for labeling!" << endl;
    }

    for (unsigned int ic=0; ic<nchans; ic++){
      TString const& channame = channames[ic];
      cout << "\t- Channel:  " << channame << endl;
      for (unsigned int is=0; is<nsqrts; is++){
        TString const& sqrtsname = sqrtsnames[is];
        cout << "\t- Data period:  " << sqrtsname << endl;
        TString cinput_main = cinputdir + "/" + (onORoffshell ? "Offshell_" : "Onshell_") + sqrtsname + "/hzz" + channame + "_" + catname;
        TString cinput_file = cinput_main + ".input.root";
        TString cinput_dc = cinput_main + ".txt";

        vector<TString> procname;
        vector<double> procrate;

        // Get process names
        ifstream tin(cinput_dc.Data());
        char line[512];
        while (string(line).find("process")==string::npos) tin.getline(line, 512);
        char* chars_array = strtok(line, " ");
        chars_array = strtok(NULL, " ");
        while (chars_array){
          procname.push_back(TString(chars_array));
          chars_array = strtok(NULL, " ");
        }
        // Get process rates
        while (string(line).find("rate")==string::npos) tin.getline(line, 512);
        chars_array = strtok(line, " ");
        chars_array = strtok(NULL, " ");
        while (chars_array){
          procrate.push_back(atof(chars_array));
          chars_array = strtok(NULL, " ");
        }
        tin.close();
        assert(procrate.size()==procname.size());

        // Construct the processes
        unordered_map<TString, process_spec> procSpecs;
        {
          std::vector<double>::iterator it_rate=procrate.begin();
          for (TString const& pname:procname){
            cout << "Process " << pname << " is found in the datacard" << endl;
            procSpecs[pname]=process_spec(pname, *it_rate);
            it_rate++;
          }
        }

        // Open the input
        TFile* finput = TFile::Open(cinput_file, "read");
        RooWorkspace* ws = (RooWorkspace*) finput->Get("w");
        unordered_map<TString, RooRealVar*> controlVars;
        TString varsToCheck[10]={
          "R", "RF", "RF_13TeV", "RV", "RV_13TeV", "R_13TeV", "CMS_zz4l_fai1", "GGsm", "kbkg_gg", "kbkg_VBF"
        };
        for (unsigned int v=0; v<10; v++){
          controlVars[varsToCheck[v]]=ws->var(varsToCheck[v]);
          if (controlVars[varsToCheck[v]]){
            if (varsToCheck[v]!="CMS_zz4l_fai1") controlVars[varsToCheck[v]]->setVal(1);
            else controlVars[varsToCheck[v]]->setVal(0);
          }
          else cerr << varsToCheck[v] << " could not be found!" << endl;
        }
        RooDataSet* data = (RooDataSet*) ws->data("data_obs");

        // Get process pdfs
        for (auto const& pname:procname){
          curdir->cd();
          cout << "Extracting the pdf and norm from process " << pname << endl;
          RooAbsPdf* pdf = ws->pdf(pname);
          RooAbsReal* norm = (RooAbsReal*) ws->factory(pname+"_norm");
          if (!norm) cout << "Warning: " << pname << "_norm is not found." << endl;
          if (!pdf) cerr << pname << " pdf could not be found." << endl;
          procSpecs[pname].pdf=pdf;
          procSpecs[pname].norm=norm;

          if (!cinputdir.Contains("SM") && !onORoffshell){ // AC onshell
            if (pname.Contains("ggZZ")){
              controlVars["kbkg_gg"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "ggH");
              controlVars["kbkg_gg"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_gg");
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              controlVars["kbkg_gg"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
              controlVars["kbkg_gg"]->setVal(1);
            }
            else if (pname.Contains("VBF")){
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "VVZZ");
              controlVars["kbkg_VBF"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_vv");
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
              controlVars["kbkg_VBF"]->setVal(1);
            }
            else if (pname.Contains("ZH") || pname.Contains("WH")){
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "VVZZ");
              controlVars["kbkg_VBF"]->setVal(1);
              controlVars["R"]->setVal(0);
              //ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_vv");
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_zz");
              controlVars["R"]->setVal(1);

              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              controlVars["kbkg_VBF"]->setVal(0);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
              controlVars["kbkg_VBF"]->setVal(1);
            }
            else if (pname.Contains("bkg_qqzz")){
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_zz");
            }
          }
          else if (!onORoffshell && (pname.Contains("ZH") || pname.Contains("WH") || pname.Contains("VBF"))){ // On-shell SM
            ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "VVZZ");
          }
          //else if (!onORoffshell && (pname.Contains("bkg_zzz") || pname.Contains("bkg_wzz") || pname.Contains("bkg_vbs"))){ // On-shell SM
          else if (!onORoffshell && (pname.Contains("bkg_zzz") || pname.Contains("bkg_wzz") || pname.Contains("bkg_vbs") || pname.Contains("bkg_gg") || pname.Contains("bkg_qqzz"))){ // On-shell SM
            //ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_vv");
            ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "bkg_zz");
          }
          else if (onORoffshell && (pname.Contains("ggZZ_offshell") || pname.Contains("VVZZ_offshell"))){
            controlVars["GGsm"]->setVal(val_GGsmALT["GGsm"]);
            controlVars["RV"]->setVal(val_GGsmALT["RV"]);
            controlVars["RF"]->setVal(val_GGsmALT["RF"]);
            ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "total_GGsmALT");
            controlVars["GGsm"]->setVal(1);
            controlVars["RV"]->setVal(1);
            controlVars["RF"]->setVal(1);

            if (!cinputdir.Contains("SM")){
              controlVars["CMS_zz4l_fai1"]->setVal(val_fai1ALT["fai1"]);
              controlVars["RV"]->setVal(val_fai1ALT["RV"]);
              controlVars["RF"]->setVal(val_fai1ALT["RF"]);
              ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D, "total_fai1ALT");
              controlVars["CMS_zz4l_fai1"]->setVal(0);
              controlVars["RV"]->setVal(1);
              controlVars["RF"]->setVal(1);
            }

          }

          ndims = extractTemplates(procSpecs[pname], data, procshape_2D, procshape_3D);
        }
        finput->Close();
        curdir->cd();
      }
    }

    if (ndims==3){
      vector<TString> varlabels;
      // KD1
      if (onORoffshell==1 || cinputdir.Contains("SM")) varlabels.push_back("m_{4l} (GeV)");
      else{ // On-shell AC
        if (catname=="Untagged") varlabels.push_back("D_{bkg}");
        else if (catname=="JJVBFTagged") varlabels.push_back("D_{bkg,m4l}^{VBF+dec}");
        else if (catname=="HadVHTagged") varlabels.push_back("D_{bkg,m4l}^{VH+dec}");
      }
      // KD2
      if (onORoffshell==1 || cinputdir.Contains("SM")){
        if (catname=="Untagged") varlabels.push_back("D_{bkg}^{kin}");
        else if (catname=="JJVBFTagged") varlabels.push_back("D_{bkg}^{VBF+dec}");
        else if (catname=="HadVHTagged") varlabels.push_back("D_{bkg}^{VH+dec}");
      }
      else{ // On-shell AC
        if (catname=="Untagged"){
          if (cinputdir.Contains("a3")) varlabels.push_back("D_{0-}^{dec}");
          else if (cinputdir.Contains("a2")) varlabels.push_back("D_{0h+}^{dec}");
          else if (cinputdir.Contains("L1")) varlabels.push_back("D_{#Lambda1}^{dec}");
        }
        else if (catname=="JJVBFTagged"){
          if (cinputdir.Contains("a3")) varlabels.push_back("D_{0-}^{VBF+dec}");
          else if (cinputdir.Contains("a2")) varlabels.push_back("D_{0h+}^{VBF+dec}");
          else if (cinputdir.Contains("L1")) varlabels.push_back("D_{#Lambda1}^{VBF+dec}");
        }
        else if (catname=="HadVHTagged"){
          if (cinputdir.Contains("a3")) varlabels.push_back("D_{0-}^{VH+dec}");
          else if (cinputdir.Contains("a2")) varlabels.push_back("D_{0h+}^{VH+dec}");
          else if (cinputdir.Contains("L1")) varlabels.push_back("D_{#Lambda1}^{VH+dec}");
        }
      }
      // KD3
      if (onORoffshell==1){
        if (cinputdir.Contains("SM")){
          if (catname=="Untagged") varlabels.push_back("D_{bsi}^{gg,dec}");
          else if (catname=="JJVBFTagged") varlabels.push_back("D_{bsi}^{VBF+dec}");
          else if (catname=="HadVHTagged") varlabels.push_back("D_{bsi}^{VH+dec}");
        }
        else{
          if (catname=="Untagged"){
            if (cinputdir.Contains("a3")) varlabels.push_back("D_{0-}^{dec}");
            else if (cinputdir.Contains("a2")) varlabels.push_back("D_{0h+}^{dec}");
            else if (cinputdir.Contains("L1")) varlabels.push_back("D_{#Lambda1}^{dec}");
          }
          else if (catname=="JJVBFTagged"){
            if (cinputdir.Contains("a3")) varlabels.push_back("D_{0-}^{VBF+dec}");
            else if (cinputdir.Contains("a2")) varlabels.push_back("D_{0h+}^{VBF+dec}");
            else if (cinputdir.Contains("L1")) varlabels.push_back("D_{#Lambda1}^{VBF+dec}");
          }
          else if (catname=="HadVHTagged"){
            if (cinputdir.Contains("a3")) varlabels.push_back("D_{0-}^{VH+dec}");
            else if (cinputdir.Contains("a2")) varlabels.push_back("D_{0h+}^{VH+dec}");
            else if (cinputdir.Contains("L1")) varlabels.push_back("D_{#Lambda1}^{VH+dec}");
          }
        }
      }
      else{ // On-shell AC
        if (catname=="Untagged"){
          if (cinputdir.Contains("a3")) varlabels.push_back("D_{CP}^{dec}");
          else if (cinputdir.Contains("a2")) varlabels.push_back("D_{int}^{dec}");
          else if (cinputdir.Contains("L1")) varlabels.push_back("D_{0h+}^{dec}");
        }
        else if (catname=="JJVBFTagged"){
          if (cinputdir.Contains("a3")) varlabels.push_back("D_{CP}^{VBF}");
          else if (cinputdir.Contains("a2")) varlabels.push_back("D_{int}^{VBF}");
          else if (cinputdir.Contains("L1")) varlabels.push_back("D_{0h+}^{VBF+dec}");
        }
        else if (catname=="HadVHTagged"){
          if (cinputdir.Contains("a3")) varlabels.push_back("D_{CP}^{VH}");
          else if (cinputdir.Contains("a2")) varlabels.push_back("D_{int}^{VH}");
          else if (cinputdir.Contains("L1")) varlabels.push_back("D_{0h+}^{VH+dec}");
        }
      }

      cout << "\t- Variable labels: ";
      for (auto& vl:varlabels) cout << vl << " ";
      cout << endl;

      curdir->cd();
      int massDim=-1;
      int kdDim=-1;
      for (unsigned int idim=0; idim<ndims; idim++){
        if (varlabels.at(idim).Contains("D_{bkg}")) kdDim=idim;
        if (varlabels.at(idim).Contains("m_{4l}")) massDim=idim;
      }

      for (unsigned int idim=0; idim<ndims; idim++){
        TString dimname = Form("_KD%i", idim+1);
        unordered_map<TString, TH1F*> procdist;
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
            const float valMassCut=300;
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
            }

            const float valKDCut=0.5;
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
          procdist[it->first]=htmp;
          cout << "\t\t- Final histogram integral: " << getHistogramIntegralAndError(procdist[it->first], 1, procdist[it->first]->GetNbinsX(), false, nullptr) << endl;
        }

        for (unsigned int ip=1; ip<proc_order.size(); ip++){
          if (!proc_order.at(ip).Contains("ALT")) procdist[proc_order.at(ip)]->Add(procdist[proc_order.at(ip-1)]);
          else{
            for (unsigned int jp=ip-1; jp>0; jp--){
              if (proc_code[jp]==0){
                procdist[proc_order.at(ip)]->Add(procdist[proc_order.at(jp)]); // Add the first bkg process and break
                break;
              }
            }
          }
          cout << proc_order.at(ip) << " integral: " << getHistogramIntegralAndError(procdist[proc_order.at(ip)], 1, procdist[proc_order.at(ip)]->GetNbinsX(), false, nullptr) << endl;
        }

        // Divide by bin width when plotting off-shell mass
        if (onORoffshell==1 && (int) idim==massDim){ for (auto proc:proc_order) divideBinWidth(procdist[proc]); }

        // Draw
        TCanvas canvas(
          TString((isEnriched ? "c_SignalEnriched_" : "c_")) + (onORoffshell ? "Offshell_" : "Onshell_") + catname + "_" + (cinputdir.Contains("a3") ? "a3" : (cinputdir.Contains("a2") ? "a2" : (cinputdir.Contains("L1") ? "L1" : "SM"))) + dimname,
          "", 8, 30, 800, 800
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

        TLegend legend(0.55, 0.90-0.10/3.*2.*float(proc_order.size()), 0.80, 0.90);
        legend.SetBorderSize(0);
        legend.SetTextFont(42);
        legend.SetTextSize(0.03);
        legend.SetLineColor(1);
        legend.SetLineStyle(1);
        legend.SetLineWidth(1);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);

        TText* text;
        TPaveText pt(0.15, 0.93, 0.85, 1, "brNDC");
        pt.SetBorderSize(0);
        pt.SetFillStyle(0);
        pt.SetTextAlign(12);
        pt.SetTextFont(42);
        pt.SetTextSize(0.045);
        text = pt.AddText(0.025, 0.45, "#font[61]{CMS}");
        text->SetTextSize(0.044);
        text = pt.AddText(0.165, 0.42, "#font[52]{Work in progress}");
        text->SetTextSize(0.0315);
        int theSqrts=13;
        TString cErgTev = Form("#font[42]{77.3 fb^{-1} (%i TeV)}", theSqrts);
        text = pt.AddText(0.82, 0.45, cErgTev);
        text->SetTextSize(0.0315);

        float ymax=-1;
        float xmin=-1, xmax=-1;
        for (unsigned int ip=0; ip<proc_order.size(); ip++){
          if (onORoffshell && idim==0){ xmin=220; xmax=1000; }
          else{
            xmin=procdist[proc_order.at(ip)]->GetXaxis()->GetBinLowEdge(1);
            xmax=procdist[proc_order.at(ip)]->GetXaxis()->GetBinUpEdge(procdist[proc_order.at(ip)]->GetNbinsX());
          }
          cout << "Process " << proc_order.at(ip) << " color: " << proc_color[ip] << endl;
          if (proc_order.at(ip).Contains("ALT")){
            if (proc_order.at(ip).Contains("GGsm")) procdist[proc_order.at(ip)]->SetLineStyle(2);
            else if (proc_order.at(ip).Contains("fai1")) procdist[proc_order.at(ip)]->SetLineStyle(7);
            procdist[proc_order.at(ip)]->SetMarkerColor(proc_color[ip]);
            procdist[proc_order.at(ip)]->SetLineColor(proc_color[ip]);
          }
          else if (proc_code.at(ip)==0){
            procdist[proc_order.at(ip)]->SetMarkerColor(proc_color[ip]);
            procdist[proc_order.at(ip)]->SetLineColorAlpha(kBlack, 0.5);
            procdist[proc_order.at(ip)]->SetFillColor(proc_color[ip]);
            procdist[proc_order.at(ip)]->SetFillStyle(1001);
          }
          else{
            procdist[proc_order.at(ip)]->SetMarkerColor(proc_color[ip]);
            procdist[proc_order.at(ip)]->SetLineColorAlpha(proc_color[ip], 0.9);
            procdist[proc_order.at(ip)]->SetFillColor(0);
            procdist[proc_order.at(ip)]->SetFillStyle(1001);
          }
          procdist[proc_order.at(ip)]->SetLineWidth(2);

          int binXlow = procdist[proc_order.at(ip)]->GetXaxis()->FindBin(xmin);
          int binXhigh = procdist[proc_order.at(ip)]->GetXaxis()->FindBin(xmax);

          procdist[proc_order.at(ip)]->GetXaxis()->SetRangeUser(xmin, xmax);
          procdist[proc_order.at(ip)]->GetXaxis()->SetNdivisions(505);
          procdist[proc_order.at(ip)]->GetXaxis()->SetLabelFont(42);
          procdist[proc_order.at(ip)]->GetXaxis()->SetLabelOffset(0.007);
          procdist[proc_order.at(ip)]->GetXaxis()->SetLabelSize(0.04);
          procdist[proc_order.at(ip)]->GetXaxis()->SetTitleSize(0.06);
          procdist[proc_order.at(ip)]->GetXaxis()->SetTitleOffset(0.9);
          procdist[proc_order.at(ip)]->GetXaxis()->SetTitleFont(42);
          procdist[proc_order.at(ip)]->GetYaxis()->SetNdivisions(505);
          procdist[proc_order.at(ip)]->GetYaxis()->SetLabelFont(42);
          procdist[proc_order.at(ip)]->GetYaxis()->SetLabelOffset(0.007);
          procdist[proc_order.at(ip)]->GetYaxis()->SetLabelSize(0.04);
          procdist[proc_order.at(ip)]->GetYaxis()->SetTitleSize(0.06);
          procdist[proc_order.at(ip)]->GetYaxis()->SetTitleOffset(1.1);
          procdist[proc_order.at(ip)]->GetYaxis()->SetTitleFont(42);
          procdist[proc_order.at(ip)]->GetYaxis()->SetTitle("Events / bin");
          for (int ix=binXlow; ix<binXhigh; ix++){
            float bc = procdist[proc_order.at(ip)]->GetBinContent(ix);
            if (bc!=0.) ymax = std::max(bc, ymax);
          }
          legend.AddEntry(procdist[proc_order.at(ip)], proc_label.at(ip), "f");
        }
        for (unsigned int ip=proc_order.size(); ip>0; ip--){
          cout << "\t- Drawing " << proc_order.at(ip-1) << endl;
          procdist[proc_order.at(ip-1)]->GetYaxis()->SetRangeUser(0, ymax*1.4);
          procdist[proc_order.at(ip-1)]->Draw((ip==proc_order.size() ? "hist" : "histsame"));
        }
        // Re-draw ALT
        for (unsigned int ip=proc_order.size(); ip>0; ip--){
          if (!proc_order.at(ip-1).Contains("ALT")) continue;
          cout << "\t- Drawing " << proc_order.at(ip-1) << endl;
          procdist[proc_order.at(ip-1)]->GetYaxis()->SetRangeUser(0, ymax*1.4);
          procdist[proc_order.at(ip-1)]->Draw("histsame");
        }
        legend.Draw();
        pt.Draw();
        canvas.RedrawAxis();
        canvas.Modified();
        canvas.Update();
        canvas.SaveAs(TString(canvas.GetName())+".pdf");
        canvas.Close();
        curdir->cd();

        for (auto it=procdist.begin(); it!=procdist.end(); it++){
          TH1F* htmp = it->second;
          cout << "\t- Deleting histogram " << htmp->GetName() << endl;
          delete htmp;
        }
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
