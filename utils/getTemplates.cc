#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include "TIterator.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TH3F.h"
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
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"

using namespace RooFit;
using namespace std;

struct process_spec{
  RooAbsPdf* pdf;
  unordered_map<string, RooAbsPdf*> pdf_shape;

  RooAbsReal* norm;
  double rate;
  TString name;
  vector<TH3F*> templates;
  unordered_map<string, pair<string, string>> systematics;

  process_spec() : pdf(0), norm(0), rate(0){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()){}
  process_spec(RooAbsPdf* pdf_, RooAbsReal* norm_, double rate_, vector<TH3F*> templates_) : pdf(pdf_), norm(norm_), rate(rate_), name(pdf->GetName()), templates(templates_){}
  process_spec(const process_spec& other) : pdf(other.pdf), norm(other.norm), rate(other.rate), templates(other.templates){}

  void setSystematic(string systname, string systtype, string systline){ systematics[systname] = pair<string, string>(systtype, systline); }
  void setShapePdf(string systname, RooAbsPdf* pdf_){ if (pdf_!=0)pdf_shape[systname] = pdf_; else cout << name << "_" << systname << " pdf does not exist!" << endl; }
  void writeTemplates(TFile* fin){ for (unsigned int it=0; it<templates.size(); it++) fin->WriteTObject(templates.at(it)); }
  void deleteTemplates(){ for (unsigned int it=0; it<templates.size(); it++) delete templates.at(it); }
  //~process_spec(){ for (unsigned int it=0; it<templates.size(); it++) delete templates.at(it); }

};

void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);
void extractTemplates(process_spec& proc, RooDataSet* data, string shapename, bool scale_width);
void extractTemplates_ggLike_fai1(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);
void extractTemplates_VVLike_fai1(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl);

void getDataTree(TString cinput){
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
  string catnames[ncats]={ "Inclusive", "VBFtagged", "VHHadrtagged", "Untagged" };
  string catname;
  for (unsigned int ic=0; ic<ncats; ic++){
    if (strinput.find(catnames[ic])!=string::npos) catname=catnames[ic];
  }

  TFile* finput = TFile::Open(cinput+".input.root", "read");
  RooWorkspace* ws = (RooWorkspace*)finput->Get("w");
  RooDataSet* data = (RooDataSet*)ws->data("data_obs");
  int nevents=data->sumEntries();
  int nvars=((const RooArgSet*)data->get())->getSize();
  double* KD = new double[nvars];
  double mass;

  TString coutput_root;
  TFile* foutput;
  coutput_root = Form("test/13TeV/CMSdata/hzz%s_%s_13TeV.root", channame.c_str(), catname.c_str());
  foutput = TFile::Open(coutput_root, "recreate");

  TTree* t = new TTree("data_obs", "");
  for (int iv=0; iv<nvars; iv++) t->Branch(Form("KD%i", iv+1), KD+iv);

  for (int ev=0; ev<nevents; ev++){
    const RooArgSet* args = (const RooArgSet*)data->get(ev);
    TIterator* coefIter = args->createIterator();
    const RooAbsArg* coef;
    unsigned int ik=0;
    while ((coef = (const RooAbsArg*)coefIter->Next())){
      const RooAbsReal* rvar = dynamic_cast<const RooAbsReal*>(coef);
      TString rname = rvar->GetName();
      if (rname.Contains("mass") || rname.Contains("Mass")){
        mass=rvar->getVal();
        if (t->GetBranchStatus("mass")==0) t->Branch("mass", &mass);
      }
      KD[ik]=rvar->getVal();
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
void getTemplates(TString cinput, double lumiScale=1, bool scale_width=true){
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
  string catnames[ncats]={ "Inclusive", "VBFtagged", "VHHadrtagged", "Untagged" };
  string catname;
  for (unsigned int ic=0; ic<ncats; ic++){
    if (strinput.find(catnames[ic])!=string::npos) catname=catnames[ic];
  }

  //TString coutput = splitinput.at(splitinput.size()-1);
  //ofstream tout((coutput+".input.txt").Data());
  TString coutput_txt = Form("test/inputs_%s_%s%s", channame.c_str(), catname.c_str(), ".txt");
  ofstream tout(coutput_txt.Data());

  TFile* finput = TFile::Open(cinput+".input.root", "read");
  RooWorkspace* ws = (RooWorkspace*)finput->Get("w");
  TString varsToCheck[6]={
    "R", "RF", "RF_13TeV", "RV", "RV_13TeV", "R_13TeV"
  };
  for (unsigned int v=0; v<6; v++){
    if (ws->var(varsToCheck[v])!=0) ((RooRealVar*)ws->var(varsToCheck[v]))->setVal(1);
  }
  ws->Print("v");

  ifstream tin((cinput+".txt").Data());
  char line[512];

  RooDataSet* data = (RooDataSet*)ws->data("data_obs");

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
  while (string(line).find("rate")==string::npos) tin.getline(line, 512);
  chars_array = strtok(line, " ");
  chars_array = strtok(NULL, " ");
  vector<double> procrate;
  while (chars_array){
    procrate.push_back(atof(chars_array));
    chars_array = strtok(NULL, " ");
  }

  // Get process pdfs
  unordered_map<const char*, process_spec> procSpecs;
  for (unsigned int ip=0; ip<procname.size(); ip++){
    cout << procname.at(ip) << ": " << procrate.at(ip) << endl;

    RooAbsPdf* pdf = ws->pdf(procname.at(ip));
    RooAbsReal* norm = 0; norm = (RooAbsReal*)ws->factory(procname.at(ip)+"_norm");
    if (norm==0) cout << "Warning: " << procname.at(ip) << "_norm is not found." << endl;

    if (pdf==0) cerr << procname.at(ip) << " pdf could not be found." << endl;
    else procSpecs[procname.at(ip).Data()]=process_spec(pdf, norm, procrate.at(ip));
  }

  // Get systemtics
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
      if (isShape || isLog){
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
            accumulate += string(procname.at(ip).Data()) + ":" + systline + " ";
          }
        }
        if (isLog) logSyst[systname] = accumulate;
      }
      else{
        RooAbsReal* systvar = (RooAbsReal*) ws->var(systname.c_str());
        double defaultVal = systvar->getVal();
        if (systvar->hasClients()){
          TIterator* clientsIter = systvar->clientIterator();
          RooAbsArg* client;
          while ((client = (RooAbsArg*)clientsIter->Next())){
            for (unsigned int ip=0; ip<procname.size(); ip++){
              RooAbsPdf* pdf = procSpecs[procname.at(ip).Data()].pdf;
              if (client->GetName() == pdf->GetName()){
                string systline="";
                for (unsigned int ip=2; ip<systdist.size(); ip++){
                  if (ip>2)systline += ":";
                  systline += systdist.at(ip);
                }
                procSpecs[procname.at(ip).Data()].setSystematic(systname, systtype, systline);
                accumulate += string(procname.at(ip).Data()) + ":" + systline + " ";
              }
            }
          }
        }
        paramSyst[systname] = accumulate;
      }

    }
  }

  // Write input file
  tout << "sqrts " << 13 << endl;
  tout << "decay " << channame << endl;
  tout << "lumi " << lumiScale << endl;
  tout << "category " << catname << endl;

  // Write channels
  for (unsigned int ip=0; ip<procname.size(); ip++) tout << "channel " << procname.at(ip) << " 1 -1 " << (procSpecs[procname.at(ip).Data()].name.Contains("bkg") ? 0 : 1) << endl;

  // Write systemtics
  for (unsigned int ip=0; ip<procname.size(); ip++){
    for (auto syst = procSpecs[procname.at(ip).Data()].systematics.begin(); syst != procSpecs[procname.at(ip).Data()].systematics.end(); ++syst){
      if (syst->second.first=="shape1"){
        tout << "systematic " << syst->first << " template ";
        tout << procname.at(ip) << ":0:1" << endl;
      }
    }
  }
  for (auto syst = logSyst.begin(); syst != logSyst.end(); ++syst) tout << "systematic " << syst->first << " lnN " << syst->second << endl;
  for (auto syst = paramSyst.begin(); syst != paramSyst.end(); ++syst) tout << "systematic " << syst->first << " template " << syst->second << endl;

  for (unsigned int ip=0; ip<procname.size(); ip++){
    TFile* foutput;
    TString coutput_root;
    
    coutput_root = Form("test/13TeV/HtoZZ%s_%s_FinalTemplates_%s_%s%s", channame.c_str(), catname.c_str(), procname.at(ip).Data(), "Nominal", ".root");
    foutput = TFile::Open(coutput_root, "recreate");
    extractTemplates(procSpecs[procname.at(ip).Data()], data, "", scale_width);
    for (unsigned int it=0; it<procSpecs[procname.at(ip).Data()].templates.size(); it++) procSpecs[procname.at(ip).Data()].templates.at(it)->Scale(1./lumiScale);
    procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
    procSpecs[procname.at(ip).Data()].deleteTemplates();
    foutput->Close();

    for (auto syst = procSpecs[procname.at(ip).Data()].systematics.begin(); syst != procSpecs[procname.at(ip).Data()].systematics.end(); ++syst){
      if (syst->second.first=="param"){
        RooRealVar* systvar = (RooRealVar*)ws->function(syst->first.c_str());
        if (systvar==0){
          cout << syst->first << " could not be found." << endl;
          continue;
        }
        for (unsigned int is=0; is<2; is++){
          systvar->setVal(double(2*is-1));
          coutput_root = Form("test/13TeV/HtoZZ%s_%s_FinalTemplates_%s_%s%s%s", channame.c_str(), catname.c_str(), procname.at(ip).Data(), syst->first.c_str(), (is==0 ? "Down" : "Up"), ".root");
          foutput = TFile::Open(coutput_root, "recreate");
          extractTemplates(procSpecs[procname.at(ip).Data()], data, "", scale_width);
          for (unsigned int it=0; it<procSpecs[procname.at(ip).Data()].templates.size(); it++) procSpecs[procname.at(ip).Data()].templates.at(it)->Scale(1./lumiScale);
          procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
          procSpecs[procname.at(ip).Data()].deleteTemplates();
          foutput->Close();
        }
      }
      else if (syst->second.first=="shape1"){
        for (unsigned int is=0; is<2; is++){
          string syst_du = syst->first + (is==0 ? "Down" : "Up");
          coutput_root = Form("test/13TeV/HtoZZ%s_%s_FinalTemplates_%s_%s%s", channame.c_str(), catname.c_str(), procname.at(ip).Data(), syst_du.c_str(), ".root");
          foutput = TFile::Open(coutput_root, "recreate");
          extractTemplates(procSpecs[procname.at(ip).Data()], data, syst_du, scale_width);
          for (unsigned int it=0; it<procSpecs[procname.at(ip).Data()].templates.size(); it++) procSpecs[procname.at(ip).Data()].templates.at(it)->Scale(1./lumiScale);
          procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
          procSpecs[procname.at(ip).Data()].deleteTemplates();
          foutput->Close();
        }
      }


    }
  }

  tin.close();
  finput->Close();
  tout.close();
}

void extractTemplates(process_spec& proc, RooDataSet* data, string shapename, bool scale_width){
  vector<TH3F*> templates;

  vector<RooRealVar*> deps;
  RooArgSet* depList = proc.pdf->getDependents(data);
  depList->Print("v");
  TIterator* coefIter = depList->createIterator();
  RooAbsArg* coef;
  while ((coef = (RooAbsArg*)coefIter->Next())) deps.push_back((RooRealVar*)coef);
  delete coefIter;
  delete depList;

  vector<RooRealVar*> pars;
  RooRealVar* fai1=0;
  RooRealVar* GGsm=0;
  RooArgSet* parList = proc.pdf->getParameters(data);
  parList->Print("v");
  coefIter = parList->createIterator();
  coef=0;
  while ((coef = (RooAbsArg*)coefIter->Next())){
    pars.push_back((RooRealVar*)coef);
    if (TString(pars.at(pars.size()-1)->GetName()).Contains("fai1")) fai1 = pars.at(pars.size()-1);
    else if (TString(pars.at(pars.size()-1)->GetName()).Contains("GGsm")) GGsm = pars.at(pars.size()-1);
  }
  delete coefIter;
  delete parList;

  RooCmdArg ycmd;
  RooCmdArg zcmd;
  if (deps.size()>1) ycmd = YVar(*(deps.at(1)));
  if (deps.size()>2) zcmd = ZVar(*(deps.at(2)));

  TString tplname = "T_";
  tplname += proc.name;

  if (proc.name.Contains("bkg")){
    TH3F* tpl=(TH3F*)proc.pdf->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd);
    double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
    double integral = 1;
    if (scale_width) integral = tpl->Integral("width");
    else integral = tpl->Integral();
    double scale = normval/integral;
    cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
    if (shapename!=""){
      delete tpl;
      tpl=(TH3F*)proc.pdf_shape[shapename]->createHistogram(tplname, *(deps.at(0)), ycmd, zcmd);
    }
    tpl->SetName(tplname);
    tpl->SetTitle("");
    tpl->Scale(scale);
    templates.push_back(tpl);
  }
  else if (proc.name.Contains("ggH") || proc.name.Contains("ttH")){
    if (fai1!=0 && GGsm==0){
      vector<TH3F*> intpl;
      for (unsigned int ifv=0; ifv<3; ifv++){
        if (ifv==0) fai1->setVal(0);
        else if (ifv==1) fai1->setVal(1);
        else if (ifv==2) fai1->setVal(0.5);

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH3F* tpl = (TH3F*)proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = 1;
        if (scale_width) integral = tpl->Integral("width");
        else integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          tpl=(TH3F*)proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      fai1->setVal(0);
      extractTemplates_ggLike_fai1(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
  }
  else if (proc.name.Contains("qqH") || proc.name.Contains("WH") || proc.name.Contains("ZH")){
    if (fai1!=0 && GGsm==0){
      vector<TH3F*> intpl;
      for (unsigned int ifv=0; ifv<5; ifv++){
        if (ifv==0) fai1->setVal(0);
        else if (ifv==1) fai1->setVal(1);
        else if (ifv==2) fai1->setVal(0.25);
        else if (ifv==3) fai1->setVal(0.5);
        else if (ifv==4) fai1->setVal(0.75);

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH3F* tpl = (TH3F*)proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm!=0) normval *= proc.norm->getVal();
        double integral = 1;
        if (scale_width) integral = tpl->Integral("width");
        else integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          tpl=(TH3F*)proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        }
        tpl->SetTitle("");
        tpl->Scale(scale);
        intpl.push_back(tpl);
      }
      fai1->setVal(0);
      extractTemplates_VVLike_fai1(tplname, intpl, templates);
      for (unsigned int it=0; it<intpl.size(); it++) delete intpl.at(it);
    }
  }

  proc.templates = templates;
}

void extractTemplates_ggLike_fai1(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl){
  double invA[3][3]={
    { 1, 0, 0 },
    { 0, 1, 0 },
    { -1, -1, 2 }
  };

  if (intpl.size()!=3) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname = Form("%s_Sig", newname.Data());
    if (ot>0) nname += Form("_ai1_%s", (ot==2 ? "1_Re" : "2"));
    outtpl.push_back((TH3F*)intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}
void extractTemplates_VVLike_fai1(const TString& newname, const vector<TH3F*>& intpl, vector<TH3F*>& outtpl){
  double A[25]={
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    9./16., 1./16., 3.*sqrt(3.)/16., 3./16., sqrt(3.)/16.,
    0.25, 0.25, 0.25, 0.25, 0.25,
    1./16., 9./16., sqrt(3.)/16., 3./16., 3.*sqrt(3.)/16.
  };
  TMatrixD matA(5, 5, A);
  TMatrixD invA = matA.Invert();

  if (intpl.size()!=5) return;
  for (unsigned int ot=0; ot<intpl.size(); ot++){
    TString nname = Form("%s_Sig", newname.Data());
    if (ot==1) nname += "_ai1_4";
    else if (ot>0) nname += Form("_ai1_%i%s", ot-1, (ot==2 || ot==4 ? "_Re" : "_PosDef"));
    outtpl.push_back((TH3F*)intpl.at(ot)->Clone(nname));
    outtpl.at(ot)->Reset("ICES");
    for (unsigned int it=0; it<intpl.size(); it++) outtpl.at(ot)->Add(intpl.at(it), invA[ot][it]);
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral() << endl;
  }
}

void splitOption(const string rawoption, string& wish, string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
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
