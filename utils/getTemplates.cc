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
  process_spec(const process_spec& other, TString newname="") : pdf(other.pdf), pdf_shape(other.pdf_shape), norm(other.norm), rate(other.rate), name(newname=="" ? other.name: newname), templates(other.templates), systematics(other.systematics){}

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
TString getSystRename(TString const& systname, TString const& strSqrts, TString const& strPeriod){
  TString res=systname;
  if (res.Contains("lumi")) res = "lumiUnc";
  else if (res == "pdf_qq") res = "pdf_qqbar";
  else if (res == "QCDscale_ggVV_bonly") res = "kbkg_gg";
  else if (res == "CMS_zz4mu_zjets") res = "CMS_hzz4l_zz4mu_zjets";
  else if (res == "CMS_zz4e_zjets") res = "CMS_hzz4l_zz4e_zjets";
  else if (res == "CMS_zz2e2mu_zjets") res = "CMS_hzz4l_zz2e2mu_zjets";
  else if (res == "Res4mu") res = "CMS_hzz4l_zz4mu_res";
  else if (res == "Res4e") res = "CMS_hzz4l_zz4e_res";
  else if (res == "Res2e2mu") res = "CMS_hzz4l_zz2e2mu_res";
  else if (res == "CMS_zz4l_smd_zjets_bkg_4mu") res = "CMS_hzz4l_zz4mu_shape_zjets";
  else if (res == "CMS_zz4l_smd_zjets_bkg_4e") res = "CMS_hzz4l_zz4e_shape_zjets";
  else if (res == "CMS_zz4l_smd_zjets_bkg_2e2mu") res = "CMS_hzz4l_zz2e2mu_shape_zjets";
  return res;
}
float getProcessRescale(TString const& procname, TString const& strSqrts, TString const& strPeriod){
  float res=1;
  if (strSqrts=="13TeV" && strPeriod=="2015"){
    float ggH[2] ={ 1.22e-2, 8.491e-3 };
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
    else if (procname.Contains("VVH")) res = (VBFH[1]+ZH[1]+WH[1])/(VBFH[0]+ZH[0]+WH[0]);
  }
  return res;
}


void getDataTree(TString cinput){
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
  string catnames[ncats]={ "Inclusive", "VBFtagged", "VHHadrtagged", "Untagged" };
  string catname=catnames[0];
  for (unsigned int ic=0; ic<ncats; ic++){
    if (strinput.find(catnames[ic])!=string::npos) catname=catnames[ic];
  }

  TFile* finput = TFile::Open(cinput+".input.root", "read");
  RooWorkspace* ws = (RooWorkspace*)finput->Get("w");
  RooDataSet* data = (RooDataSet*)ws->data("data_obs");
  int nevents=data->sumEntries();
  int nvars=((const RooArgSet*)data->get())->getSize();
  float* KD = new float[nvars];
  float mass;

  TString coutput_data = "Decompilation/Data";
  gSystem->Exec("mkdir -p " + coutput_data);
  TString coutput_root;
  TFile* foutput;
  coutput_root = Form("%s/hzz%s_%s_%s.root", coutput_data.Data(), channame.c_str(), catname.c_str(), strSqrtsPeriod.Data());
  foutput = TFile::Open(coutput_root, "recreate");

  TTree* t = new TTree("data_obs", "");
  //for (int iv=0; iv<nvars; iv++) t->Branch(Form("KD%i", iv+1), KD+iv);

  unsigned int nKDs=0;
  for (int ev=0; ev<nevents; ev++){
    const RooArgSet* args = (const RooArgSet*)data->get(ev);
    TIterator* coefIter = args->createIterator();
    const RooAbsArg* coef;
    unsigned int ik=0;
    unsigned int ikd=0;
    while ((coef = (const RooAbsArg*)coefIter->Next())){
      const RooAbsReal* rvar = dynamic_cast<const RooAbsReal*>(coef);
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
void getTemplates(TString cinput, double lumiScale=1, bool scale_width=true, bool copy_ggH_to_VVH=false, bool rescale_xsec=false){
  TString strSqrtsPeriod, strSqrts, strPeriod;
  getSqrtsPeriod(cinput, strSqrtsPeriod, strSqrts, strPeriod);

  getDataTree(cinput);

  string strinput = cinput.Data();
  vector<string> splitinput;
  splitOptionRecursive(strinput, splitinput, '/');
  strinput = splitinput.back(); splitinput.clear();
  splitOptionRecursive(strinput, splitinput, '/');

  TString coutput_templates = "Decompilation/Templates";
  TString coutput_inputs = "Decompilation/Inputs";

  gSystem->Exec("mkdir -p " + coutput_templates);
  gSystem->Exec("mkdir -p " + coutput_inputs);

  const unsigned int nchans=3;
  string channames[nchans]={ "4mu", "4e", "2e2mu" };
  string channame;
  for (unsigned int ic=0; ic<nchans; ic++){
    if (strinput.find(channames[ic])!=string::npos) channame=channames[ic];
  }

  const unsigned int ncats=4;
  string catnames[ncats]={ "Inclusive", "JJVBFTagged", "HadVHTagged", "Untagged" };
  string catname = catnames[0];
  for (unsigned int ic=0; ic<ncats; ic++){
    if (strinput.find(catnames[ic])!=string::npos) catname=catnames[ic];
  }

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
    RooAbsReal* norm = 0; norm = (RooAbsReal*)ws->factory(procname.at(ip)+"_norm");
    if (norm==0) cout << "Warning: " << procname.at(ip) << "_norm is not found." << endl;

    if (pdf==0) cerr << procname.at(ip) << " pdf could not be found." << endl;
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
          TIterator* clientsIter = systvar->clientIterator();
          RooAbsArg* client;
          while ((client = (RooAbsArg*)clientsIter->Next())){
            for (unsigned int ip=0; ip<procname.size(); ip++){
              RooAbsPdf* pdf = procSpecs[procname.at(ip).Data()].pdf;
              if (client->GetName() == pdf->GetName()){
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
  for (unsigned int ip=0; ip<procname.size(); ip++) tout << "channel " << procname.at(ip) << " 1 -1 " << (procSpecs[procname.at(ip).Data()].name.Contains("bkg") ? 0 : 1) << endl;

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
    systName = getSystRename(systName, strSqrts, strPeriod);
    tout << "systematic " << systName << " template " << syst->second << endl;
  }
  for (auto syst = logSyst.begin(); syst != logSyst.end(); ++syst){
    TString systName = syst->first.c_str();
    systName = getSystRename(systName, strSqrts, strPeriod);
    tout << "systematic " << systName << " lnN " << syst->second << endl;
  }
  for (auto syst = paramSyst.begin(); syst != paramSyst.end(); ++syst){
    TString systName = syst->first.c_str();
    systName = getSystRename(systName, strSqrts, strPeriod);
    tout << "systematic " << systName << " template " << syst->second << endl;
  }

  for (unsigned int ip=0; ip<procname.size(); ip++){
    cout << "Attempting to extract tpls for process " << procname.at(ip) << endl;

    TString const& theProcName = procname.at(ip);
    TString theProcNameLower = theProcName; theProcNameLower.ToLower();
    float tplscale=1;
    if (!theProcNameLower.Contains("zjets")) tplscale=1./lumiScale;
    if (rescale_xsec) tplscale *= getProcessRescale(theProcName, strSqrts, strPeriod);

    TFile* foutput;
    TString coutput_root;
    
    coutput_root = Form("%s/HtoZZ%s_%s_FinalTemplates_%s_%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procname.at(ip).Data(), "Nominal", ".root");
    foutput = TFile::Open(coutput_root, "recreate");
    extractTemplates(procSpecs[procname.at(ip).Data()], data, "", scale_width);
    for (unsigned int it=0; it<procSpecs[procname.at(ip).Data()].templates.size(); it++) procSpecs[procname.at(ip).Data()].templates.at(it)->Scale(tplscale);
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
          systvar->setVal(double(2*is)-1);
          cout << "Setting param systematic " << systvar->GetName() << " to " << double(2*is)-1 << endl;
          TString systName = syst->first.c_str();
          systName = getSystRename(systName, strSqrts, strPeriod);
          coutput_root = Form("%s/HtoZZ%s_%s_FinalTemplates_%s_%s%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procname.at(ip).Data(), systName.Data(), (is==0 ? "Down" : "Up"), ".root");
          foutput = TFile::Open(coutput_root, "recreate");
          extractTemplates(procSpecs[procname.at(ip).Data()], data, "", scale_width);
          for (unsigned int it=0; it<procSpecs[procname.at(ip).Data()].templates.size(); it++) procSpecs[procname.at(ip).Data()].templates.at(it)->Scale(tplscale);
          procSpecs[procname.at(ip).Data()].writeTemplates(foutput);
          procSpecs[procname.at(ip).Data()].deleteTemplates();
          foutput->Close();
        }
      }
      else if (syst->second.first=="shape1"){
        for (unsigned int is=0; is<2; is++){
          string syst_du = syst->first + (is==0 ? "Down" : "Up");
          TString systName = syst->first.c_str();
          systName = getSystRename(systName, strSqrts, strPeriod);
          coutput_root = Form("%s/HtoZZ%s_%s_FinalTemplates_%s_%s%s%s", coutput_templates.Data(), channame.c_str(), catname.c_str(), procname.at(ip).Data(), systName.Data(), (is==0 ? "Down" : "Up"), ".root");
          foutput = TFile::Open(coutput_root, "recreate");
          extractTemplates(procSpecs[procname.at(ip).Data()], data, syst_du, scale_width);
          for (unsigned int it=0; it<procSpecs[procname.at(ip).Data()].templates.size(); it++) procSpecs[procname.at(ip).Data()].templates.at(it)->Scale(tplscale);
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
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral("width") << endl;
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
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral("width") << endl;
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
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral("width") << endl;
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
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral("width") << endl;
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
    cout << outtpl.at(ot)->GetName() << " integral = " << outtpl.at(ot)->Integral("width") << endl;
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

void extractTemplates(process_spec& proc, RooDataSet* data, string shapename, bool scale_width){
  vector<TH3F*> templates;

  vector<RooRealVar*> deps;
  if (!proc.pdf){
    cerr << "ERROR: PDF IS NULL!" << endl;
    assert(0);
  }
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
  else if (proc.name.Contains("gg") || proc.name.Contains("tt")){
    float RVdefval=0;
    if (RV){
      cout << "\t- Setting RV to 0" << endl;
      RVdefval = RV->getVal();
      RV->setVal(0);
    }
    else{
      cout << "Could not find RV" << endl;
    }
    if (fai1 && !GGsm){
      vector<TH3F*> intpl;
      for (unsigned int ifv=0; ifv<3; ifv++){
        fai1->setVal(0.5*float(ifv));

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH3F* tpl = (TH3F*)proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
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
    else if (!fai1 && GGsm){
      vector<TH3F*> intpl;
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
        TH3F* tpl = (TH3F*) proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
        double integral = 1;
        if (scale_width) integral = tpl->Integral("width");
        else integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          tpl=(TH3F*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
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
      vector<TH3F*> intpl;
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
        TH3F* tpl = (TH3F*) proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
        double integral = 1;
        if (scale_width) integral = tpl->Integral("width");
        else integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          tpl=(TH3F*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
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
  else if (proc.name.Contains("qqH") || proc.name.Contains("VBF") || proc.name.Contains("WH") || proc.name.Contains("ZH") || proc.name.Contains("VV")){
    float RFdefval=0;
    if (RF){
      cout << "\t- Setting RF to 0" << endl;
      RFdefval = RF->getVal();
      RF->setVal(0);
    }
    else{
      cout << "Could not find RF" << endl;
    }
    if (fai1 && !GGsm){
      vector<TH3F*> intpl;
      for (unsigned int ifv=0; ifv<5; ifv++){
        fai1->setVal(0.25*float(ifv));

        TString theName = Form("%s_%i", tplname.Data(), ifv);
        TH3F* tpl = (TH3F*)proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
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
    else if (!fai1 && GGsm){
      vector<TH3F*> intpl;
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
        TH3F* tpl = (TH3F*) proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
        double integral = 1;
        if (scale_width) integral = tpl->Integral("width");
        else integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          tpl=(TH3F*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
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
      vector<TH3F*> intpl;
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
        TH3F* tpl = (TH3F*) proc.pdf->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
        double normval = proc.rate; if (proc.norm) normval *= proc.norm->getVal();
        double integral = 1;
        if (scale_width) integral = tpl->Integral("width");
        else integral = tpl->Integral();
        double scale = normval/integral;
        cout << "Scaling template " << tplname << " by " << normval << " / " << integral << endl;
        if (shapename!=""){
          delete tpl;
          tpl=(TH3F*) proc.pdf_shape[shapename]->createHistogram(theName, *(deps.at(0)), ycmd, zcmd);
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

  proc.templates = templates;
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
