#ifndef DJETCONSTRAINEDSYSTEMATIC_CC
#define DJETCONSTRAINEDSYSTEMATIC_CC


#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <fstream>

#include "TROOT.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"


#include "DjetConstrainedSystematic.h"

using namespace std;

DjetConstrainedSystematic::DjetConstrainedSystematic(string fileNominal, string fileQCDUp, string fileQCDDown, string filePDFUp, string filePDFDown, int useDjet, int ACfai1){
  DjetCat = useDjet;
  AnomCouplIndex = ACfai1;

  nACterms_ggZZ_signal = 1;
  nACterms_ggZZ_interf = 1;
  nACterms_VBF_signal = 1;
  nACterms_VBF_interf = 1;
  if (AnomCouplIndex==1){
    nACterms_ggZZ_signal = 3;
    nACterms_ggZZ_interf = 2;
    nACterms_VBF_signal = 5;
    nACterms_VBF_interf = 3;
  }

  for (int s=0; s<10; s++){
    for (int c=0; c<100; c++){
      djet_nondjet_ratio_signal[s][c]=0;
      djet_nondjet_ratio_interf[s][c]=0;
    }
    djet_nondjet_ratio_bkg[s]=0;
  }



  TString strDjet[2]={ "_nonDjet.root", "_Djet.root" };
  const int nMaxSyst=5;
  for (int s=0; s<nMaxSyst; s++){
    TString strSystFile;
    if (s==0) strSystFile = fileNominal;
    else if (s==1) strSystFile = fileQCDUp;
    else if (s==2) strSystFile = fileQCDDown;
    else if (s==3) strSystFile = filePDFUp;
    else if (s==4) strSystFile = filePDFDown;

    for (int j=0; j<2; j++){
      TString cinput = strSystFile;
      cinput.Append(strDjet[j]);
      TFile* finput=0;
      finput = new TFile(cinput, "read");
      if (finput==0 || finput->IsZombie()) continue;

      TH2F* htemp = 0;
      for (int al=0; al<nACterms_ggZZ_signal; al++){
        TString hname = "T_2D_1";
        if (AnomCouplIndex==1 && al>0) hname.Append(Form("_mZZ2_%i", al));
        else if (AnomCouplIndex!=0 && al>0) hname.Append(Form("_ACTerm_%i", al));
        htemp = (TH2F*)finput->Get(hname);
        if (j==0) djet_nondjet_ratio_signal[s][al] = htemp->Integral("width");
        else if (djet_nondjet_ratio_signal[s][al]!=0) djet_nondjet_ratio_signal[s][al] = htemp->Integral("width") / djet_nondjet_ratio_signal[s][al];
        else{
          cerr << "Invalid ggZZ signal integral!!!" << endl; djet_nondjet_ratio_signal[s][al]=0;
        }
        if (htemp!=0) delete htemp;
      }
      for (int al=0; al<nACterms_ggZZ_interf; al++){
        TString hname = "T_2D_4";
        if (AnomCouplIndex==1 && al>0) hname.Append(Form("_mZZ2_%i", al));
        else if (AnomCouplIndex!=0 && al>0) hname.Append(Form("_ACTerm_%i", al));
        htemp = (TH2F*)finput->Get(hname);
        if (j==0) djet_nondjet_ratio_interf[s][al] = htemp->Integral("width");
        else if (djet_nondjet_ratio_interf[s][al]!=0) djet_nondjet_ratio_interf[s][al] = htemp->Integral("width") / djet_nondjet_ratio_interf[s][al];
        else{
          cerr << "Invalid ggZZ interf integral!!!" << endl; djet_nondjet_ratio_interf[s][al]=0;
        }
        if (htemp!=0) delete htemp;
      }
      for (int al=0; al<nACterms_VBF_signal; al++){
        TString hname = "T_2D_VBF_1";
        if (AnomCouplIndex==1 && al>0) hname.Append(Form("_mZZ2_%i", al));
        else if (AnomCouplIndex!=0 && al>0) hname.Append(Form("_ACTerm_%i", al));
        htemp = (TH2F*)finput->Get(hname);
        if (j==0) djet_nondjet_ratio_signal[s+nMaxSyst][al] = htemp->Integral("width");
        else if (djet_nondjet_ratio_signal[s+nMaxSyst][al]!=0) djet_nondjet_ratio_signal[s+nMaxSyst][al] = htemp->Integral("width") / djet_nondjet_ratio_signal[s+nMaxSyst][al];
        else{
          cerr << "Invalid VBF signal integral!!!" << endl; djet_nondjet_ratio_signal[s+nMaxSyst][al]=0;
        }
        if (htemp!=0) delete htemp;
      }
      for (int al=0; al<nACterms_VBF_interf; al++){
        TString hname = "T_2D_VBF_4";
        if (AnomCouplIndex==1 && al>0) hname.Append(Form("_mZZ2_%i", al));
        else if (AnomCouplIndex!=0 && al>0) hname.Append(Form("_ACTerm_%i", al));
        htemp = (TH2F*)finput->Get(hname);
        if (j==0) djet_nondjet_ratio_interf[s+nMaxSyst][al] = htemp->Integral("width");
        else if (djet_nondjet_ratio_interf[s+nMaxSyst][al]!=0) djet_nondjet_ratio_interf[s+nMaxSyst][al] = htemp->Integral("width") / djet_nondjet_ratio_interf[s+nMaxSyst][al];
        else{
          cerr << "Invalid VBF interf integral!!!" << endl; djet_nondjet_ratio_interf[s+nMaxSyst][al]=0;
        }
        if (htemp!=0) delete htemp;
      }
      htemp = (TH2F*)finput->Get("T_2D_2");
      if (j==0) djet_nondjet_ratio_bkg[s] = htemp->Integral("width");
      else if (djet_nondjet_ratio_bkg[s]!=0) djet_nondjet_ratio_bkg[s] = htemp->Integral("width") / djet_nondjet_ratio_bkg[s];
      else{
        cerr << "Invalid ggZZ bkg integral!!!" << endl; djet_nondjet_ratio_bkg[s]=0;
      }
      if (htemp!=0) delete htemp;

      htemp = (TH2F*)finput->Get("T_2D_VBF_2");
      if (j==0) djet_nondjet_ratio_bkg[s+nMaxSyst] = htemp->Integral("width");
      else if (djet_nondjet_ratio_bkg[s+nMaxSyst]!=0) djet_nondjet_ratio_bkg[s+nMaxSyst] = htemp->Integral("width") / djet_nondjet_ratio_bkg[s+nMaxSyst];
      else{
        cerr << "Invalid VBF bkg integral!!!" << endl; djet_nondjet_ratio_bkg[s+nMaxSyst]=0;
      }
      if (htemp!=0) delete htemp;

      finput->Close();
    }
  }
}


DjetConstrainedSystematic::~DjetConstrainedSystematic()
{
  //destructor

}


double DjetConstrainedSystematic::getConstrainedRatio(string GGsmAmplitude, string systematic, string coupling){
  int iSyst, iCoupl;
  int GGsmTerm = matchStringToIndex(GGsmAmplitude, systematic, coupling, iSyst, iCoupl);

  if(iSyst<0){
    cerr << "iSyst DOES NOT EXIST!!!" << endl; return 0;
  }

  if (DjetCat==0) return 0;
  else if (DjetCat==2) return -1;

  if (GGsmTerm==1) return djet_nondjet_ratio_signal[iSyst][iCoupl];
  else if (GGsmTerm==2) return djet_nondjet_ratio_bkg[iSyst];
  else if (GGsmTerm==4) return djet_nondjet_ratio_interf[iSyst][iCoupl];
  else{
    cerr << "Failed to find the item requested!!!" << endl; return 0;
  }
}


int DjetConstrainedSystematic::matchStringToIndex(string GGsmAmplitude, string systematic, string coupling, int& iSyst, int& iCoupl){
  int GGsmTerm = -1;
  iSyst = -1;
  iCoupl = -1;

  if (GGsmAmplitude == "signal") GGsmTerm = 1;
  else if (GGsmAmplitude == "interf") GGsmTerm = 4;
  else if (GGsmAmplitude == "bkg") GGsmTerm = 2;

  if (systematic=="ggZZNominal") iSyst = 0;
  else if (systematic=="ggZZQCDUp") iSyst = 1;
  else if (systematic=="ggZZQCDDown") iSyst = 2;
  else if (systematic=="ggZZPDFUp") iSyst = 3;
  else if (systematic=="ggZZPDFDown") iSyst = 4;
  else if (systematic=="VBFNominal") iSyst = 5;
  else if (systematic=="VBFUp") iSyst = 8; // Receives up from PDFUp file
  else if (systematic=="VBFDown") iSyst = 9; // Receives down from PDFDown file

  if (AnomCouplIndex==0) iCoupl=0;
  else if (coupling=="mZZ2_1" || coupling=="ACTerm_1") iCoupl = 1;
  else if (coupling=="mZZ2_2" || coupling=="ACTerm_2") iCoupl = 2;
  else if (coupling=="mZZ2_3" || coupling=="ACTerm_3") iCoupl = 3;
  else if (coupling=="mZZ2_4" || coupling=="ACTerm_4") iCoupl = 4;

  return GGsmTerm;
}

#endif
