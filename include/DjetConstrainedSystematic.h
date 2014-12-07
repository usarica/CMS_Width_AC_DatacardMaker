#ifndef DJETCONSTRAINEDSYSTEMATIC_H
#define DJETCONSTRAINEDSYSTEMATIC_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

using namespace std;

class DjetConstrainedSystematic
{

public:

  DjetConstrainedSystematic(string fileNominal, string fileQCDUp, string fileQCDDown, string filePDFUp, string filePDFDown, int useDjet, int ACfai1);
  ~DjetConstrainedSystematic();

  double getConstrainedRatio(string GGsmAmplitude, string systematic, string coupling);

private:

  int DjetCat;
  int AnomCouplIndex;

  int nACterms_ggZZ_signal;
  int nACterms_ggZZ_interf;
  int nACterms_VBF_signal;
  int nACterms_VBF_interf;

  double djet_nondjet_ratio_signal[10][100];
  double djet_nondjet_ratio_interf[10][100];
  double djet_nondjet_ratio_bkg[10];

  int matchStringToIndex(string GGsmAmplitude, string systematic, string coupling, int& iSyst, int& iCoupl);
};

#endif
