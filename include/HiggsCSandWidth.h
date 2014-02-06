#ifndef HIGGSCSANDWIDTH_H
#define HIGGSCSANDWIDTH_H

#define PI 3.14159

#define  ID_ggToH  1
#define  ID_VBF    2
#define  ID_WH     3
#define  ID_ZH     4
#define  ID_ttH    5
#define  ID_Total  0 

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"


/**********************************************************/
/*            Class for Higgs Width and CS                */
/*                                                        */
/*  All numbers for CS and width are taken from official  */
/*  numbers on Higgs CS Twiki (Spring 2011)               */
/*                                                        */
/*  Cross Sections are given in pb                        */
/*  Widths are given in GeV                               */
/*                                                        */
/*  These numbers are taken into memory and a simple      */
/*  linear interpolation is done.                         */
/*                                                        */
/*  For any invalid process or mH out of range, -1 will   */
/*  be returned.                                          */
/*                                                        */
/*    Written by:                                         */
/*         Matt Snowball                                  */
/*         University of Florida                          */
/*         snowball@phys.ufl.edu                          */
/*                                                        */
/*       Last Update: April 5, 2012                       */
/*                                                        */
/**********************************************************/



class HiggsCSandWidth
{

 public:

  HiggsCSandWidth(std::string fileLoc = "include/txtFiles");
  ~HiggsCSandWidth();

  double HiggsCS(int ID, double mH, double sqrts);
  double HiggsCSErrPlus(int ID, double mH, double sqrts);
  double HiggsCSErrMinus(int ID, double mH, double sqrts);
  double HiggsCSscaleErrPlus(int ID, double mH, double sqrts);
  double HiggsCSscaleErrMinus(int ID, double mH, double sqrts);
  double HiggsCSpdfErrPlus(int ID, double mH, double sqrts);
  double HiggsCSpdfErrMinus(int ID, double mH, double sqrts);
  double HiggsWidth(int ID,double mH);
  double HiggsBR(int ID,double mH);
  double getInterpXS(int sqrts, int ID, double mH, int maxI, double mhArray[][6], double varArray[][6]);

 private:

  std::string fileName;
  
  double scratchMass;
  double mass_BR[217];
  double BR[26][217];

  double mass_XS_7tev[197][6];
  double CS_7tev[197][6];
  double CSerrPlus_7tev[197][6];
  double CSerrMinus_7tev[197][6];
  double CSscaleErrPlus_7tev[197][6];
  double CSscaleErrMinus_7tev[197][6];
  double CSpdfErrPlus_7tev[197][6];
  double CSpdfErrMinus_7tev[197][6];

  double mass_XS_8tev[223][6];
  double CS_8tev[223][6];
  double CSerrPlus_8tev[223][6];
  double CSerrMinus_8tev[223][6];
  double CSscaleErrPlus_8tev[223][6];
  double CSscaleErrMinus_8tev[223][6];
  double CSpdfErrPlus_8tev[223][6];
  double CSpdfErrMinus_8tev[223][6];

  double mass_XS_14tev[223][6];
  double CS_14tev[223][6];
  double CSerrPlus_14tev[223][6];
  double CSerrMinus_14tev[223][6];
  double CSscaleErrPlus_14tev[223][6];
  double CSscaleErrMinus_14tev[223][6];
  double CSpdfErrPlus_14tev[223][6];
  double CSpdfErrMinus_14tev[223][6];

  int N_BR;
  int N_CS_7tev[6];  
  int N_CS_8tev[6];  
  int N_CS_14tev[6];  



};

#endif
