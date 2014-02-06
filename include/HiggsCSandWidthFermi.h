#ifndef HIGGSCSANDWIDTHFERMI_H
#define HIGGSCSANDWIDTHFERMI_H

#define PI 3.14159

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>


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



class HiggsCSandWidthFermi
{

 public:

  HiggsCSandWidthFermi();
  ~HiggsCSandWidthFermi();

  double HiggsWidth(int ID,double mH,bool spline);
  double HiggsBR(int ID,double mH,bool spline);


 private:

  double BR[12][681];
  double mass_BR[681];

  int N_BR;
  
  std::string FileLoc;


};

#endif
