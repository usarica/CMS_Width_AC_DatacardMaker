#ifndef HIGGSCSANDWIDTHFERMI_CC
#define HIGGSCSANDWIDTHFERMI_CC


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

#include "HiggsCSandWidthFermi.h"

using namespace std;

HiggsCSandWidthFermi::HiggsCSandWidthFermi()
{

  N_BR = 681;

  ifstream file;
 
  // Read Widths into memory
  FileLoc = "include/txtFiles/Higgs_BR_Fermiophobic.txt"; //directory of input file
  const char* BranchRatioFileLoc = FileLoc.c_str(); 
  file.open(BranchRatioFileLoc);
  for(int k = 0; k < N_BR; k++){

    file >> mass_BR[k] >> BR[7][k] >> BR[8][k] >> BR[9][k] >> BR[10][k] >> BR[11][k] >> BR[0][k];


  }
  file.close();


}


HiggsCSandWidthFermi::~HiggsCSandWidthFermi()
{
  //destructor

}


// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidthFermi::HiggsWidth(int ID, double mH, bool spline){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /**************************************************/



  double TotalWidth = 0;
  double PartialWidth = 0;
  double Width = 0;
  int i = 0;
  double closestMass = 0;
  double tmpLow1, tmpHigh1, deltaX, deltaY1, slope1;
  double deltaY2, tmpLow2, tmpHigh2, slope2, step;


  // If ID is unavailable return -1                                           
  if((ID > 11 || ID < 8) && ID != 0){return 0;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 80 || mH > 250){return 0;}
  else{

    //Find index and closest higgs mass for which we have numbers
    step = 0.25; i = (int)((mH - 80)/step); closestMass = (double)(step*i + 80);
    
      tmpLow1 = BR[ID][i]*BR[0][i];                                                                                                                        
      tmpHigh1 = BR[ID][i+1]*BR[0][i+1];                                                                                                                   


      tmpLow2 = BR[0][i];
      tmpHigh2 = BR[0][i+1];
      deltaX = mH - closestMass;

      deltaY1 = tmpHigh1 - tmpLow1;
      slope1 = deltaY1/step;


      deltaY2 = tmpHigh2 - tmpLow2;
      slope2 = deltaY2/step;


      if(!spline)
	{
	  // For partial widths                                                                                                                 
	  if(deltaX == 0){ PartialWidth = tmpLow1;
	    TotalWidth = tmpLow2;}
	  else{ PartialWidth = slope1*deltaX + tmpLow1;
	    TotalWidth = slope2*deltaX + tmpLow2;}
	  // For total width  
	  if( ID == 0 ){ Width = TotalWidth; }
	  else{ Width = PartialWidth;}
	}
      else if(spline)
	{
	  if( ID == 0 )
	    {
	      if(i < 1){i = 1;}
	      if(i+2 >= N_BR){i = N_BR - 3;}

	      int indexWFF = 4;
	      double xmhWFF[indexWFF], sigWFF[indexWFF];
	      xmhWFF[0]=mass_BR[i-1];xmhWFF[1]=mass_BR[i];xmhWFF[2]=mass_BR[i+1];xmhWFF[3]=mass_BR[i+2];
	      sigWFF[0]=BR[ID][i-1]; sigWFF[1]=BR[ID][i]; sigWFF[2]=BR[ID][i+1]; sigWFF[3]=BR[ID][i+2];
	      //gROOT->Reset();
	      TGraph *graphWFF = new TGraph(indexWFF, xmhWFF, sigWFF);
	      TSpline3 *gsWFF = new TSpline3("gsW",graphWFF);
	      gsWFF->Draw();
	      Width = gsWFF->Eval(mH);
	      delete gsWFF;
	      delete graphWFF;
	    }
	  else{
	      if(i < 1){i = 1;}
	      if(i+2 >= N_BR){i = N_BR - 3;}

	      int indexWFF = 4;
	      double xmhWFF[indexWFF], sigWFF[indexWFF];
	      xmhWFF[0]=mass_BR[i-1];xmhWFF[1]=mass_BR[i];xmhWFF[2]=mass_BR[i+1];xmhWFF[3]=mass_BR[i+2];
	      sigWFF[0]=BR[0][i-1]; sigWFF[1]=BR[0][i]; sigWFF[2]=BR[0][i+1]; sigWFF[3]=BR[0][i+2];
	      //gROOT->Reset();
	      TGraph *graphWFF = new TGraph(indexWFF, xmhWFF, sigWFF);
	      TSpline3 *gsWFF = new TSpline3("gsWFF",graphWFF);
	      gsWFF->Draw();
	      PartialWidth = gsWFF->Eval(mH);
	      delete gsWFF;
	      delete graphWFF;
   
	      int indexPWFF = 4;
	      double xmhPWFF[indexPWFF], sigPWFF[indexPWFF];
	      xmhPWFF[0]=mass_BR[i-1];xmhPWFF[1]=mass_BR[i];xmhPWFF[2]=mass_BR[i+1];xmhPWFF[3]=mass_BR[i+2];
	      sigPWFF[0]=BR[ID][i-1]; sigPWFF[1]=BR[ID][i]; sigPWFF[2]=BR[ID][i+1]; sigPWFF[3]=BR[ID][i+2];
	      //gROOT->Reset();
	      TGraph *graphPWFF = new TGraph(indexPWFF, xmhPWFF, sigPWFF);
	      TSpline3 *gsPWFF = new TSpline3("gsPWFF",graphPWFF);
	      gsPWFF->Draw();
	      PartialWidth *= gsPWFF->Eval(mH);
	      delete gsPWFF;
	      delete graphPWFF;

	      Width = PartialWidth;
	  
	  }
	  
	}
      
  }
  
  return Width;
  
} 



double HiggsCSandWidthFermi::HiggsBR(int ID, double mH, bool spline){


  /***********************IDs************************/
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /**************************************************/




  double PartialBR = 0;
  double BranchRatio = 0;
  int i = 0;
  double closestMass = 0;
  double tmpLow1, tmpHigh1, deltaX, deltaY1, slope1;
  double step;


  // If ID is unavailable return -1                                           
  if((ID > 11 || ID < 8) && ID != 0){return 0;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 80 || mH > 250){return 0;}
  else{

    //Find index and closest higgs mass for which we have numbers
    step = 0.25; i = (int)((mH - 80)/step); closestMass = (double)(step*i + 80);
    
    tmpLow1 = BR[ID][i];
    tmpHigh1 = BR[ID][i+1];


    deltaX = mH - closestMass;

      deltaY1 = tmpHigh1 - tmpLow1;
      slope1 = deltaY1/step;


      if(!spline)
	{
	  if(deltaX == 0){ PartialBR = tmpLow1;}
	  else{ PartialBR = slope1*deltaX + tmpLow1;}
	}
      else if(spline)
	{
	  if(i < 1){i = 1;}
	  if(i+2 >= N_BR){i =N_BR - 3;}

	  int indexBRFF = 4;
	  double xmhBRFF[indexBRFF], sigBRFF[indexBRFF];
	  xmhBRFF[0]=mass_BR[i-1];xmhBRFF[1]=mass_BR[i];xmhBRFF[2]=mass_BR[i+1];xmhBRFF[3]=mass_BR[i+2];
	  sigBRFF[0]=BR[ID][i-1]; sigBRFF[1]=BR[ID][i]; sigBRFF[2]=BR[ID][i+1]; sigBRFF[3]=BR[ID][i+2];
	  //gROOT->Reset();
	  TGraph *graphBRFF = new TGraph(indexBRFF, xmhBRFF, sigBRFF);
	  TSpline3 *gsBRFF = new TSpline3("gsBRFF",graphBRFF);
	  gsBRFF->Draw();
	  PartialBR = gsBRFF->Eval(mH);
	  delete gsBRFF;
	  delete graphBRFF;
	}

      BranchRatio = PartialBR;

  }

  return BranchRatio;

} 





#endif
