#ifndef HIGGSCSANDWIDTHSM4_CC
#define HIGGSCSANDWIDTHSM4_CC


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

#include "HiggsCSandWidthSM4.h"

using namespace std;

HiggsCSandWidthSM4::HiggsCSandWidthSM4()
{

  N_BR = 102;
  N_CS = 109;
  N_CSE = 175;
  N_BR_gg = 6;

  ifstream file;
  double scratchHgg;
  // Read Widths into memory
  file.open("include/txtFiles/Higgs_BR_SM4.txt");
  for(int k = 0; k < N_BR; k++){

    file >> mass_BR[k] >> BR[0][k] >> BR[1][k] >> BR[2][k] >> BR[3][k] >> BR[4][k] >> BR[5][k] >> BR[6][k] >> BR[7][k] >> scratchHgg >> BR[9][k]
	 >> BR[10][k] >> BR[11][k] >> BR[12][k] >> BR[13][k] >> BR[14][k] >> BR[15][k] >> BR[16][k] >> BR[17][k];


  }
  file.close();

  // Read Widths into memory
  file.open("include/txtFiles/Higgs_BR_SM4_Hgg.txt");
  for(int k = 0; k < N_BR_gg; k++){

    file >> mass_BR_gg[k] >> BR[8][k];


  }
  file.close();

  // Read CS into memory
  file.open("include/txtFiles/HiggsCS_Official_SM4.txt");//directory of input file
  for(int k = 0; k < N_CS; k++){

    file >> mass_XS[k] >> CS[ID_ggToH][k];// >> CS[ID_VBF][k] >> CS[ID_WH][k] >> CS[ID_ZH][k] >> CS[ID_ttH][k] >> CS[ID_Total][k];

    //cout << scratchMass << "  " << CS[ID_ggToH][k] << endl;

  }
  file.close();

  file.open("include/txtFiles/HiggsCS_Error_Official_SM4.txt");//directory of input file                       
  for(int k = 0; k < N_CSE; k++){

    file >> scratchMass >> CSerrPlus[ID_ggToH][k] >> CSerrMinus[ID_ggToH][k] >> CSscaleErrPlus[ID_ggToH][k] >> CSscaleErrMinus[ID_ggToH][k]
	 >> CSpdfErrPlus[ID_ggToH][k] >> CSpdfErrMinus[ID_ggToH][k];

  }
  file.close();


}


HiggsCSandWidthSM4::~HiggsCSandWidthSM4()
{
  //destructor

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)
double HiggsCSandWidthSM4::HiggsCS(int ID, double mH, double sqrts, bool spline){

  /**********IDs*************/ 
  /*     ggToH = 1          */
  /*       VBF = 2          */ 
  /*        WH = 3          */ 
  /*        ZH = 4          */
  /*       ttH = 5          */
  /*     Total = 0          */
  /**************************/
 
  int i = 0;
  double closestMass = 0;
  double reqCS = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                                                                                                
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}
  // If Ecm is not 7 TeV return -1
  if(sqrts != 7){return -1;}
 

  // If mH is out of range return -1                                           
  // else find what array number to read         
  if( mH < 100 || mH > 1000){return 0;}
  else{


    //Find index and closest higgs mass for which we have numbers
    if(mH <= 140){step = 10; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 140 && mH <= 150 ){step = 5; i = (int)(4 + (mH-140)/step); closestMass = (step*(i-4) + 140);}
    if(mH > 150 && mH <= 190 ){step = 1; i = (int)(6 + (mH-150)/step); closestMass = (int)(step*(i-6) + 150);}
    if(mH > 190 && mH <= 200 ){step = 5; i = (int)(46 + (mH-190)/step); closestMass = (int)(step*(i-46) + 190);}
    if(mH > 200 && mH <= 330 ){step = 10; i = (int)(48 + (mH-200)/step); closestMass = (int)(step*(i-48) + 200);}
    if(mH > 330 && mH <= 335 ){step = 5; i = (int)(61 + (mH-330)/step); closestMass = (int)(step*(i-61) + 330);}
    if(mH > 335 && mH <= 355 ){step = 1; i = (int)(62 + (mH-335)/step); closestMass = (int)(step*(i-62) + 335);}
    if(mH > 355 && mH <= 360 ){step = 5; i = (int)(82 + (mH-355)/step); closestMass = (int)(step*(i-82) + 355);}
    if(mH > 360 && mH <= 370 ){step = 10; i = (int)(83 + (mH-360)/step); closestMass = (int)(step*(i-83) + 360);}
    if(mH > 370 && mH <= 380 ){step = 5; i = (int)(84 + (mH-370)/step); closestMass = (int)(step*(i-84) + 370);}
    if(mH > 380 && mH <= 500 ){step = 10; i = (int)(86 + (mH-380)/step); closestMass = (int)(step*(i-86) + 380);}
    if(mH > 500 && mH <= 1000 ){step = 50; i = (int)(98 + (mH-500)/step); closestMass = (int)(step*(i-98) + 500);}


      tmpLow = CS[ID][i];
      tmpHigh = CS[ID][i+1];

      deltaX = mH - closestMass;

      deltaY = tmpHigh - tmpLow;
      slope = deltaY/step;

      if(!spline)
	{
	  if(deltaX == 0){ reqCS = tmpLow;}
	  else{ reqCS = slope*deltaX + tmpLow;}
	}
      else if(spline)
	{
	  if(i < 1){i = 1;}
	  if(i+2 >= N_CS){i = N_CS - 3;}
	  
	  int indexSM4 = 4;
	  double xmhSM4[indexSM4], sigSM4[indexSM4];
	  xmhSM4[0]=mass_XS[i-1];xmhSM4[1]=mass_XS[i];xmhSM4[2]=mass_XS[i+1];xmhSM4[3]=mass_XS[i+2];
	  sigSM4[0]=CS[ID][i-1]; sigSM4[1]=CS[ID][i]; sigSM4[2]=CS[ID][i+1]; sigSM4[3]=CS[ID][i+2];
          
	  TGraph *graphSM4 = new TGraph(indexSM4, xmhSM4, sigSM4);
	  TSpline3 *gsSM4 = new TSpline3("gsSM4",graphSM4);
	  gsSM4->Draw();
	  reqCS = gsSM4->Eval(mH);
	  delete gsSM4;
	  delete graphSM4;
	}
     

    }

  return reqCS;

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)                   
double HiggsCSandWidthSM4::HiggsCSErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  int i = 0;
  double closestMass = 0;
  double reqCSerr = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                                                                                    
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}
  // If Ecm is not 7 TeV return -1                                                                                                
  if(sqrts != 7){return -1;}
 
  // If mH is out of range return -1                                                                        
  // else find what array number to read                                          
  if( mH < 100 || mH > 600){return 0;}
  else{

    if(mH <= 110 ){step = 5; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(2 + (mH - 110)/step); closestMass = (step*(i-2) + 110);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(62 + (mH - 140)/step); closestMass = (int)(step*(i-62) + 140);}
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(82 + (mH - 160)/step); closestMass = (int)(step*(i-82) + 160);}
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(147 + (mH - 290)/step); closestMass = (int)(step*(i-147) + 290);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(159 + (mH-350)/step); closestMass = (int)(step*(i-159) + 350);}
    if(mH > 400){step = 20; i = (int)(164 + (mH-400)/step); closestMass = (int)(step*(i-164) + 400);}




    tmpLow = CSerrPlus[ID][i];
    tmpHigh = CSerrPlus[ID][i+1];

    deltaX = mH - closestMass;

    deltaY = tmpHigh - tmpLow;
    slope = deltaY/step;

    // For partial widths                                                                                                      
    if(deltaX == 0){ reqCSerr = tmpLow;}
    else{ reqCSerr = slope*deltaX + tmpLow;}
    reqCSerr *= .01; //Account for percentage
  }

  return reqCSerr;

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)      
double HiggsCSandWidthSM4::HiggsCSErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  int i = 0;
  double closestMass = 0;
  double reqCSerr = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                                                                                       
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}

  // If Ecm is not 7 TeV return -1                                                                                           
  if(sqrts != 7){return -1;}
  
  // If mH is out of range return -1                                                                           
  // else find what array number to read                                                                 
  if( mH < 100 || mH > 600){return 0;}
  else{

    if(mH <= 110 ){step = 5; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(2 + (mH - 110)/step); closestMass = (step*(i-2) + 110);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(62 + (mH - 140)/step); closestMass = (int)(step*(i-62) + 140);}
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(82 + (mH - 160)/step); closestMass = (int)(step*(i-82) + 160);}
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(147 + (mH - 290)/step); closestMass = (int)(step*(i-147) + 290);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(159 + (mH-350)/step); closestMass = (int)(step*(i-159) + 350);}
    if(mH > 400){step = 20; i = (int)(164 + (mH-400)/step); closestMass = (int)(step*(i-164) + 400);}




    tmpLow = CSerrMinus[ID][i];
    tmpHigh = CSerrMinus[ID][i+1];

    deltaX = mH - closestMass;

    deltaY = tmpHigh - tmpLow;
    slope = deltaY/step;

    // For partial widths                                               
    if(deltaX == 0){ reqCSerr = tmpLow;}
    else{ reqCSerr = slope*deltaX + tmpLow;}
    reqCSerr *= .01; //Account for percentage                                                                                                   
  }

  return reqCSerr;

}

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)          
double HiggsCSandWidthSM4::HiggsCSscaleErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  int i = 0;
  double closestMass = 0;
  double reqCSerr = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                                                         
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}

  // If Ecm is not 7 TeV return -1                                                
  if(sqrts != 7){return -1;}

  // If mH is out of range return -1                                                         
  // else find what array number to read                                                      
  if( mH < 100 || mH > 600){return 0;}
  else{

    if(mH <= 110 ){step = 5; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(2 + (mH - 110)/step); closestMass = (step*(i-2) + 110);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(62 + (mH - 140)/step); closestMass = (int)(step*(i-62) + 140);}
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(82 + (mH - 160)/step); closestMass = (int)(step*(i-82) + 160);}
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(147 + (mH - 290)/step); closestMass = (int)(step*(i-147) + 290);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(159 + (mH-350)/step); closestMass = (int)(step*(i-159) + 350);}
    if(mH > 400){step = 20; i = (int)(164 + (mH-400)/step); closestMass = (int)(step*(i-164) + 400);}




    tmpLow = CSscaleErrPlus[ID][i];
    tmpHigh = CSscaleErrPlus[ID][i+1];

    deltaX = mH - closestMass;

    deltaY = tmpHigh - tmpLow;
    slope = deltaY/step;

    // For partial widths                                              
    if(deltaX == 0){ reqCSerr = tmpLow;}
    else{ reqCSerr = slope*deltaX + tmpLow;}
    reqCSerr *= .01; //Account for percentage                 
  }

  return reqCSerr;

}

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)    
double HiggsCSandWidthSM4::HiggsCSscaleErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  int i = 0;
  double closestMass = 0;
  double reqCSerr = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                     
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}

  // If Ecm is not 7 TeV return -1                                                               
  if(sqrts != 7){return -1;}
 
  // If mH is out of range return -1                        
  // else find what array number to read                              
  if( mH < 100 || mH > 600){return 0;}
  else{

    if(mH <= 110 ){step = 5; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(2 + (mH - 110)/step); closestMass = (step*(i-2) + 110);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(62 + (mH - 140)/step); closestMass = (int)(step*(i-62) + 140);}
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(82 + (mH - 160)/step); closestMass = (int)(step*(i-82) + 160);}
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(147 + (mH - 290)/step); closestMass = (int)(step*(i-147) + 290);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(159 + (mH-350)/step); closestMass = (int)(step*(i-159) + 350);}
    if(mH > 400){step = 20; i = (int)(164 + (mH-400)/step); closestMass = (int)(step*(i-164) + 400);}



    tmpLow = CSscaleErrMinus[ID][i];
    tmpHigh = CSscaleErrMinus[ID][i+1];

    deltaX = mH - closestMass;

    deltaY = tmpHigh - tmpLow;
    slope = deltaY/step;

    // For partial widths                                             
    if(deltaX == 0){ reqCSerr = tmpLow;}
    else{ reqCSerr = slope*deltaX + tmpLow;}
    reqCSerr *= .01; //Account for percentage                                                                                          
  }

  return reqCSerr;

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)                  
double HiggsCSandWidthSM4::HiggsCSpdfErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  int i = 0;
  double closestMass = 0;
  double reqCSerr = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                                                                           
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}

  // If Ecm is not 7 TeV return -1                                                                                         
  if(sqrts != 7){return -1;}
 


  // If mH is out of range return -1                                                                                  
  // else find what array number to read                                                              
  if( mH < 100 || mH > 600){return 0;}
  else{

    if(mH <= 110 ){step = 5; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(2 + (mH - 110)/step); closestMass = (step*(i-2) + 110);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(62 + (mH - 140)/step); closestMass = (int)(step*(i-62) + 140);}
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(82 + (mH - 160)/step); closestMass = (int)(step*(i-82) + 160);}
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(147 + (mH - 290)/step); closestMass = (int)(step*(i-147) + 290);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(159 + (mH-350)/step); closestMass = (int)(step*(i-159) + 350);}
    if(mH > 400){step = 20; i = (int)(164 + (mH-400)/step); closestMass = (int)(step*(i-164) + 400);}




    tmpLow = CSpdfErrPlus[ID][i];
    tmpHigh = CSpdfErrPlus[ID][i+1];

    deltaX = mH - closestMass;

    deltaY = tmpHigh - tmpLow;
    slope = deltaY/step;

    // For partial widths                    
    if(deltaX == 0){ reqCSerr = tmpLow;}
    else{ reqCSerr = slope*deltaX + tmpLow;}
    reqCSerr *= .01; //Account for percentage                                                                                             
  }

  return reqCSerr;

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV (numbers are for 7 TeV only in this version)         
double HiggsCSandWidthSM4::HiggsCSpdfErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  int i = 0;
  double closestMass = 0;
  double reqCSerr = 0;
  double tmpLow = 0, tmpHigh = 0;
  double deltaX = 0, deltaY = 0;
  double slope = 0;
  double step = 0;


  // If ID is unavailable return -1                           
  if(ID > ID_ggToH || ID < ID_ggToH){return 0;}

  // If Ecm is not 7 TeV return -1                                                                 
  if(sqrts != 7){return -1;}
  
  // If mH is out of range return -1                                                              
  // else find what array number to read                            
  if( mH < 100 || mH > 1000){return 0;}
  else{


    //Find index and closest higgs mass for which we have numbers
    if(mH <= 110 ){step = 5; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(2 + (mH - 110)/step); closestMass = (step*(i-2) + 110);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(62 + (mH - 140)/step); closestMass = (int)(step*(i-62) + 140);}
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(82 + (mH - 160)/step); closestMass = (int)(step*(i-82) + 160);}
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(147 + (mH - 290)/step); closestMass = (int)(step*(i-147) + 290);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(159 + (mH-350)/step); closestMass = (int)(step*(i-159) + 350);}
    if(mH > 400){step = 20; i = (int)(164 + (mH-400)/step); closestMass = (int)(step*(i-164) + 400);}




    tmpLow = CSpdfErrMinus[ID][i];
    tmpHigh = CSpdfErrMinus[ID][i+1];

    deltaX = mH - closestMass;


    deltaY = tmpHigh - tmpLow;
    slope = deltaY/step;

    // For partial widths                                                          
    if(deltaX == 0){ reqCSerr = tmpLow;}
    else{ reqCSerr = slope*deltaX + tmpLow;}
    reqCSerr *= .01; //Account for percentage                                                                                    
  }

  return reqCSerr;

}



// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidthSM4::HiggsWidth(int ID, double mH, bool spline){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                   H->4e/4mu = 12               */
  /*                    H->2e2mu = 13               */
  /*             H->4l(e/mu/tau) = 14               */
  /*                       H->4q = 15               */
  /*                     H->2l2q = 16               */
  /*                       H->4f = 17               */
  /**************************************************/



  double TotalWidth = 0;
  double PartialWidth = 0;
  double Width = 0;
  int i = 0;
  double closestMass = 0;
  double tmpLow1, tmpHigh1, deltaX, deltaY1, slope1;
  double deltaY2, tmpLow2, tmpHigh2, slope2, step;


  // If ID is unavailable return -1                                           
  if(ID > 17 || ID < 0){return 0;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 100 || mH > 1000){return 0;}
  else{

    if( ID == 8 )
      {
	if(mH <= 150){step = 10; i = (int)((mH-100)/step); closestMass = (int)(step*i + 100);}
	else{ return 0;}
      }
    else{
      //Find index and closest higgs mass for which we have numbers
      if(mH <= 140){step = 10; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
      if(mH > 140 && mH <= 150 ){step = 5; i = (int)(4 + (mH-140)/step); closestMass = (step*(i-4) + 140);}
      if(mH > 150 && mH <= 190 ){step = 1; i = (int)(6 + (mH-150)/step); closestMass = (int)(step*(i-6) + 150);}
      if(mH > 190 && mH <= 200 ){step = 5; i = (int)(46 + (mH-190)/step); closestMass = (int)(step*(i-46) + 190);}
      if(mH > 200 && mH <= 330 ){step = 10; i = (int)(48 + (mH-200)/step); closestMass = (int)(step*(i-48) + 200);}
      if(mH > 330 && mH <= 335 ){step = 5; i = (int)(61 + (mH-330)/step); closestMass = (int)(step*(i-61) + 330);}
      if(mH > 335 && mH <= 340 ){step = 1; i = (int)(62 + (mH-335)/step); closestMass = (int)(step*(i-62) + 335);}
      if(mH > 340 && mH <= 345 ){step = 5; i = (int)(67 + (mH-340)/step); closestMass = (int)(step*(i-67) + 340);}
      if(mH > 345 && mH <= 355 ){step = 1; i = (int)(68 + (mH-345)/step); closestMass = (int)(step*(i-68) + 345);}
      if(mH > 355 && mH <= 360 ){step = 5; i = (int)(78 + (mH-355)/step); closestMass = (int)(step*(i-78) + 355);}
      if(mH > 360 && mH <= 370 ){step = 10; i = (int)(79 + (mH-360)/step); closestMass = (int)(step*(i-79) + 360);}
      if(mH > 370 && mH <= 380 ){step = 5; i = (int)(80 + (mH-370)/step); closestMass = (int)(step*(i-80) + 370);}
      if(mH > 380 && mH <= 500 ){step = 10; i = (int)(82 + (mH-380)/step); closestMass = (int)(step*(i-82) + 380);}
      if(mH > 500 && mH <= 800 ){step = 50; i = (int)(94 + (mH-500)/step); closestMass = (int)(step*(i-94) + 500);}
      if(mH > 800 && mH <= 1000 ){step = 200; i = (int)(100 + (mH-800)/step); closestMass = (int)(step*(i-100) + 800);}
    }


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
	      int indexWSM4 = 4;
	      double xmhWSM4[indexWSM4], sigWSM4[indexWSM4];
	      xmhWSM4[0]=mass_BR[i-1];xmhWSM4[1]=mass_BR[i];xmhWSM4[2]=mass_BR[i+1];xmhWSM4[3]=mass_BR[i+2];
	      sigWSM4[0]=BR[ID][i-1]; sigWSM4[1]=BR[ID][i]; sigWSM4[2]=BR[ID][i+1]; sigWSM4[3]=BR[ID][i+2];
	      
	      TGraph *graphWSM4 = new TGraph(indexWSM4, xmhWSM4, sigWSM4);
	      TSpline3 *gsWSM4 = new TSpline3("gsWSM4",graphWSM4);
	      gsWSM4->Draw();
	      Width = gsWSM4->Eval(mH);
	      delete gsWSM4;
	      delete graphWSM4;
	    }
	  else{
	    if(i < 1){i = 1;}
	    if(i+2 >= N_BR){i = N_BR - 3;}
	    if(ID == 8 && i+2 >= N_BR){i = N_BR_gg - 3;}

	    int indexWSM4 = 4;
	    double xmhWSM4[indexWSM4], sigWSM4[indexWSM4];
	    xmhWSM4[0]=mass_BR[i-1];xmhWSM4[1]=mass_BR[i];xmhWSM4[2]=mass_BR[i+1];xmhWSM4[3]=mass_BR[i+2];
	    sigWSM4[0]=BR[0][i-1]; sigWSM4[1]=BR[0][i]; sigWSM4[2]=BR[0][i+1]; sigWSM4[3]=BR[0][i+2];
	    
	    TGraph *graphWSM4 = new TGraph(indexWSM4, xmhWSM4, sigWSM4);
	    TSpline3 *gsWSM4 = new TSpline3("gsWSM4",graphWSM4);
	    gsWSM4->Draw();
	    PartialWidth = gsWSM4->Eval(mH);
	    delete gsWSM4;
	    delete graphWSM4;
	    
	    int indexPWSM4 = 4;
	    double xmhPWSM4[indexPWSM4], sigPWSM4[indexPWSM4];
	    xmhPWSM4[0]=mass_BR[i-1];xmhPWSM4[1]=mass_BR[i];xmhPWSM4[2]=mass_BR[i+1];xmhPWSM4[3]=mass_BR[i+2];
	    sigPWSM4[0]=BR[ID][i-1]; sigPWSM4[1]=BR[ID][i]; sigPWSM4[2]=BR[ID][i+1]; sigPWSM4[3]=BR[ID][i+2];
	    
	    TGraph *graphPWSM4 = new TGraph(indexPWSM4, xmhPWSM4, sigPWSM4);
	    TSpline3 *gsPWSM4 = new TSpline3("gsPWSM4",graphPWSM4);
	    gsPWSM4->Draw();
	    PartialWidth *= gsPWSM4->Eval(mH);
	    delete gsPWSM4;
	    delete graphPWSM4;
	    
	    Width = PartialWidth;
	    
	  }
	}
  
  }
  
  return Width;
  
  
}


double HiggsCSandWidthSM4::HiggsBR(int ID, double mH, bool spline){


  /***********************IDs************************/
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                   H->4e/4mu = 12               */
  /*                    H->2e2mu = 13               */
  /*             H->4l(e/mu/tau) = 14               */
  /*                       H->4q = 15               */
  /*                     H->2l2q = 16               */
  /*                       H->4f = 17               */
  /**************************************************/



  double PartialBR = 0;
  double BranchRatio = 0;
  int i = 0;
  double closestMass = 0;
  double tmpLow1, tmpHigh1, deltaX, deltaY1, slope1;
  double step;


  // If ID is unavailable return -1                                           
  if(ID > 17 || ID < 1){return 0;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 100 || mH > 1000){return 0;}
  else{

    if( ID == 8 )
      {
	if(mH <= 150){step = 10; i = (int)((mH-100)/step); closestMass = (int)(step*i + 100);}
	else{ return 0;}
      }
    else{
      //Find index and closest higgs mass for which we have numbers
      if(mH <= 140){step = 10; i = (int)((mH - 100)/step); closestMass = (int)(step*i + 100);}
      if(mH > 140 && mH <= 150 ){step = 5; i = (int)(4 + (mH-140)/step); closestMass = (step*(i-4) + 140);}
      if(mH > 150 && mH <= 190 ){step = 1; i = (int)(6 + (mH-150)/step); closestMass = (int)(step*(i-6) + 150);}
      if(mH > 190 && mH <= 200 ){step = 5; i = (int)(46 + (mH-190)/step); closestMass = (int)(step*(i-46) + 190);}
      if(mH > 200 && mH <= 330 ){step = 10; i = (int)(48 + (mH-200)/step); closestMass = (int)(step*(i-48) + 200);}
      if(mH > 330 && mH <= 335 ){step = 5; i = (int)(61 + (mH-330)/step); closestMass = (int)(step*(i-61) + 330);}
      if(mH > 335 && mH <= 340 ){step = 1; i = (int)(62 + (mH-335)/step); closestMass = (int)(step*(i-62) + 335);}
      if(mH > 340 && mH <= 345 ){step = 5; i = (int)(67 + (mH-340)/step); closestMass = (int)(step*(i-67) + 340);}
      if(mH > 345 && mH <= 355 ){step = 1; i = (int)(68 + (mH-345)/step); closestMass = (int)(step*(i-68) + 345);}
      if(mH > 355 && mH <= 360 ){step = 5; i = (int)(78 + (mH-355)/step); closestMass = (int)(step*(i-78) + 355);}
      if(mH > 360 && mH <= 370 ){step = 10; i = (int)(79 + (mH-360)/step); closestMass = (int)(step*(i-79) + 360);}
      if(mH > 370 && mH <= 380 ){step = 5; i = (int)(80 + (mH-370)/step); closestMass = (int)(step*(i-80) + 370);}
      if(mH > 380 && mH <= 500 ){step = 10; i = (int)(82 + (mH-380)/step); closestMass = (int)(step*(i-82) + 380);}
      if(mH > 500 && mH <= 800 ){step = 50; i = (int)(94 + (mH-500)/step); closestMass = (int)(step*(i-94) + 500);}
      if(mH > 800 && mH <= 1000 ){step = 200; i = (int)(100 + (mH-800)/step); closestMass = (int)(step*(i-100) + 800);}
    }


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
          if(i+2 >= N_BR){i = N_BR - 3;}
          if(ID == 8 && i+2 >= N_BR){i = N_BR_gg - 3;}

	  int indexBRSM4 = 4;
	  double xmhBRSM4[indexBRSM4], sigBRSM4[indexBRSM4];
	  xmhBRSM4[0]=mass_BR[i-1];xmhBRSM4[1]=mass_BR[i];xmhBRSM4[2]=mass_BR[i+1];xmhBRSM4[3]=mass_BR[i+2];
	  sigBRSM4[0]=BR[ID][i-1]; sigBRSM4[1]=BR[ID][i]; sigBRSM4[2]=BR[ID][i+1]; sigBRSM4[3]=BR[ID][i+2];
	  
	  TGraph *graphBRSM4 = new TGraph(indexBRSM4, xmhBRSM4, sigBRSM4);
	  TSpline3 *gsBRSM4 = new TSpline3("gsBRSM4",graphBRSM4);
	  gsBRSM4->Draw();
	  PartialBR = gsBRSM4->Eval(mH);
	  delete gsBRSM4;
	  delete graphBRSM4;
	}

      BranchRatio = PartialBR;
      
  }
  
  
  return BranchRatio;
  
} 


double HiggsCSandWidthSM4::HiggsBRErr_Hff(int ID, double mH, double sqrts)
{

  /***********************IDs************************/
  /*                   H->tautau = 2                */
  /*                   H->gamgam = 8                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /**************************************************/


  double sigma = 0.1;
  double err = 0;
  double BR_ = HiggsBR(1,mH,true)+HiggsBR(2,mH,true)+HiggsBR(3,mH,true)+HiggsBR(4,mH,true)+HiggsBR(5,mH,true)+HiggsBR(6,mH,true);

  // If ID is unavailable return -1                                           
  if(ID != 2 && ID != 8 && ID != 10 && ID != 11){return 0;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 100 || mH > 1000){return 0;}
  else{

    if( ID == 2 )
      {
	err = kappaFunc1(sigma,BR_);
      }
    else 
      {
	err = kappaFunc2(sigma,BR_);
      }
  
  }

  return err;

}

  

double HiggsCSandWidthSM4::HiggsBRErr_HVV(int ID, double mH, double sqrts)
{

  /***********************IDs************************/
  /*                   H->tautau = 2                */
  /*                   H->gamgam = 8                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /**************************************************/


  double sigma = 0.5;
  double err = 0;
  double BR_ = HiggsBR(10,mH,true)+HiggsBR(11,mH,true);

  // If ID is unavailable return -1                                           
  if(ID != 2 && ID != 8 && ID != 10 && ID != 11){return 0;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 100 || mH > 1000){return 0;}
  else{

    if( ID == 10 || ID == 11)
      {
	err = kappaFunc1(sigma,BR_);
      }
    else 
      {
	err = kappaFunc2(sigma,BR_);
      }
  
  }

  return err;

}



double HiggsCSandWidthSM4::HiggsBRErr_Hgluglu(int ID, double mH, double sqrts)
{

  /***********************IDs************************/
  /*                   H->tautau = 2                */
  /*                   H->gamgam = 8                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /**************************************************/


  double sigma = 0.05;
  double err = 0;
  double BR_ = HiggsBR(7,mH,true);

  // If ID is unavailable return -1                                           
  if(ID != 2 && ID != 8 && ID != 10 && ID != 11){return 0;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 100 || mH > 1000){return 0;}
  else{
	err = kappaFunc2(sigma,BR_);
  }

  return err;

}


double HiggsCSandWidthSM4::HiggsBRErr_Hgamgam(int ID, double mH, double sqrts)
{

  /***********************IDs************************/
  /*                   H->gamgam = 8                */
  /**************************************************/


  double sigma = 0.15;
  double err = 0;
  double BR_ = HiggsBR(8,mH,true);

  // If ID is unavailable return -1                                           
  if(ID != 8){return 0;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 100 || mH > 1000){return 0;}
  else{
	err = kappaFunc1(sigma,BR_);
  }

  return err;

}



double HiggsCSandWidthSM4::kappaFunc1(double sigma, double BR_)
{

  double kappa = 1 + sigma*(1-BR_);

  return kappa;

}


double HiggsCSandWidthSM4::kappaFunc2(double sigma, double BR_)
{

  double kappa = 1/(1+sigma*BR_);

  return kappa;

}



#endif
