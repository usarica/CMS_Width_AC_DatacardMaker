#ifndef HIGGSCSANDWIDTH_CC
#define HIGGSCSANDWIDTH_CC


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


#include "HiggsCSandWidth.h"

using namespace std;

HiggsCSandWidth::HiggsCSandWidth(std::string fileLoc)
{

  N_BR = 217;

  N_CS_7tev[ID_ggToH] = 197;
  N_CS_7tev[ID_VBF] = 197;
  N_CS_7tev[ID_WH]  = 152;
  N_CS_7tev[ID_ZH]  = 152;
  N_CS_7tev[ID_ttH] = 152;

  N_CS_8tev[ID_ggToH] = 223;
  N_CS_8tev[ID_VBF] = 223;
  N_CS_8tev[ID_WH]  = 178;
  N_CS_8tev[ID_ZH]  = 178;
  N_CS_8tev[ID_ttH] = 178;

  N_CS_14tev[ID_ggToH] = 50;
  N_CS_14tev[ID_VBF] = 50;
  N_CS_14tev[ID_WH]  = 33;
  N_CS_14tev[ID_ZH]  = 33;
  N_CS_14tev[ID_ttH] = 33;

  ifstream file;
  // ---------------- Read BR into memory ------------------ //         
  fileName = fileLoc+"/HiggsBR_7TeV_Official.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_BR; k++){

    file >> mass_BR[k] >> BR[0][k] >> BR[1][k] >> BR[2][k] >> BR[3][k] >> BR[4][k] >> BR[5][k] >> BR[6][k] >> BR[7][k] >> BR[8][k] >> BR[9][k]
	 >> BR[10][k] >> BR[11][k] >> BR[12][k] >> BR[13][k] >> BR[14][k] >> BR[15][k] >> BR[16][k] >> BR[17][k] >> BR[18][k] >> BR[19][k] >> BR[20][k]
	 >> BR[21][k] >> BR[22][k] >> BR[23][k] >> BR[24][k] >> BR[25][k];


  }
  file.close();

  // ---------------- Read 8 TeV CS into memory ------------------ //         
  fileName = fileLoc+"/7TeV-ggH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_7tev[ID_ggToH]; k++){

    file >> mass_XS_7tev[k][ID_ggToH] >> CS_7tev[k][ID_ggToH] >> CSerrPlus_7tev[k][ID_ggToH] >> CSerrMinus_7tev[k][ID_ggToH] 
	 >> CSscaleErrPlus_7tev[k][ID_ggToH] >> CSscaleErrMinus_7tev[k][ID_ggToH] >> CSpdfErrPlus_7tev[k][ID_ggToH] >> CSpdfErrMinus_7tev[k][ID_ggToH];
  
  }
  file.close();

  fileName = fileLoc+"/7TeV-vbfH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_7tev[ID_VBF]; k++){

    file >> mass_XS_7tev[k][ID_VBF] >> CS_7tev[k][ID_VBF] >> CSerrPlus_7tev[k][ID_VBF] >> CSerrMinus_7tev[k][ID_VBF] >> CSscaleErrPlus_7tev[k][ID_VBF]
	 >> CSscaleErrMinus_7tev[k][ID_VBF] >> CSpdfErrPlus_7tev[k][ID_VBF] >> CSpdfErrMinus_7tev[k][ID_VBF];

  }
  file.close();

  fileName = fileLoc+"/7TeV-ttH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_7tev[ID_ttH]; k++){

    file >> mass_XS_7tev[k][ID_ttH] >> CS_7tev[k][ID_ttH] >> CSerrPlus_7tev[k][ID_ttH] >> CSerrMinus_7tev[k][ID_ttH] >> CSscaleErrPlus_7tev[k][ID_ttH]
	 >> CSscaleErrMinus_7tev[k][ID_ttH] >> CSpdfErrPlus_7tev[k][ID_ttH] >> CSpdfErrMinus_7tev[k][ID_ttH];

  }
  file.close();

  fileName = fileLoc+"/7TeV-ZH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_7tev[ID_ZH]; k++){

    file >> mass_XS_7tev[k][ID_ZH] >> CS_7tev[k][ID_ZH] >> CSerrPlus_7tev[k][ID_ZH] >> CSerrMinus_7tev[k][ID_ZH] >> CSscaleErrPlus_7tev[k][ID_ZH]
	 >> CSscaleErrMinus_7tev[k][ID_ZH] >> CSpdfErrPlus_7tev[k][ID_ZH] >> CSpdfErrMinus_7tev[k][ID_ZH];
  }
  file.close();

  fileName = fileLoc+"/7TeV-WH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_7tev[ID_WH]; k++){

    file >> mass_XS_7tev[k][ID_WH] >> CS_7tev[k][ID_WH] >> CSerrPlus_7tev[k][ID_WH] >> CSerrMinus_7tev[k][ID_WH] >> CSscaleErrPlus_7tev[k][ID_WH]
	 >> CSscaleErrMinus_7tev[k][ID_WH] >> CSpdfErrPlus_7tev[k][ID_WH] >> CSpdfErrMinus_7tev[k][ID_WH];
  }
  file.close();

  // ---------------- Read 8 TeV CS into memory ------------------ //         
  fileName = fileLoc+"/8TeV-ggH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_8tev[ID_ggToH]; k++){

    file >> mass_XS_8tev[k][ID_ggToH] >> CS_8tev[k][ID_ggToH] >> CSerrPlus_8tev[k][ID_ggToH] >> CSerrMinus_8tev[k][ID_ggToH] 
	 >> CSscaleErrPlus_8tev[k][ID_ggToH] >> CSscaleErrMinus_8tev[k][ID_ggToH] >> CSpdfErrPlus_8tev[k][ID_ggToH] >> CSpdfErrMinus_8tev[k][ID_ggToH];
  
  }
  file.close();

  fileName = fileLoc+"/8TeV-vbfH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_8tev[ID_VBF]; k++){

    file >> mass_XS_8tev[k][ID_VBF] >> CS_8tev[k][ID_VBF] >> CSerrPlus_8tev[k][ID_VBF] >> CSerrMinus_8tev[k][ID_VBF] >> CSscaleErrPlus_8tev[k][ID_VBF]
	 >> CSscaleErrMinus_8tev[k][ID_VBF] >> CSpdfErrPlus_8tev[k][ID_VBF] >> CSpdfErrMinus_8tev[k][ID_VBF];

  }
  file.close();

  fileName = fileLoc+"/8TeV-ttH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_8tev[ID_ttH]; k++){

    file >> mass_XS_8tev[k][ID_ttH] >> CS_8tev[k][ID_ttH] >> CSerrPlus_8tev[k][ID_ttH] >> CSerrMinus_8tev[k][ID_ttH] >> CSscaleErrPlus_8tev[k][ID_ttH]
	 >> CSscaleErrMinus_8tev[k][ID_ttH] >> CSpdfErrPlus_8tev[k][ID_ttH] >> CSpdfErrMinus_8tev[k][ID_ttH];

  }
  file.close();

  fileName = fileLoc+"/8TeV-ZH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_8tev[ID_ZH]; k++){

    file >> mass_XS_8tev[k][ID_ZH] >> CS_8tev[k][ID_ZH] >> CSerrPlus_8tev[k][ID_ZH] >> CSerrMinus_8tev[k][ID_ZH] >> CSscaleErrPlus_8tev[k][ID_ZH]
	 >> CSscaleErrMinus_8tev[k][ID_ZH] >> CSpdfErrPlus_8tev[k][ID_ZH] >> CSpdfErrMinus_8tev[k][ID_ZH];
  }
  file.close();

  fileName = fileLoc+"/8TeV-WH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_8tev[ID_WH]; k++){

    file >> mass_XS_8tev[k][ID_WH] >> CS_8tev[k][ID_WH] >> CSerrPlus_8tev[k][ID_WH] >> CSerrMinus_8tev[k][ID_WH] >> CSscaleErrPlus_8tev[k][ID_WH]
	 >> CSscaleErrMinus_8tev[k][ID_WH] >> CSpdfErrPlus_8tev[k][ID_WH] >> CSpdfErrMinus_8tev[k][ID_WH];
  }
  file.close();


  // ---------------- Read 14 TeV CS into memory ------------------ //         
  fileName = fileLoc+"/14TeV-ggH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_14tev[ID_ggToH]; k++){

    file >> mass_XS_14tev[k][ID_ggToH] >> CS_14tev[k][ID_ggToH] >> CSerrPlus_14tev[k][ID_ggToH] >> CSerrMinus_14tev[k][ID_ggToH] 
	 >> CSscaleErrPlus_14tev[k][ID_ggToH] >> CSscaleErrMinus_14tev[k][ID_ggToH] >> CSpdfErrPlus_14tev[k][ID_ggToH] >> CSpdfErrMinus_14tev[k][ID_ggToH];
  
  }
  file.close();

  fileName = fileLoc+"/14TeV-vbfH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_14tev[ID_VBF]; k++){

    file >> mass_XS_14tev[k][ID_VBF] >> CS_14tev[k][ID_VBF] >> CSerrPlus_14tev[k][ID_VBF] >> CSerrMinus_14tev[k][ID_VBF] >> CSscaleErrPlus_14tev[k][ID_VBF]
	 >> CSscaleErrMinus_14tev[k][ID_VBF] >> CSpdfErrPlus_14tev[k][ID_VBF] >> CSpdfErrMinus_14tev[k][ID_VBF];

  }
  file.close();

  fileName = fileLoc+"/14TeV-ttH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_14tev[ID_ttH]; k++){

    file >> mass_XS_14tev[k][ID_ttH] >> CS_14tev[k][ID_ttH] >> CSerrPlus_14tev[k][ID_ttH] >> CSerrMinus_14tev[k][ID_ttH] >> CSscaleErrPlus_14tev[k][ID_ttH]
	 >> CSscaleErrMinus_14tev[k][ID_ttH] >> CSpdfErrPlus_14tev[k][ID_ttH] >> CSpdfErrMinus_14tev[k][ID_ttH];

  }
  file.close();

  fileName = fileLoc+"/14TeV-ZH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_14tev[ID_ZH]; k++){

    file >> mass_XS_14tev[k][ID_ZH] >> CS_14tev[k][ID_ZH] >> CSerrPlus_14tev[k][ID_ZH] >> CSerrMinus_14tev[k][ID_ZH] >> CSscaleErrPlus_14tev[k][ID_ZH]
	 >> CSscaleErrMinus_14tev[k][ID_ZH] >> CSpdfErrPlus_14tev[k][ID_ZH] >> CSpdfErrMinus_14tev[k][ID_ZH];
  }
  file.close();

  fileName = fileLoc+"/14TeV-WH.txt";
  file.open(fileName.c_str());
  for(int k = 0; k < N_CS_14tev[ID_WH]; k++){

    file >> mass_XS_14tev[k][ID_WH] >> CS_14tev[k][ID_WH] >> CSerrPlus_14tev[k][ID_WH] >> CSerrMinus_14tev[k][ID_WH] >> CSscaleErrPlus_14tev[k][ID_WH]
	 >> CSscaleErrMinus_14tev[k][ID_WH] >> CSpdfErrPlus_14tev[k][ID_WH] >> CSpdfErrMinus_14tev[k][ID_WH];
  }
  file.close();


}


HiggsCSandWidth::~HiggsCSandWidth()
{
  //destructor

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV 
double HiggsCSandWidth::HiggsCS(int ID, double mH, double sqrts){

  /**********IDs*************/ 
  /*     ggToH = 1          */
  /*       VBF = 2          */ 
  /*        WH = 3          */ 
  /*        ZH = 4          */
  /*       ttH = 5          */
  /*     Total = 0          */
  /**************************/
 
  double val = -1;

  // If ID is unavailable return -1                                                                                                
  if(ID > ID_ttH || ID < ID_Total) return -1;
  // If Ecm is not 7 or 8 TeV return -1
  if(sqrts != 7 && sqrts != 8 && sqrts != 14) return -1;
  //Don't interpolate btw 0 and numbers for mH300
  if(ID > ID_VBF && mH > 300) return 0;

  // If mH is out of range return -1                                           
  // else find what array number to read         
  if( mH < 90 || mH > 1000){ return -1;}
  else{
    
    if(sqrts == 7)
      {
	if(ID == 0)
	  {
	    for(int i = 1; i <= 5; i++)
	      {
		val += getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CS_7tev);
	      }
	  }
	else val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CS_7tev);
      }
    else if (sqrts == 8)
      {
	if(ID == 0)
	  {
	    for(int i = 1; i <= 5; i++)
	      {
		val += getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CS_8tev);
	      }
	  }
	else val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CS_8tev);
      }
    else if(sqrts == 14)
      {
        if(ID == 0)
          {
            for(int i =1; i <= 5; i++)
	      { 
		val += getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CS_14tev);
	      } 
          }
	else val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CS_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCS --- unknown sqrts! Choose 7,8, or 14." << endl; return -1;}
  }  
  
  return val;
}
  
  

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV 
double HiggsCSandWidth::HiggsCSErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;


  // If ID is unavailable return -1                                                                                    
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  if(ID == ID_Total){return 0;}
  // If Ecm is not 7 or 8 TeV return -1                                                                                                
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH300                        
  if(ID > ID_VBF && mH > 300){return 0;}

  // If mH is out of range return -1                                                                        
  // else find what array number to read                                          
  if( mH < 90 || mH > 1000){return -1;}
  else{

    if(sqrts == 7)
      {
	val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSerrPlus_7tev);
      }
    else if (sqrts == 8)
      {
	val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSerrPlus_8tev);
      }
    else if(sqrts == 14)
      {
	val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSerrPlus_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCSErrPlus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}
  }

  val *= .01;
  
  return val;
  
}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                                                                                       
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  if(ID == ID_Total){return 0;}
  // If Ecm is not 7 or 8 TeV return -1                                                                                           
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH300                                                        
  if(ID > ID_VBF && mH > 300){return 0;}


  // If mH is out of range return -1                                                                           
  // else find what array number to read                                                                 
  if( mH < 90 || mH > 1000){return -1;}
  else{

    if(sqrts == 7)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSerrMinus_7tev);
      }
    else if (sqrts == 8)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSerrMinus_8tev);
      }
    else if(sqrts == 14)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSerrMinus_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCSErrMinus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }

  val *= .01;

  return val;

}

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSscaleErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/


  double val = 0;

  // If ID is unavailable return -1                                                         
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  if(ID == ID_Total){return 0;}
  // If Ecm is not 7 or 8 TeV return -1                                                
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH300                                           
  if(ID > ID_VBF && mH > 300){return 0;}

  // If mH is out of range return -1                                                         
  // else find what array number to read                                                      
  if( mH < 90 || mH > 1000){return -1;}
  else{
    
    
    if(sqrts == 7)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSscaleErrPlus_7tev);
      }
    else if (sqrts == 8)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSscaleErrPlus_8tev);
      }
    else if(sqrts == 14)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSscaleErrPlus_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCSscaleErrPlus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }
  
  
  val *= .01; //Account for percentage  
  
  return val;
  
}

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSscaleErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                     
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  if(ID == ID_Total){return 0;}
  // If Ecm is not 7 or 8 TeV return -1                                                               
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH300                                   
  if(ID > ID_VBF && mH > 300){return 0;}


  // If mH is out of range return -1                        
  // else find what array number to read                              
  if( mH < 90 || mH > 1000){return -1;}
  else{

    if(sqrts == 7)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSscaleErrMinus_7tev);
      }
    else if (sqrts == 8)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSscaleErrMinus_8tev);
      }
    else if(sqrts == 14)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSscaleErrMinus_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCSscaleErrMinus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}
  }
  
  val *= .01; //Account for percentage  
  
  return val;

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSpdfErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                                                                           
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  if(ID == ID_Total){return 0;}
  // If Ecm is not 7 or 8 TeV return -1                                                                                         
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH300                                                  
  if(ID > ID_VBF && mH > 300){return 0;}

  // If mH is out of range return -1                                                                                  
  // else find what array number to read                                                              
  if( mH < 90 || mH > 1000){return -1;}
  else{

    if(sqrts == 7)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSpdfErrPlus_7tev);
      }
    else if (sqrts == 8)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSpdfErrPlus_8tev);
      }
    else if(sqrts == 14)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSpdfErrPlus_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCSpdfErrPlus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }
  
  val *= .01; //Account for percentage  
  
  return val;


}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSpdfErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                           
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  if(ID == ID_Total){return 0;}
  // If Ecm is not 7 or 8 TeV return -1                                                                 
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH300             
  if(ID > ID_VBF && mH > 300){return 0;}


  // If mH is out of range return -1                                                              
  // else find what array number to read                            
  if( mH < 90 || mH > 1000){return -1;}
  else{

    if(sqrts == 7)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSpdfErrMinus_7tev);
      }
    else if (sqrts == 8)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSpdfErrMinus_8tev);
      }
    else if(sqrts == 14)
      {
        val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSpdfErrMinus_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCSpdfErrMinus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }
  
  val *= .01; //Account for percentage  
  
  return val;

}



// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidth::HiggsWidth(int ID, double mH){


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
  /*                       H->4e = 12               */
  /*                    H->2e2mu = 13               */
  /*              H->4lep (e,mu) = 14               */
  /*          H->4lep (e,mu,tau) = 15               */
  /*                H->e+nu e-nu = 16               */
  /*               H->e+nu mu-nu = 17               */
  /*    H->2l2nu(l=e,mu)(nu=any) = 18               */
  /* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
  /*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
  /* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
  /* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
  /*  H->2nu2q (nu=any)(q=udcsb) = 23               */
  /*            H->4q (q=udcsb) = 24                */
  /*      H->4f (f=any fermion) = 25                */
  /**************************************************/



  double PartialWidth = 0;
  double Width = 0;
  int i = 0;
  double step;

  // If ID is unavailable return -1                                           
  if(ID > 25 || ID < 0){return -1;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 90 || mH > 1000){return -1;}
  else{

    //Find index and closest higgs mass for which we have numbers
    if(mH <=110 ){step = 5; i = (int)((mH - 90)/step);}
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(4 + (mH-110)/step);}
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(64 + (mH-140)/step); }
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(84 + (mH-160)/step); }
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(149 + (mH-290)/step); }
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(161 + (mH-350)/step); }
    if(mH > 400 && mH <= 600 ){step = 20; i = (int)(166 + (mH-400)/step); }
    if(mH > 600){step = 10; i = (int)(176 + (mH-600)/step); }


    if( ID == 0 )
      {
	if(i < 1){i = 1;}
	if(i+2 >= N_BR){i = N_BR - 3;}
	int indexW = 4;
	double xmhW[indexW], sigW[indexW];
	xmhW[0]=mass_BR[i-1];xmhW[1]=mass_BR[i];xmhW[2]=mass_BR[i+1];xmhW[3]=mass_BR[i+2];
	sigW[0]=BR[ID][i-1]; sigW[1]=BR[ID][i]; sigW[2]=BR[ID][i+1]; sigW[3]=BR[ID][i+2];
	
	TGraph *graphW = new TGraph(indexW, xmhW, sigW);
	TSpline3 *gsW = new TSpline3("gsW",graphW);
	gsW->Draw();
	Width = gsW->Eval(mH);
	delete gsW;
	delete graphW;
      }
    else{
      if(i < 1){i = 1;}
      if(i+2 >= N_BR){i = N_BR - 3;}
      
      int indexW = 4;
      double xmhW[indexW], sigW[indexW];
      xmhW[0]=mass_BR[i-1];xmhW[1]=mass_BR[i];xmhW[2]=mass_BR[i+1];xmhW[3]=mass_BR[i+2];
      sigW[0]=BR[0][i-1]; sigW[1]=BR[0][i]; sigW[2]=BR[0][i+1]; sigW[3]=BR[0][i+2];
      
      TGraph *graphW = new TGraph(indexW, xmhW, sigW);
      TSpline3 *gsW = new TSpline3("gsW",graphW);
      gsW->Draw();
      PartialWidth = gsW->Eval(mH);
      delete gsW;
      delete graphW;
      
      int indexPW = 4;
      double xmhPW[indexPW], sigPW[indexPW];
      xmhPW[0]=mass_BR[i-1];xmhPW[1]=mass_BR[i];xmhPW[2]=mass_BR[i+1];xmhPW[3]=mass_BR[i+2];
      sigPW[0]=BR[ID][i-1]; sigPW[1]=BR[ID][i]; sigPW[2]=BR[ID][i+1]; sigPW[3]=BR[ID][i+2];
      
      TGraph *graphPW = new TGraph(indexPW, xmhPW, sigPW);
      TSpline3 *gsPW = new TSpline3("gsPW",graphPW);
      gsPW->Draw();
      PartialWidth *= gsPW->Eval(mH);
      delete gsPW;
      delete graphPW;
      
      Width = PartialWidth;
      
    }
    
  }
  
  return Width;
  
} 


// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidth::HiggsBR(int ID, double mH){


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
  /*                       H->4e = 12               */
  /*                    H->2e2mu = 13               */
  /*              H->4lep (e,mu) = 14               */
  /*          H->4lep (e,mu,tau) = 15               */
  /*                H->e+nu e-nu = 16               */
  /*               H->e+nu mu-nu = 17               */
  /*    H->2l2nu(l=e,mu)(nu=any) = 18               */
  /* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
  /*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
  /* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
  /* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
  /*  H->2nu2q (nu=any)(q=udcsb) = 23               */
  /*            H->4q (q=udcsb) = 24                */
  /*      H->4f (f=any fermion) = 25                */
  /**************************************************/



  double PartialBR = 0;
  double BranchRatio = 0;
  int i = 0;
  double step;

  // If ID is unavailable return -1                                           
  if(ID > 25 || ID < 1){return -1;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < 90 || mH > 1000){return -1;}
  else{

    //Find index and closest higgs mass for which we have numbers
    if(mH <=110 ){step = 5; i = (int)((mH - 90)/step); }
    if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(4 + (mH-110)/step); }
    if(mH > 140 && mH <= 160 ){step = 1; i = (int)(64 + (mH-140)/step); }
    if(mH > 160 && mH <= 290 ){step = 2; i = (int)(84 + (mH-160)/step); }
    if(mH > 290 && mH <= 350 ){step = 5; i = (int)(149 + (mH-290)/step);}
    if(mH > 350 && mH <= 400 ){step = 10; i = (int)(161 + (mH-350)/step);}
    if(mH > 400 && mH <= 600 ){step = 20; i = (int)(166 + (mH-400)/step); }
    if(mH > 600){step = 10; i = (int)(176 + (mH-600)/step); }

    
    if(i < 1){i = 1;}
    if(i+2 >= N_BR){i = N_BR - 3;}
    int indexBR = 4;
    double xmhBR[indexBR], sigBR[indexBR];
    xmhBR[0]=mass_BR[i-1];xmhBR[1]=mass_BR[i];xmhBR[2]=mass_BR[i+1];xmhBR[3]=mass_BR[i+2];
    sigBR[0]=BR[ID][i-1]; sigBR[1]=BR[ID][i]; sigBR[2]=BR[ID][i+1]; sigBR[3]=BR[ID][i+2];
    
    TGraph *graphBR = new TGraph(indexBR, xmhBR, sigBR);
    TSpline3 *gsBR = new TSpline3("gsBR",graphBR);
    gsBR->Draw();
    PartialBR = gsBR->Eval(mH);
    delete gsBR;
    delete graphBR;
    
    BranchRatio = PartialBR;
    
  }

  return BranchRatio;

} 


double HiggsCSandWidth::getInterpXS(int sqrts, int ID, double mH, int maxI, double mhArray[][6], double varArray[][6])
{

  using namespace std;

  int i = 0;
  double reqCS = 0;
  double step = 0;
  int index = 4;
  double xmh[index], sig[index];
  
  if(sqrts == 7)
    {
      if(ID == ID_ggToH || ID == ID_VBF)
	{
	  if(mH <= 110 ){step = 5; i = (int)((mH - 90)/step); }
	  if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(4 + (mH - 110)/step); }
	  if(mH > 140 && mH <= 160 ){step = 1; i = (int)(64 + (mH - 140)/step); }
	  if(mH > 160 && mH <= 290 ){step = 2; i = (int)(84 + (mH - 160)/step); }
	  if(mH > 290 && mH <= 350 ){step = 5; i = (int)(149 + (mH - 290)/step); }
	  if(mH > 350 && mH <= 400 ){step = 10; i = (int)(161 + (mH-350)/step); }
	  if(mH > 400){step = 20; i = (int)(166 + (mH-400)/step); }
	}
      else if(ID == ID_WH || ID == ID_ZH || ID == ID_ttH)
	{
	  if(mH <= 110 ){step = 5; i = (int)((mH - 90)/step); }
          if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(4 + (mH - 110)/step);}
          if(mH > 140 && mH <= 160 ){step = 1; i = (int)(64 + (mH - 140)/step); }
          if(mH > 160 && mH <= 290 ){step = 2; i = (int)(84 + (mH - 160)/step); }
          if(mH > 290 && mH <= 300 ){step = 5; i = (int)(149 + (mH - 290)/step); }
          if(mH > 300) return 0;	  
	}
    }
  else if (sqrts == 8)
    {
      if(ID == ID_ggToH)
	{
	  if(mH <= 110 ){step = 1; i = (int)((mH - 80)/step); }
	  if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(30 + (mH - 110)/step); }
	  if(mH > 140 && mH <= 160 ){step = 1; i = (int)(90 + (mH - 140)/step); }
	  if(mH > 160 && mH <= 290 ){step = 2; i = (int)(110 + (mH - 160)/step); }
	  if(mH > 290 && mH <= 350 ){step = 5; i = (int)(175 + (mH - 290)/step); }
	  if(mH > 350 && mH <= 400 ){step = 10; i = (int)(187 + (mH - 350)/step); }
	  if(mH > 400 && mH <= 1000 ){step = 20; i = (int)(192 + (mH - 400)/step); }
	}
      else if(ID == ID_VBF)
	{
	  if(mH <= 110 ){step = 1; i = (int)((mH - 80)/step); }
	  if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(30 + (mH - 110)/step); }
	  if(mH > 140 && mH <= 160 ){step = 1; i = (int)(90 + (mH - 140)/step); }
	  if(mH > 160 && mH <= 290 ){step = 2; i = (int)(110 + (mH - 160)/step); }
	  if(mH > 290 && mH <= 350 ){step = 5; i = (int)(175 + (mH - 290)/step); }
	  if(mH > 350 && mH <= 400 ){step = 10; i = (int)(187 + (mH-350)/step); }
	  if(mH > 400){step = 20; i = (int)(192 + (mH-400)/step); }
	}
      else if(ID == ID_WH || ID == ID_ZH || ID == ID_ttH)
	{
	  if(mH <= 110 ){step = 1; i = (int)((mH - 80)/step); }
	  if(mH > 110 && mH <= 140 ){step = 0.5; i = (int)(30 + (mH - 110)/step);}
	  if(mH > 140 && mH <= 160 ){step = 1; i = (int)(90 + (mH - 140)/step); }
	  if(mH > 160 && mH <= 290 ){step = 2; i = (int)(110 + (mH - 160)/step); }
	  if(mH > 290 && mH <= 300 ){step = 5; i = (int)(175 + (mH - 290)/step); }
	  if(mH > 300) return 0;
	}
      else{ return 0;}
    }
  else if(sqrts == 14)
    {
      if(ID == ID_ggToH || ID == ID_VBF)
	{
	  if(mH <= 200 ){step = 5; i = (int)((mH - 90)/step); }
	  if(mH > 200 && mH <= 300 ){step = 10; i = (int)(22 + (mH - 200)/step); }
	  if(mH > 300 && mH <= 400 ){step = 20; i = (int)(32 + (mH - 300)/step); }
	  if(mH > 400 && mH <= 1000 ){step = 50; i = (int)(37 + (mH - 400)/step); }
	}
      else if(ID == ID_WH || ID == ID_ZH || ID == ID_ttH)
	{
	  if(mH <= 200 ){step = 5; i = (int)((mH - 90)/step); }
	  if(mH > 200 && mH <= 300 ){step = 10; i = (int)(22 + (mH - 200)/step); }
	}
      else{ return 0;}
    }
  else{cout << "HiggsCSandWidth --- unknown sqrts! Choose 7,8, or 14." << endl; return -1;}
 
  //Do the interpolation
  if(i < 1){i = 1;}
  if(i+2 >= maxI){i = maxI - 3;}
  xmh[0]=mhArray[i-1][ID];  xmh[1]=mhArray[i][ID];  xmh[2]=mhArray[i+1][ID];  xmh[3]=mhArray[i+2][ID];
  sig[0]=varArray[i-1][ID]; sig[1]=varArray[i][ID]; sig[2]=varArray[i+1][ID]; sig[3]=varArray[i+2][ID];
  
  TGraph *graph = new TGraph(index, xmh, sig);
  TSpline3 *gs = new TSpline3("gs",graph);
  gs->Draw();
  reqCS = gs->Eval(mH);
  delete gs;
  delete graph;

  return reqCS;  

}






#endif
