
void compareScans(bool mev=true){
  gStyle->SetOptTitle(0);
  const int nfiles = 1;
  //TString files[]={"cards_07_04_14_update_2D_8TeV/HCG/220_2l2nu/","cards_07_04_14_update_2D_8TeV/HCG/220_2l2nu/","cards_07_04_14_update_2D_8TeV/HCG/220/","cards_07_04_14_update_2D_8TeV/HCG/220/","cards_07_04_14_update_2D_8TeV/HCG/only2l2nu/","cards_07_04_14_update_2D_8TeV/HCG/only2l2nu/",};
TString files[]={"cards_Paper/HCG/220/","cards_07_04_14_update_2D_8TeV/HCG/220_2l2nu/","cards_07_04_14_update_2D_8TeV/HCG/220/","cards_07_04_14_update_2D_8TeV/HCG/220/","cards_07_04_14_update_2D_8TeV/HCG/only2l2nu/","cards_07_04_14_update_2D_8TeV/HCG/only2l2nu/",};
  //TString files[]={"cards_03_17_Moriond_093_1DDgg/HCG/220/","cards_03_17_Moriond_093_1Dm4l/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_1DDgg/HCG/220/","cards_03_17_Moriond_093_1Dm4l/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_1_2D/HCG/220_noSyst/","cards_03_17_Moriond_1_2D/HCG/220/"}//Unblind
  //TString files[]={"cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220_postFit/"};
  //Combined 4l-2l2n
  //TString files[]={"cards_03_17_Combined/HCG/220/","cards_03_17_Combined/HCG/220/","cards_03_17_Combined/HCG/220_1/","cards_03_17_Combined/HCG/220_noSyst/"};
  //2D fits
  //TString files[]={"cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_1_2D/HCG/220/","cards_03_17_Moriond_1_2D/HCG/220_noSyst/"}
  //TString files[]={"cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_1_2D/HCG/220_noSyst/"}
  //TString files[]={"cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_1_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_1_2D/HCG/220_noSyst/","cards_03_17_Moriond_1_2D/HCG/220/"}//Unblind
  //TString files[]={"cards_03_05_Unblind_2D/HCG/220/","cards_03_05_Unblind_093_2D/HCG/220/","cards_03_05_Unblind_2D/HCG/220_noSyst/","cards_03_05_Unblind_2D/HCG/220/"}//Unblind
  //TString files[]={"cards_03_11_Approval_093_2D/HCG/220_all/","cards_03_11_Approval_093_2D/HCG/220_all/"};//Approval
  //2D channel splitting
  //TString files[]={"cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220_2e2mu/","cards_03_17_Moriond_093_2D/HCG/220_4e/","cards_03_17_Moriond_093_2D/HCG/220_4mu/","cards_03_17_Moriond_093_2D/HCG/220/","cards_03_17_Moriond_093_2D/HCG/220_2e2mu/","cards_03_17_Moriond_093_2D/HCG/220_4e/","cards_03_17_Moriond_093_2D/HCG/220_4mu/"};
  //1D(m4l) fits
  //TString files[]={"cards_03_17_Moriond_093_1Dm4l/HCG/220/","cards_03_17_Moriond_093_1Dm4l/HCG/220/","cards_03_17_Moriond_1_1Dm4l/HCG/220/","cards_03_17_Moriond_1_1Dm4l/HCG/220_noSyst/"}
  //1D(Dgg)
  //TString files[]={"cards_03_17_Moriond_093_1DDgg/HCG/220/","cards_03_17_Moriond_093_1DDgg/HCG/220/","cards_03_17_Moriond_1_1DDgg/HCG/220/","cards_03_17_Moriond_1_1DDgg/HCG/220_noSyst/"}
  //Alternative K hyp.
  //TString files[]={"cards_03_05_Unblind_093_2D/HCG/220/","cards_03_05_Unblind_093_2D/HCG/220_noSyst/","cards_093_AlternativeKBKG/HCG/220/","cards_093_AlternativeKBKG/HCG/220_noSyst/","cards_093_AlternativeKBKG_28/HCG/220/","cards_093_AlternativeKBKG_28/HCG/220/"};
  //TString files[]={"cards_KBKG_noUnc/HCG/220/","cards_KBKG_noUnc/HCG/220/","cards_AltKBKG_noUnc/HCG/220/","cards_AltKBKG_noUnc/HCG/220_noSyst/","cards_AltKBKG28_noUnc/HCG/220/","cards_AltKBKG28_noUnc/HCG/220_noSyst/"};
  //low r scans
  //TString files[]={"cards_03_17_Combined/HCG/220_lowR/","cards_03_17_Combined/HCG/220/","cards_03_17_Combined/HCG/220_1/","cards_03_17_Combined/HCG/220_noSyst/"};
  //TString files[]={"cards_lowRScan_mu1_2D/HCG/220_all/","cards_lowRScan_mu1_2D/HCG/220_2e2mu/","cards_lowRScan_mu1_2D/HCG/220_4e/","cards_lowRScan_mu1_2D/HCG/220_4mu/","cards_lowRScanObs_2D/HCG/220_all/"};
  //TString files[]={"cards_lowRScanObs_2D/HCG/220_all/","cards_lowRScanObs_2D/HCG/220_2e2mu/","cards_lowRScanObs_2D/HCG/220_4e/","cards_lowRScanObs_2D/HCG/220_4mu/"};
  //TString files[]={"cards_03_17_lowRscan2D_093/HCG/220_all/","cards_03_17_lowRscan2D_093/HCG/220_2e2mu/","cards_03_17_lowRscan2D_093/HCG/220_4e/","cards_03_17_lowRscan2D_093/HCG/220_4mu/"}
  //13TeV projections and lumi
  //TString files[]={"cards_03_11_Approval_093_2D/HCG/220_all/","cards_03_11_13TeV_2013lumi_093_2D/HCG/220/","cards_03_11_13TeV_100fblumi_093_2D/HCG/220/"};
  //TString files[]={"cards_03_11_Approval_093_2D/HCG/220_all/","cards_03_11_8TeV_100fblumi_093_2D/HCG/220/","cards_03_11_8TeV_300fblumi_093_2D/HCG/220/","cards_03_11_8TeV_3000fblumi_093_2D/HCG/220/"};

  //int colors[]={kBlack,kBlack,kRed+1,kBlue,kGreen+2,kBlack,kRed+1,kBlue,kGreen+2,kYellow+2,kYellow+2,kYellow+3,kBlack,kBlue,kRed+1,kGreen+2};
  //int colors[]={kBlack,kRed+1,kBlue,kGreen+2,kBlack,kRed+1,kBlue,kGreen+2,kYellow+2,kYellow+2,kYellow+3,kBlack,kBlue,kRed+1,kGreen+2};
  //int colors[]={kBlack,kRed+1,kBlue,kBlack,kRed+1,kBlue,kGreen+2,kYellow+2,kYellow+2,kYellow+3,kBlack,kBlue,kRed+1,kGreen+2};
  //int colors[]={kBlack,kBlack,kYellow+2,kYellow+2,kPink-1,kPink-1,kPink-4,kPink-4,kYellow-1,kYellow-1};
  //int colors[]={kBlack,kRed+1,kBlue,kGreen+2,kYellow+2,kYellow+2,kYellow+3,kBlack,kBlue,kRed+1,kGreen+2};
  //int colors[]={kBlack,kTeal-5,kRed,kPink+5,kRed-7,kYellow-1};
  //int colors[]={kRed-7,kRed,kRed+2,kRed+4};
  int colors[] = {kBlue,kBlue,kBlue-4,kBlue-4};

  //fit plots
  TString grnames[]={"Observed","Expected","Observed 4l","Expected 4l"};
  //TString grnames2[] = {"Observed","Expected"};
  //TString grnames[]={"Expected #mu=#mu_{obs}","Expected #mu=1","Expected #mu=1 w/o syst","Observed #mu=1"};
  //TString grnames[]={"Observed","Expected #mu=#mu_{obs}","Expected #mu=1 w/o syst","Observed #mu=1"};
  //TString grnames[]={"Observed","Expected #mu=#mu_{obs}","Expected #mu=#mu_{obs} w/o syst","Observed #mu=1"};
  //TString grnames[]={"Observed","Expected #mu=#mu_{obs}","Expected #mu=#mu_{obs} w/o syst","Observed #mu=1","Observed","Expected #mu=#mu_{obs}","Expected #mu=#mu_{obs} w/o syst","Observed #mu=1"};
  //TString grnames[]={"Observed","Expected","Observed allRaw","Expected allRaw","Observed allk3a","Expected allk3a","Observed smooth+k3a","Expected smooth+k3a"};
  //channels plots
  //TString grnames[]={"Observed #mu=#mu_{obs} 4l","Observed #mu=#mu_{obs} 2e2#mu","Observed #mu=#mu_{obs} 4e","Observed #mu=#mu_{obs} 4#mu","Expected #mu=#mu_{obs} 4l","Expected #mu=#mu_{obs} 2e2#mu","Expected #mu=#mu_{obs} 4e","Expected #mu=#mu_{obs} 4#mu"};
  //TString grnames[]={"Expected #mu=#mu_{obs} 4l","Expected #mu=#mu_{obs} 2e2#mu","Expected #mu=#mu_{obs} 4e","Expected #mu=#mu_{obs} 4#mu"};
  //Alternative K plots
  //TString grnames[]={"Observed #mu=#mu_{obs}","Expected #mu=#mu_{obs}","Observed #mu=#mu_{obs} K=1/2.8","Expected #mu=#mu_{obs} K=1/2.8","Observed #mu=#mu_{obs} K=2.8","Expected #mu=#mu_{obs} K=2.8","Observed #mu=#mu_{obs} K=2.8 w/o unc","Observed #mu=#mu_{obs} K=1/2.8 w/o unc"};
  //13 TeV
  //TString grnames[]={"Expected 2013","Expected 13TeV 17.912/fb","Expected 13TeV 100/fb","Expected #mu=1 w/o syst","Observed #mu=1"};
  //TString grnames[]={"Expected 2014, 19.712/fb","Expected 100/fb","Expected 300/fb","Expected 3000/fb","Observed #mu=1"};


  //TString plotLabel = "H#rightarrow ZZ#rightarrow 4l+2l#nu";
  TString plotLabel = "H #rightarrow ZZ";
  //TString plotLabel = "H#rightarrow ZZ#rightarrow 4l_{105.6-140.6GeV}+4l_{220+GeV}+2l2#nu_{220+GeV}";

  //TString toyPlotname = "toyParallel2D_093_95/toPlot.root";
  //TFile toyPlotFile = "toPlot";

  //tell this flag which are obsered
  bool obs[] = {1,0,1,0,1,0,1,0,1,0,1,0};
  double mass[] = {125.6,125.6,125.6,125.6,125.6,125.6,125.6,125.6,125.6,125.6,125.6,125.6};
  int maxwidth = 12.0;
  float upLim =11.;
  //bool mev=true;
  bool printpval=false;
  bool uncBand =false;
  bool toyPlot =false;
  if(uncBand)toyPlot=false;
  TString outString = "testCombined_11_04_14";//"03_17_2DchanExp_093";

  //values for 1DDgg_093, expected mu=0.93
  //double limits95[]={6.05994,7.97529,12.1601,18.8235,26.9906};
  //double limits68[]={2.40192,3.22084,5.71285,10.862,16.2785};
  //values for 1Dm4l_093, expected mu=0.93
  //double limits95[]={ 8.00674,10.8502,16.7568,26.2496,35.4528};
  //double limits68[]={3.09662,4.28213,7.78485,14.8901,22.1333 }:
  //values for 2D_093, expected mu=0.93
  double limits95[]={ 5.01781,6.99495,10.6292,17.0564,24.4662};
  double limits68[]={ 1.65029,2.7713,4.90192,9.99193,15.0928 };
  //values for 4l+2l2n Combination (from Chris)
  //double limits95[]={4.382,5.973,9.09,13.879,18.823};
  //double limits68[]={1.594,2.567,4.925,8.409,11.752};

  // gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);
  float gglimit = float(maxwidth);

  TGraph *g[nfiles];

  TLegend *leg = new TLegend(0.27,0.7,0.75,0.9);
  //TLegend *leg = new TLegend(0.18,0.55,0.5,0.93);
  //TLegend *leg = new TLegend(0.55,0.41,0.81,0.71);
  //leg->SetX1(0.22);
  //leg->SetX2(0.5);
  //leg->SetY1(0.7);
  //leg->SetY2(0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  TLegendEntry *tentry = leg->AddEntry((TObject*)0,plotLabel.Data(),"");
  tentry->SetTextSize(0.05);
  tentry->SetTextFont(62);
  leg->AddEntry((TObject*)0, "","");

  double limitObs =0;

  if(mev){
    for (int imev=0;imev<5;imev++){
      limits95[imev]=limits95[imev]*4.15;
      limits68[imev]=limits68[imev]*4.15;
    }
    maxwidth=maxwidth*4.15;
    outString+="_MeV";
    gglimit*=4.15;
  }

  for(int i=0;i<nfiles;i++){
    char boh[200];
    //printf("%d\n",i);
    //if(i==3)mass=240;
    //TString filepath;filepath.Form("HCG/%d/",mass);
    TString obsString = "exp";
    //if(i)obsString="NoSyst";
    if(obs[i])obsString="obs";
    int nDi =2;
    //printf("%d\n",i);
    sprintf(boh,"%shiggsCombine%dD_%s.MultiDimFit.mH%.1f.root", files[i].Data(),nDi,obsString.Data(),mass[i]);
    TFile *f1=TFile::Open(boh);
    //TFile *f1=TFile::Open("cards_Paper/HCG/220/higgsCombineObs_Combined_nLL_scan.MultiDimFit.mH125.6.root");
    TTree *t1=(TTree*)f1->Get("limit");
    t1->Draw("2*deltaNLL:CMS_zz4l_GGsm", "deltaNLL > 0","PL");
    TGraph *gr0 = (TGraph*)gROOT->FindObject("Graph")->Clone();
    gr0->SetName(grnames[i].Data());
    gr0->SetLineWidth(3);
    gr0->SetLineColor(colors[i]);
    if(i == 1)gr0->SetLineStyle(9);
    else if (i == 2)gr0->SetLineStyle(2);
    else if (i == 3)gr0->SetLineStyle(10);
    gr0->SetTitle(grnames[i].Data());
    leg->AddEntry(gr0);
    g[i]=(TGraph*)gr0->Clone();
    double *y = gr0->GetY();
    double *x = gr0->GetX();
    int ipol=-1,ipol68=-1;
    for(int ipo=0;ipo<gr0->GetN();ipo++){
      if(y[ipo]<3.84&&y[ipo+1]>3.84){
	ipol = ipo;
      }
      if(y[ipo]<1&&y[ipo+1]>1){
	ipol68 = ipo;
      }
      if(mev)g[i]->SetPoint(ipo,x[ipo]*4.15,y[ipo]);
    }
    double fact=1;
    if(mev)fact=4.15;
    double a =  (y[ipol+1]-y[ipol])/(x[ipol+1]-x[ipol])/fact;
    double b = y[ipol]-a*x[ipol]*fact;
    if(obs[i] && limitObs<1)limitObs=(3.84-b)/a;
    printf("%s limit@95CL %.2f\n",grnames[i].Data(),(3.84-b)/a);
    a =  (y[ipol68+1]-y[ipol68])/(x[ipol68+1]-x[ipol68])/fact;
    b = y[ipol68]-a*x[ipol68]*fact;
    printf("%s limit@68CL %.2f\n",grnames[i].Data(),(1-b)/a);
  }


  //  leg->AddEntry(gr0, "Expected - no Syst","l");
  //  leg->AddEntry(gr1, "Expected","l");
  //  leg->Draw();

  float lumi7TeV=5.1;
  float lumi8TeV=19.7;

  //TFile *fToy = TFile::Open(toyPlotname.Data());
  //TCanvas *cToy = fToy->Get("toPlot");
  //TH1F *htoy = (TH1F*)toPlot->FindObject("Toys");

  TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
  c1->cd();
  g[0]->Draw("AL");
  g[0]->GetXaxis()->SetTitle("#Gamma/#Gamma_{SM}");
  if(mev)  g[0]->GetXaxis()->SetTitle("#Gamma_{H} (MeV)");
  g[0]->GetXaxis()->SetTitleOffset(0.8);
  g[0]->GetYaxis()->SetTitle("-2 #Delta lnL");
  g[0]->GetYaxis()->SetTitleSize(0.05);
  g[0]->GetXaxis()->SetTitleSize(0.05);
  g[0]->GetXaxis()->SetLabelSize(0.04);
  g[0]->GetYaxis()->SetLabelSize(0.04);
  if(gglimit<5)upLim=1.01;
  g[0]->GetYaxis()->SetRangeUser(0.,upLim);//12
  g[0]->GetXaxis()->SetRangeUser(0.,gglimit);

  if(uncBand){
    TLine *l2_95=new TLine();
    //l2_95->SetLineStyle(9);
    l2_95->SetLineWidth(20);
    l2_95->SetLineColor(kYellow);
    l2_95->DrawLine(limits95[0],3.84,TMath::Min(limits95[4],gglimit),3.84);
    l2_95->Draw();
    
    TLine *l2_68=new TLine();
    //l2_68->SetLineStyle(9);
    l2_68->SetLineWidth(20);
    l2_68->SetLineColor(kGreen);
    l2_68->DrawLine(limits95[1],3.84,limits95[3],3.84);
    l2_68->Draw("same");
    
    TLine *l1_95=new TLine();
    //l1_95->SetLineStyle(9);
    l1_95->SetLineWidth(20);
    l1_95->SetLineColor(kYellow);
    l1_95->DrawLine(limits68[0],1,TMath::Min(limits68[4],gglimit),1);
    l1_95->Draw("same");
    
    TLine *l1_68=new TLine();
    //l1_68->SetLineStyle(9);
    l1_68->SetLineWidth(20);
    l1_68->SetLineColor(kGreen);
    l1_68->DrawLine(limits68[1],1,limits68[3],1);
    leg->AddEntry(l1_68,"68% CL","l");
    leg->AddEntry(l1_95,"95% CL","l");
    l1_68->Draw("same");

  }
//   //p-value 0.13;
  if(printpval){
    TPaveText *pval = new TPaveText(0.54,0.35,0.85,0.45,"brNDC");
    pval->SetBorderSize(0);
    pval->SetTextAlign(12);
    pval->SetFillStyle(0);
    pval->SetTextFont(42);
    pval->SetTextSize(0.03);
    pval->AddText(0,0,"p-value @ 95%CL = 0.12");
    pval->Draw();
  }

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  //TText *text = pt->AddText(0.01,0.5,"CMS Preliminary");
  TText *text = pt->AddText(0.01,0.5,"CMS");
  text->SetTextFont(62);
  text->SetTextSize(0.03146853);
  text = pt->AddText(0.15,0.6,Form("#sqrt{s} = 7 TeV, L = %.1f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi7TeV,lumi8TeV));
  //text = pt->AddText(0.5,0.5,Form("#sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi8TeV));
  text->SetTextFont(62);
  text->SetTextSize(0.03146853);
  pt->Draw();  

  TPaveText *oneSig = new TPaveText(0.85,0.18,0.9,0.22,"NDC");
  oneSig->SetFillColor(0);
  oneSig->SetTextFont(42);
  oneSig->SetTextColor(kRed);
  oneSig->SetBorderSize(0);
  oneSig->AddText("1#sigma"); 
  //  oneSig->Draw();

  TPaveText *twoSig = new TPaveText(0.85,0.44,0.9,0.48,"NDC");
  twoSig->SetFillColor(0);
  twoSig->SetTextFont(42);
  twoSig->SetTextColor(kRed);
  twoSig->SetBorderSize(0);
  twoSig->AddText("2#sigma"); 
  //  twoSig->Draw();

  TLine *l1=new TLine();
  l1->SetLineStyle(9);
  l1->SetLineWidth(1);
  l1->SetLineColor(kBlack);
  l1->DrawLine(0.0,1.0,gglimit,1.0);
  l1->Draw("same");
  TLine *l2=new TLine();
  l2->SetLineStyle(9);
  l2->SetLineWidth(1);
  l2->SetLineColor(kBlack);
  l2->DrawLine(0.0,3.84,gglimit,3.84);
  l2->Draw("same");

  for(int i=0;i<nfiles;i++)g[i]->Draw("LSAME");
  leg->Draw("SAME");

  if(uncBand){
    TGraph *medians = new TGraph(2);
    medians->SetPoint(0,limits95[2],3.84);
    medians->SetPoint(1,limits68[2],1.0);
    medians->SetFillStyle(0);
    medians->SetMarkerStyle(30);
    medians->SetMarkerSize(1.5);
    medians->SetMarkerColor(kBlue);
    medians->SetLineColor(0);
    //medians->Draw("PSAME");  
    //leg->AddEntry(medians,"Expected median","p");
  }
  if(toyPlot){
    TPad *p = new TPad("p","p",0.6,4.5/12.0,0.9,7.3/12.0);
    p->SetMargin(0.05,0,0.05,0);
    p->Draw();
    p->cd();
    htoy->SetStats(0);
    htoy->Draw();
    //TText *pvalt = new TText(limitObs+1,20,"p-value @ 95%CL = 0.12");
    TText *pvalt = new TText(12,200,"p-value @ 95%CL = 0.12");
    pvalt->SetTextSize(0.08);
    TArrow *ar2 = new TArrow(limitObs,0.0,limitObs,htoy->GetBinContent(htoy->FindBin(limitObs)),0.02,"<|");
    ar2->SetAngle(30);
    ar2->Draw();  
    pvalt->Draw();
  }

  TString saveString;
  saveString.Form("%s.C",outString.Data());
  c1->SaveAs(saveString.Data());
  //saveString.Form("%s.gif",outString.Data());
  //c1->SaveAs(saveString.Data());
  saveString.Form("%s.eps",outString.Data());
  c1->SaveAs(saveString.Data());
  saveString.Form("%s.pdf",outString.Data());
  c1->SaveAs(saveString.Data());
  saveString.Form("%s.png",outString.Data());
  c1->SaveAs(saveString.Data());
  saveString.Form("%s.root",outString.Data());
  c1->SaveAs(saveString.Data());

}
