void compareScans(){
  gStyle->SetOptTitle(0);
  const int nfiles = 6;
  //2D fits
  //TString files[]={"cards_03_05_Unblind_093_2D/HCG/220/","cards_03_05_Unblind_2D/HCG/220/","cards_03_05_Unblind_093_2D/HCG/220/","cards_03_05_Unblind_2D/HCG/220_noSyst/","cards_03_05_Unblind_2D/HCG/220/"}
  //1D(m4l) fits
  //TString files[]={"cards_03_05_Unblind_093_1Dm4l/HCG/220/","cards_03_05_Unblind_1Dm4l/HCG/220/","cards_03_05_Unblind_093_1Dm4l/HCG/220/","cards_03_05_Unblind_1Dm4l/HCG/220_noSyst/","cards_03_05_Unblind_1Dm4l/HCG/220/"}
  //1D(Dgg)
  //TString files[]={"cards_03_09_bkgFix_093_1DDgg/HCG/220/","cards_03_09_bkgFix_1_1DDgg/HCG/220/","cards_03_09_bkgFix_093_1DDgg/HCG/220/","cards_03_09_bkgFix_1_1DDgg/HCG/220_noSyst/"};
  //Alternative K hyp.
  TString files[]={"cards_03_05_Unblind_093_2D/HCG/220/","cards_03_05_Unblind_093_2D/HCG/220_noSyst/","cards_093_AlternativeKBKG/HCG/220/","cards_093_AlternativeKBKG/HCG/220_noSyst/","cards_093_AlternativeKBKG_28/HCG/220/","cards_093_AlternativeKBKG_28/HCG/220_noSyst/"};
  //low r scans
  //TString files[]={"cards_lowRScan_mu1_2D/HCG/220_all/","cards_lowRScan_mu1_2D/HCG/220_2e2mu/","cards_lowRScan_mu1_2D/HCG/220_4e/","cards_lowRScan_mu1_2D/HCG/220_4mu/","cards_lowRScanObs_2D/HCG/220_all/"}

  //int colors[]={kBlack,kRed+1,kBlue,kGreen+2,kYellow+2,kYellow+2,kYellow+3,kBlack,kBlue,kRed+1,kGreen+2};
  int colors[]={kBlack,kBlack,kYellow+2,kYellow+2,kPink-1,kPink-1};

  //fit plots
  //TString grnames[]={"Observed","Expected #mu=#mu_{exp}","Expected #mu=#mu_{obs}","Expected #mu=1 w/o syst","Observed #mu=1"};
  //channels plots
  //TString grnames[]={"Observed #mu=1","Observed #mu=1, 2e2#mu","Observed #mu=1 4e","Observed #mu=1 4#mu","Observed #mu=0.93"};
  //Alternative K plots
  TString grnames[]={"Observed #mu=#mu_{obs}","Expected #mu=#mu_{obs}","Observed #mu=#mu_{obs} K=1/2.8","Expected #mu=#mu_{obs} K=1/2.8","Observed #mu=#mu_{obs} K=2.8","Expected #mu=#mu_{obs} K=2.8"};

  //tell this flag which are obsered
  bool obs[] = {1,0,0,0,1,0,0};
  int mass = 220;
  int maxwidth = 30.0;
  bool blind = true;
  bool uncBand = false;


  //values for 1DDgg_093, expected mu=0.93
  double limits95[]={6.1548, 7.98996, 12.193, 19.1693, 26.382};
  double limits68[]={2.3321, 3.17787, 5.65613, 10.973, 16.4882x};
  //values for 2D_093, expected mu=0.93
  //double limits95[]={ 5.48667, 7.1966,10.8518,17.2703,24.0398};
  //double limits68[]={2.07673,2.81889,4.97083,9.94989, 14.86};

  // gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);
  float gglimit = float(maxwidth);

  TMultiGraph *mg = new TMultiGraph();

  TGraph *g[nfiles];

  TLegend *leg = new TLegend(0.25,0.73,0.5,0.93);
  //TLegend *leg = c1->BuildLegend();
  leg->SetX1(0.22);
  leg->SetX2(0.5);
  leg->SetY1(0.7);
  leg->SetY2(0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  for(int i=0;i<nfiles;i++){
    char boh[200];
    //printf("%d\n",i);
    //if(i==3)mass=240;
    //TString filepath;filepath.Form("HCG/%d/",mass);
    TString obsString = "exp";
    if(obs[i])obsString="obs";
    int nDi =2;
    //printf("%d\n",i);
    sprintf(boh,"%shiggsCombine%dD_%s.MultiDimFit.mH%d.root", files[i].Data(),nDi,obsString.Data(),mass);
    TFile *f1=TFile::Open(boh);
    TTree *t1=(TTree*)f1->Get("limit");
    t1->Draw("2*deltaNLL:CMS_zz4l_GGsm", "deltaNLL > 0","PL");
    TGraph *gr0 = (TGraph*)gROOT->FindObject("Graph")->Clone();
    gr0->SetName(grnames[i].Data());
    gr0->SetLineWidth(2.5);
    gr0->SetLineColor(colors[i]);
    if(!obs[i])gr0->SetLineStyle(2);
    gr0->SetTitle(grnames[i].Data());
    leg->AddEntry(gr0);
    g[i]=(TGraph*)gr0->Clone();
    mg->Add(gr0,"l");
    double *y = gr0->GetY();
    double *x = gr0->GetX();
    int ipol=-1;
    for(int ipo=0;ipo<gr0->GetN();ipo++){
      if(y[ipo]<3.84&&y[ipo+1]>3.84){
	ipol = ipo;
	break;
      }
    }
    double a =  (y[ipol+1]-y[ipol])/(x[ipol+1]-x[ipol]);
    double b = y[ipol]-a*x[ipol];
    printf("%d) limit %.5f\n",i,(3.84-b)/a);
  }


  //  leg->AddEntry(gr0, "Expected - no Syst","l");
  //  leg->AddEntry(gr1, "Expected","l");
  //  leg->Draw();

  float lumi7TeV=5.1;
  float lumi8TeV=19.7;


  TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
  c1->cd();
  g[0]->Draw("AL");
  g[0]->GetXaxis()->SetTitle("#Gamma/#Gamma_{SM}");
  g[0]->GetYaxis()->SetTitle("-2 #Delta lnL");
  g[0]->GetXaxis()->SetLabelSize(0.04);
  g[0]->GetYaxis()->SetLabelSize(0.04);
  g[0]->GetYaxis()->SetRangeUser(0.,12.);//12
  g[0]->GetXaxis()->SetRangeUser(0.,gglimit);

  if(uncBand){
    TLine *l2_95=new TLine();
    //l2_95->SetLineStyle(9);
    l2_95->SetLineWidth(10);
    l2_95->SetLineColor(kYellow);
    l2_95->DrawLine(limits95[0],3.84,TMath::Min(limits95[4],gglimit),3.84);
    if(uncBand)l2_95->Draw();
    
    TLine *l2_68=new TLine();
    //l2_68->SetLineStyle(9);
    l2_68->SetLineWidth(10);
    l2_68->SetLineColor(kGreen);
    l2_68->DrawLine(limits95[1],3.84,limits95[3],3.84);
    if(uncBand)l2_68->Draw("same");
    
    TLine *l1_95=new TLine();
    //l1_95->SetLineStyle(9);
    l1_95->SetLineWidth(10);
    l1_95->SetLineColor(kYellow);
    l1_95->DrawLine(limits68[0],1,TMath::Min(limits68[4],gglimit),1);
    if(uncBand)l1_95->Draw("same");
    
    TLine *l1_68=new TLine();
    //l1_68->SetLineStyle(9);
    l1_68->SetLineWidth(10);
    l1_68->SetLineColor(kGreen);
    l1_68->DrawLine(limits68[1],1,limits68[3],1);
    if(uncBand)l1_68->Draw("same");
  }

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  TText *text = pt->AddText(0.01,0.5,"CMS");
  // text = pt->AddText(0.2,0.6,Form("#sqrt{s} = 7 TeV, L = %.1f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi7TeV,lumi8TeV));
  text = pt->AddText(0.2,0.6,Form("#sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi8TeV));
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
  l1->SetLineWidth(2);
  l1->SetLineColor(kRed);
  l1->DrawLine(0.0,1.0,gglimit,1.0);
  l1->Draw("same");
  TLine *l2=new TLine();
  l2->SetLineStyle(9);
  l2->SetLineWidth(2);
  l2->SetLineColor(kRed);
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
    //medians->SetMarkerSize(2);
    medians->SetMarkerColor(kBlue);
    medians->Draw("PSAME");  
  }

  c1->SaveAs("compareK_2.gif");
  c1->SaveAs("compareK_2.eps");
  c1->SaveAs("compareK_2.pdf");
  c1->SaveAs("compareK_2.png");

}
