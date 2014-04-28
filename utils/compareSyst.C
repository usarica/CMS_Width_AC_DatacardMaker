void compareSyst(){
  gStyle->SetOptTitle(0);
  const int nfiles = 16;
  TString files[]={"lumi_7TeV","lumi_8TeV","pdf_gg","CMS_eff_m","CMS_eff_e","CMS_hzz2e2mu_Zjets","CMS_hzz4mu_Zjets","CMS_hzz4e_Zjets","pdf_qqbar","QCDscale_qqH","CMS_widthH_kbkg","QCDscale_ggH","CMS_zz4l_VBFscale_syst","CMS_zz4l_ZXshape_syst","QCDscale_VV","EWKcorr_VV"};
  //int colors[]={kBlack,kRed+1,kBlue,kGreen+2,kYellow+2,kYellow+2,kYellow+3,kBlack,kBlue,kRed+1,kGreen+2};
  //int colors[]={kBlack,kTeal-5,kRed,kPink+5,kRed-7,kYellow-1};
  //int colors[]={kRed-7,kRed,kRed+2,kRed+4};
  int colors[] = {kPink+5,kRed,kBlue,kTeal-5,kBlack,kGreen+2,kYellow+2,kPink-1,kPink-4,kYellow+3,kRed-7,kRed+2,kRed+3,kOrange,kCyan,kMagenta};

  //fit plots
  TString grnames[]={"lumi_7TeV","lumi_8TeV","pdf_gg","CMS_eff_m","CMS_eff_e","CMS_hzz2e2mu_Zjets","CMS_hzz4mu_Zjets","CMS_hzz4e_Zjets","pdf_qqbar","QCDscale_qqH","CMS_widthH_kbkg","QCDscale_ggH","CMS_zz4l_VBFscale_syst","CMS_zz4l_ZXshape_syst","QCDscale_VV","EWKcorr_VV"};


  //TString plotLabel = "H#rightarrow ZZ#rightarrow 4l+2l2#nu";
  TString plotLabel = "H#rightarrow ZZ#rightarrow 4l";
  //TString plotLabel = "H#rightarrow ZZ#rightarrow 4l_{low}+4l_{high}+2l2#nu_{high}";

  TString toyPlotname = "toyParallel2D_093_95/toPlot.root";
  //TFile toyPlotFile = "toPlot";
  double limits95[]={ 5.01781,6.99495,10.6292,17.0564,24.4662};
  double limits68[]={ 1.65029,2.7713,4.90192,9.99193,15.0928 };
  //tell this flag which are obsered
  int maxwidth = 30.0;
  bool mev=true;
  bool printpval=false;
  bool uncBand =false;
  bool toyPlot =false;
  if(uncBand)toyPlot=false;
  TString outString = "testSyst";//"03_17_2DchanExp_093";

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);
  float gglimit = float(maxwidth);

  TGraph *g[nfiles];

  TLegend *leg = new TLegend(0.22,0.7,0.5,0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  TLegendEntry *tentry = leg->AddEntry((TObject*)0,plotLabel.Data(),"");
  tentry->SetTextSize(0.04);
  //leg->AddEntry((TObject*)0, "","");

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
    int nDi =2;
    //printf("%d\n",i);
    sprintf(boh,"higgsCombine2D_exp_%s.MultiDimFit.mH220.root", files[i].Data());
    TFile *f1=TFile::Open(boh);
    TTree *t1=(TTree*)f1->Get("limit");
    t1->Draw("2*deltaNLL:CMS_zz4l_GGsm", "deltaNLL > 0","PL");
    TGraph *gr0 = (TGraph*)gROOT->FindObject("Graph")->Clone();
    gr0->SetName(grnames[i].Data());
    gr0->SetLineWidth(2.5);
    gr0->SetLineColor(colors[i]);
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
    printf("%s limit@95CL %.2f\n",grnames[i].Data(),(3.84-b)/a);
    a =  (y[ipol+1]-y[ipol])/(x[ipol+1]-x[ipol])/fact;
    b = y[ipol]-a*x[ipol]*fact;
    printf("%s limit@68CL %.2f\n",grnames[i].Data(),(1-b)/a);
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
  if(mev)  g[0]->GetXaxis()->SetTitle("#Gamma [MeV]");
  g[0]->GetYaxis()->SetTitle("-2 #Delta lnL");
  g[0]->GetYaxis()->SetTitleSize(0.05);
  g[0]->GetXaxis()->SetTitleSize(0.05);
  g[0]->GetXaxis()->SetLabelSize(0.04);
  g[0]->GetYaxis()->SetLabelSize(0.04);
  float upLim =20.;
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
  TText *text = pt->AddText(0.01,0.5,"CMS ");
  // text = pt->AddText(0.2,0.6,Form("#sqrt{s} = 7 TeV, L = %.1f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi7TeV,lumi8TeV));
  text = pt->AddText(0.5,0.5,Form("#sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi8TeV));
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
    medians->SetMarkerSize(1.5);
    medians->SetMarkerColor(kBlue);
    medians->SetLineColor(0);
    //medians->Draw("PSAME");  
    //leg->AddEntry(medians,"Expected median","p");
  }
  if(toyPlot){
    TFile *fToy = TFile::Open(toyPlotname.Data());
    TCanvas *cToy = fToy->Get("toPlot");
    TH1F *htoy = (TH1F*)toPlot->FindObject("Toys");
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
