void compareScans(){

  const int nfiles = 3;
  //String files[]={"cards_test2D","cards_test2D_mysyst","cards_test2D_allsyst"};
  TString files[]={"cards_vbfsyst/HCG/220/","cards_vbfsyst/HCG/220/","cards_vbfsyst/HCG/220_noSyst/","cards_freezing2"};
  int colors[]={kBlack,kGreen+2,kBlue,kRed+1,kYellow+3};
  TString grnames[]={"Observed","Expected","Expected w/o syst"};

  bool obs[] = {1,0,0,0};
  int mass = 220;
  int maxwidth = 30;
  bool blind = true;

  // gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);
  float gglimit = float(maxwidth);

  TMultiGraph *mg = new TMultiGraph();

  for(int i=0;i<nfiles;i++){
    char boh[200];
    //if(i==3)mass=240;
    //TString filepath;filepath.Form("HCG/%d/",mass);
    TString obsString = "exp";
    if(obs[i])obsString="obs";
    sprintf(boh,"%shiggsCombine2D_%s.MultiDimFit.mH%d.root", files[i].Data(),obsString.Data(),mass);
    if(i==2)    sprintf(boh,"%shiggsCombine1D_%s.MultiDimFit.mH%d.root", files[i].Data(),obsString.Data(),mass);

    TFile *f1=TFile::Open(boh);
    TTree *t1=(TTree*)f1->Get("limit");
    t1->Draw("2*deltaNLL:CMS_zz4l_GGsm", "deltaNLL > 0","PL");
    TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
    gr0->SetName(grnames[i].Data());
    gr0->SetLineWidth(2.5);
    gr0->SetLineColor(colors[i]);
    gr0->SetLineStyle(2);
    gr0->SetTitle("");
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
    printf("%d) limit %.1f\n",i,(3.84-b)/a);
  }
  
  TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
  c1->cd();
  mg->Draw("AL");
  mg->GetXaxis()->SetTitle("#Gamma/#Gamma_{SM}");
  mg->GetYaxis()->SetTitle("-2 #Delta lnL");
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetRangeUser(0.,12.);
  mg->GetXaxis()->SetRangeUser(0.,gglimit);

  //TLegend *leg = new TLegend(0.25,0.73,0.5,0.93);
  TLegend *leg = c1->BuildLegend();
  leg->SetX1(0.22);
  leg->SetX2(0.5);
  leg->SetY1(0.7);
  leg->SetY2(0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  //  leg->AddEntry(gr0, "Expected - no Syst","l");
  //  leg->AddEntry(gr1, "Expected","l");
  //  leg->Draw();

  float lumi7TeV=5.1;
  float lumi8TeV=19.7;

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

  //c1->SaveAs("can_scan1D_ggsm.C");
  //c1->SaveAs("can_scan1D_ggsm.root");
  //c1->SaveAs("can_scan1D_ggsm.eps");
  c1->SaveAs("2DobsExp.gif");
  c1->SaveAs("2DobsExp.eps");
  c1->SaveAs("2DobsExp.png");

}
