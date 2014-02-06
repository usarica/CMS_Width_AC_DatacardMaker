void plotScan1D(int mass = 240, int maxwidth = 30, bool blind = true){

  // gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);
  float gglimit = float(maxwidth);

  char boh[200];
  sprintf(boh,"higgsCombine1D_exp.MultiDimFit.mH%d.root", mass);
  TFile *f=new TFile(boh);
  TTree *t=(TTree*)f->Get("limit");
  t->Draw("2*deltaNLL:CMS_zz4l_GGsm", "deltaNLL > 0","PL");

  TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
  gr0->SetName("Exp1D");
  gr0->SetLineWidth(2);
  gr0->SetLineColor(kBlack);
  gr0->SetLineStyle(2);
  gr0->SetTitle("");

  if(!blind)
    {
      sprintf(boh,"higgsCombine1D_observed.MultiDimFit.mH%d.root", mass);
      TFile *f1=new TFile(boh);
      TTree *t1=(TTree*)f1->Get("limit");
      t1->Draw("2*deltaNLL:CMS_zz4l_GGsm", "deltaNLL > 0","PL");

      TGraph *gr1 = (TGraph*) gROOT->FindObject("Graph")->Clone();
      gr1->SetName("Obs1D");
      gr1->SetLineWidth(2);
      gr1->SetLineColor(kBlack);
      gr1->SetTitle("");
    }

  TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
  c1->cd();
  gr0->GetXaxis()->SetTitle("#Gamma/#Gamma_{SM}");
  gr0->GetYaxis()->SetTitle("-2 #Delta lnL");
  gr0->GetXaxis()->SetLabelSize(0.04);
  gr0->GetYaxis()->SetLabelSize(0.04);
  gr0->GetYaxis()->SetRangeUser(0.,12.);
  gr0->GetXaxis()->SetRangeUser(0.,gglimit);
  gr0->Draw("AL");
  if(!blind) gr1->Draw("L");

  TLegend *leg = new TLegend(0.25,0.73,0.5,0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  leg->AddEntry(gr0, "Expected","l");
  if(!blind) leg->AddEntry(gr1,"Observed","l");
  leg->Draw();

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
  oneSig->Draw();

  TPaveText *twoSig = new TPaveText(0.85,0.44,0.9,0.48,"NDC");
  twoSig->SetFillColor(0);
  twoSig->SetTextFont(42);
  twoSig->SetTextColor(kRed);
  twoSig->SetBorderSize(0);
  twoSig->AddText("2#sigma"); 
  twoSig->Draw();

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
  c1->SaveAs("can_scan1D_ggsm.gif");
}
