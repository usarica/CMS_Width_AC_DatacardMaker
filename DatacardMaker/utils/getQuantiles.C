void getQuantiles(TString str = "higgsCombineToys.root")
{

  float limitObs =32.98/4.15;
  TFile* f = new TFile(str,"OPEN");
  TTree* limit = f->Get("limit");
  double values;
  limit->SetBranchAddress("limit",&values);
  TH1F* h = new TH1F("Toys","Toys",100,0.,50.);
  const Int_t nq = 100;
  int nTrueE;
  for(int i = 0; i < limit->GetEntries(); i++)
    {
      limit->GetEntry(i);
      if(values > 1. && values <= 50)
	{
	  h->Fill(values);
	  nTrueE++;
	}
    }
  
  Double_t xq[nq];  // position where to compute the quantiles in [0,1]
  Double_t yq[nq];  // array to contain the quantiles
  for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
  h->GetQuantiles(nq,yq,xq);
  cout<<"entries: "<<nTrueE<<endl;  
  //show the original histogram in the top pad
  TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,700,900);
  c1->Divide(1,2);
  c1->cd(1);
  h->Draw();
  
  // show the quantiles in the bottom pad
  c1->cd(2);
  gPad->SetGrid();
  TGraph *gr = new TGraph(nq,xq,yq);
  gr->SetTitle("Quantiles");
  gr->SetMarkerStyle(21);
  gr->Draw("al");

  cout << "95: " << gr->Eval(0.025) << " " << gr->Eval(0.975) << endl;
  cout << "68: " << gr->Eval(0.16) << " " << gr->Eval(0.84) << endl;
  cout << "50: " << gr->Eval(0.5)  << endl;

  cout<<gr->Eval(0.025)<<","<<gr->Eval(0.16)<<","<<gr->Eval(0.5)<<","<<gr->Eval(0.84)<<","<<gr->Eval(0.975) << endl;
  //  double pval=0;
  for(double p =0; p<1;p=p+0.001){
    if(gr->Eval(p)>limitObs){
      cout<<"p-value "<<p<<endl;
      //pval=p;
      break;
    }
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *cToPlot = new TCanvas("toPlot","toPlot");
  cToPlot->cd();
  h->Draw();
  h->GetXaxis()->SetRange(3.0,70.0);
  cToPlot->SaveAs("toPlotComb.root");
  cToPlot->SaveAs("toPlotComb.png");
  cToPlot->SaveAs("toPlotComb.pdf");
  cToPlot->SaveAs("toPlotComb.eps");
}


