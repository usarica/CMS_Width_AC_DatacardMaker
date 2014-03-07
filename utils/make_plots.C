make_plots()
{

  bool sumChannels = false;

  gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
  gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
  gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");

  float lumi7TeV=5.1;
  float lumi8TeV=19.7;

  TCanvas* c1 = new TCanvas("c1","c1",800,800);

  TFile* TWOeTWOmu_f = new TFile("hzz4l_2e2muS_8TeV.input.root","OPEN");
  RooWorkspace* TWOeTWOmu_w = TWOeTWOmu_f->Get("w");
  RooAbsPdf* Zjets_2e2mu = TWOeTWOmu_w->pdf("bkg_zjets");
  RooAbsPdf* qqZZ_2e2mu = TWOeTWOmu_w->pdf("bkg_qqzz");
  RooRealVar* TWOeTWOmu_mu = TWOeTWOmu_w->var("CMS_zz4l_mu");
  RooRealVar* TWOeTWOmu_GGsm = TWOeTWOmu_w->var("CMS_zz4l_GGsm");
  RooRealVar* TWOeTWOmu_mass = TWOeTWOmu_w->var("CMS_zz4l_widthMass");
  RooRealVar* TWOeTWOmu_kd = TWOeTWOmu_w->var("CMS_zz4l_widthKD");

  RooDataSet* data_2e2mu = TWOeTWOmu_w->data("data_obs")->Clone("data_2e2mu");
  RooDataSet* data = data_2e2mu->Clone("data_full");

  ifstream TWOeTWOmu_card;
  TWOeTWOmu_card.open("hzz4l_2e2muS_8TeV.txt");
  char line[256];
  for (int n_line=0; n_line < 7; n_line++)
    {
      TWOeTWOmu_card.getline(line,256);
    }
  TWOeTWOmu_card.getline(line,256);
  for(int n_line=0; n_line < 3; n_line++)
    {
      TWOeTWOmu_card.getline(line,256);
    }
  TWOeTWOmu_card.getline(line,256);
  TWOeTWOmu_card.getline(line,256);
  TWOeTWOmu_card.getline(line,256);
  TWOeTWOmu_card.getline(line,256);
  TWOeTWOmu_card.getline(line,256);
  //cout << line << endl;
  
  RooRealVar* qqZZ_2e2mu_norm = new RooRealVar("qqZZ_2e2mu_norm","qqZZ_norm",0.);
  RooRealVar* ZX_2e2mu_norm = new RooRealVar("ZX_2e2mu_norm","ZX_norm",0.);
  
  char* chars_array = strtok(line," ");
  int count = 0;
  while(chars_array)
    {
      if (count == 3) qqZZ_2e2mu_norm->setVal(atof(chars_array));
      else if (count == 4) ZX_2e2mu_norm->setVal(atof(chars_array));
      chars_array = strtok(NULL," ");
      count++;
    }
  
  
  TH2* Zjets_2e2mu_2Dhist = Zjets_2e2mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  Zjets_2e2mu_2Dhist->Draw("COLZ");
  gStyle->SetOptStat(0);
  c1->SaveAs("Zjets.C");
  c1->SaveAs("Zjets.root");
  c1->SaveAs("Zjets.eps");
  c1->SaveAs("Zjets.png");
  c1->SaveAs("Zjets.pdf");
  //Zjets_2e2mu_2Dhist->GetXaxis()->SetLimits(220,800);
  Zjets_2e2mu_2Dhist->SetTitle("Zjets 2e2mu;m_{4l};D_{gg}");
  TH1D* Zjets_2e2mu_mass = Zjets_2e2mu_2Dhist->ProjectionX("Zjets_2e2mu_mass");
  TH1D* Zjets_2e2mu_kd = Zjets_2e2mu_2Dhist->ProjectionY("Zjets_2e2mu_kd");
  Zjets_2e2mu_mass->Scale(ZX_2e2mu_norm->getVal()/Zjets_2e2mu_mass->Integral());
  Zjets_2e2mu_kd->Scale(ZX_2e2mu_norm->getVal()/Zjets_2e2mu_kd->Integral());
  Zjets_2e2mu_mass->SetLineColor(1);
  Zjets_2e2mu_mass->SetLineWidth(2);
  Zjets_2e2mu_mass->SetFillColor(kGreen-5);
  Zjets_2e2mu_kd->SetLineColor(1);
  Zjets_2e2mu_kd->SetLineWidth(2);
  Zjets_2e2mu_kd->SetFillColor(kGreen-5);

  TH2* qqZZ_2e2mu_2Dhist = qqZZ_2e2mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  qqZZ_2e2mu_2Dhist->Draw("COLZ");
  gStyle->SetOptStat(0);
  c1->SaveAs("qqZZ.C");
  c1->SaveAs("qqZZ.root");
  c1->SaveAs("qqZZ.eps");
  c1->SaveAs("qqZZ.png");
  c1->SaveAs("qqZZ.pdf");
  //qqZZ_2e2mu_2Dhist->GetXaxis()->SetLimits(220,800);
  qqZZ_2e2mu_2Dhist->SetTitle("qqZZ 2e2mu;m_{4l};D_{gg}");
  TH1D* qqZZ_2e2mu_mass = qqZZ_2e2mu_2Dhist->ProjectionX("qqZZ_2e2mu_mass");
  TH1D* qqZZ_2e2mu_kd = qqZZ_2e2mu_2Dhist->ProjectionY("qqZZ_2e2mu_kd");
  qqZZ_2e2mu_mass->Scale(qqZZ_2e2mu_norm->getVal()/qqZZ_2e2mu_mass->Integral());
  qqZZ_2e2mu_kd->Scale(qqZZ_2e2mu_norm->getVal()/qqZZ_2e2mu_kd->Integral());
  qqZZ_2e2mu_mass->SetLineColor(1);
  qqZZ_2e2mu_mass->SetLineWidth(2);
  qqZZ_2e2mu_mass->SetFillColor(kAzure-9);
  qqZZ_2e2mu_kd->SetLineColor(1);
  qqZZ_2e2mu_kd->SetLineWidth(2);
  qqZZ_2e2mu_kd->SetFillColor(kAzure-9);


  RooAbsPdf* ggZZ_2e2mu = TWOeTWOmu_w->pdf("ggzz");
  RooFormulaVar* ggZZ_2e2mu_norm = TWOeTWOmu_w->function("ggzz_norm");
  TWOeTWOmu_mu->setVal(1.0);
  TWOeTWOmu_GGsm->setVal(1.0);

  TH2* ggZZ_2e2mu_2Dhist = ggZZ_2e2mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  ggZZ_2e2mu_2Dhist->Draw("COLZ");
  gStyle->SetOptStat(0);
  c1->SaveAs("ggZZ_SM.C");
  c1->SaveAs("ggZZ_SM.root");
  c1->SaveAs("ggZZ_SM.eps");
  c1->SaveAs("ggZZ_SM.png");
  c1->SaveAs("ggZZ_SM.pdf");
  //ggZZ_2e2mu_2Dhist->GetXaxis()->SetLimits(220,800);
  ggZZ_2e2mu_2Dhist->SetTitle("ggZZ 2e2mu;m_{4l};D_{gg}");
  TH1D* ggZZ_2e2mu_mass = ggZZ_2e2mu_2Dhist->ProjectionX("ggZZ_2e2mu_mass");
  TH1D* ggZZ_2e2mu_kd = ggZZ_2e2mu_2Dhist->ProjectionY("ggZZ_2e2mu_kd");
  ggZZ_2e2mu_mass->Scale(ggZZ_2e2mu_norm->getVal()/ggZZ_2e2mu_mass->Integral());
  ggZZ_2e2mu_kd->Scale(ggZZ_2e2mu_norm->getVal()/ggZZ_2e2mu_kd->Integral());
  cout << ggZZ_2e2mu_norm->getVal() << " " << ggZZ_2e2mu_mass->Integral() << " " << ggZZ_2e2mu_kd->Integral() << endl;
  ggZZ_2e2mu_mass->SetLineColor(kOrange+10);
  ggZZ_2e2mu_mass->SetLineWidth(2);
  ggZZ_2e2mu_mass->SetFillColor(0);
  ggZZ_2e2mu_kd->SetLineColor(kOrange+10);
  ggZZ_2e2mu_kd->SetLineWidth(2);
  ggZZ_2e2mu_kd->SetFillColor(0);

  TWOeTWOmu_mu->setVal(1.0);
  TWOeTWOmu_GGsm->setVal(25.0);

  TH2* BSM_2e2mu_2Dhist = ggZZ_2e2mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  BSM_2e2mu_2Dhist->Draw("COLZ");
  gStyle->SetOptStat(0);
  c1->SaveAs("ggZZ_BSM.C");
  c1->SaveAs("ggZZ_BSM.root");
  c1->SaveAs("ggZZ_BSM.eps");
  c1->SaveAs("ggZZ_BSM.png");
  c1->SaveAs("ggZZ_BSM.pdf");
  //BSM_2e2mu_2Dhist->GetXaxis()->SetLimits(220,800);
  BSM_2e2mu_2Dhist->SetTitle("BSM 2e2mu;m_{4l};D_{gg}");
  TH1D* BSM_2e2mu_mass = BSM_2e2mu_2Dhist->ProjectionX("BSM_2e2mu_mass");
  TH1D* BSM_2e2mu_kd = BSM_2e2mu_2Dhist->ProjectionY("BSM_2e2mu_kd");
  BSM_2e2mu_mass->Scale(ggZZ_2e2mu_norm->getVal()/BSM_2e2mu_mass->Integral());
  BSM_2e2mu_kd->Scale(ggZZ_2e2mu_norm->getVal()/BSM_2e2mu_kd->Integral());
  cout << ggZZ_2e2mu_norm->getVal() << " " << BSM_2e2mu_mass->Integral() << " " << BSM_2e2mu_kd->Integral() << endl;
  BSM_2e2mu_mass->SetLineColor(kOrange+10);
  BSM_2e2mu_mass->SetLineWidth(2);
  BSM_2e2mu_mass->SetLineStyle(2);
  BSM_2e2mu_mass->SetFillColor(0);
  BSM_2e2mu_kd->SetLineColor(kOrange+10);
  BSM_2e2mu_kd->SetLineWidth(2);
  BSM_2e2mu_kd->SetLineStyle(2);
  BSM_2e2mu_kd->SetFillColor(0);

  RooAbsPdf* VBF_2e2mu = TWOeTWOmu_w->pdf("vbf_offshell");
  RooFormulaVar* VBF_2e2mu_norm = TWOeTWOmu_w->function("vbf_offshell_norm");
  TWOeTWOmu_mu->setVal(1.0);
  TWOeTWOmu_GGsm->setVal(1.0);

  TH2* VBF_2e2mu_2Dhist = VBF_2e2mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  VBF_2e2mu_2Dhist->Draw("COLZ");
  gStyle->SetOptStat(0);
  c1->SaveAs("VBF_SM.C");
  c1->SaveAs("VBF_SM.root");
  c1->SaveAs("VBF_SM.eps");
  c1->SaveAs("VBF_SM.png");
  c1->SaveAs("VBF_SM.pdf");
  //VBF_2e2mu_2Dhist->GetXaxis()->SetLimits(220,800);
  VBF_2e2mu_2Dhist->SetTitle("VBF 2e2mu;m_{4l};D_{gg}");
  TH1D* VBF_2e2mu_mass = VBF_2e2mu_2Dhist->ProjectionX("VBF_2e2mu_mass");
  TH1D* VBF_2e2mu_kd = VBF_2e2mu_2Dhist->ProjectionY("VBF_2e2mu_kd");
  VBF_2e2mu_mass->Scale(VBF_2e2mu_norm->getVal()/VBF_2e2mu_mass->Integral());
  VBF_2e2mu_kd->Scale(VBF_2e2mu_norm->getVal()/VBF_2e2mu_kd->Integral());
  cout << VBF_2e2mu_norm->getVal() << " " << VBF_2e2mu_mass->Integral() << " " << VBF_2e2mu_kd->Integral() << endl;
  VBF_2e2mu_mass->SetLineColor(kViolet-2);
  VBF_2e2mu_mass->SetLineWidth(2);
  VBF_2e2mu_mass->SetFillColor(0);
  VBF_2e2mu_kd->SetLineColor(kViolet-2);
  VBF_2e2mu_kd->SetLineWidth(2);
  VBF_2e2mu_kd->SetFillColor(0);

  TWOeTWOmu_mu->setVal(1.0);
  TWOeTWOmu_GGsm->setVal(25.0);

  TH2* BSM_VBF_2e2mu_2Dhist = VBF_2e2mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  BSM_VBF_2e2mu_2Dhist->Draw("COLZ");
  gStyle->SetOptStat(0);
  c1->SaveAs("VBF_BSM.C");
  c1->SaveAs("VBF_BSM.root");
  c1->SaveAs("VBF_BSM.eps");
  c1->SaveAs("VBF_BSM.png");
  c1->SaveAs("VBF_BSM.pdf");
  //BSM_VBF_2e2mu_2Dhist->GetXaxis()->SetLimits(220,800);
  BSM_VBF_2e2mu_2Dhist->SetTitle("BSM VBF 2e2mu;m_{4l};D_{gg}");
  TH1D* BSM_VBF_2e2mu_mass = BSM_VBF_2e2mu_2Dhist->ProjectionX("BSM_VBF_2e2mu_mass");
  TH1D* BSM_VBF_2e2mu_kd = BSM_VBF_2e2mu_2Dhist->ProjectionY("BSM_VBF_2e2mu_kd");
  BSM_VBF_2e2mu_mass->Scale(VBF_2e2mu_norm->getVal()/BSM_VBF_2e2mu_mass->Integral());
  BSM_VBF_2e2mu_kd->Scale(VBF_2e2mu_norm->getVal()/BSM_VBF_2e2mu_kd->Integral());
  cout << VBF_2e2mu_norm->getVal() << " " << BSM_VBF_2e2mu_mass->Integral() << " " << BSM_VBF_2e2mu_kd->Integral() << endl;
  BSM_VBF_2e2mu_mass->SetLineColor(kViolet-2);
  BSM_VBF_2e2mu_mass->SetLineWidth(2);
  BSM_VBF_2e2mu_mass->SetLineStyle(2);
  BSM_VBF_2e2mu_mass->SetFillColor(0);
  BSM_VBF_2e2mu_kd->SetLineColor(kViolet-2);
  BSM_VBF_2e2mu_kd->SetLineWidth(2);
  BSM_VBF_2e2mu_kd->SetLineStyle(2);
  BSM_VBF_2e2mu_kd->SetFillColor(0);

  TFile* FOURmu_f = new TFile("hzz4l_4muS_8TeV.input.root","OPEN");
  RooWorkspace* FOURmu_w = FOURmu_f->Get("w");
  RooAbsPdf* Zjets_4mu = FOURmu_w->pdf("bkg_zjets");
  RooAbsPdf* qqZZ_4mu = FOURmu_w->pdf("bkg_qqzz");
  RooRealVar* FOURmu_mu = FOURmu_w->var("CMS_zz4l_mu");
  RooRealVar* FOURmu_GGsm = FOURmu_w->var("CMS_zz4l_GGsm");


  RooDataSet* data_4mu = FOURmu_w->data("data_obs")->Clone("data_4mu");
  if(sumChannels)data->append(*data_4mu);

  ifstream FOURmu_card;
  FOURmu_card.open("hzz4l_4muS_8TeV.txt");
  char line[256];
  for (int n_line=0; n_line < 7; n_line++)
    {
      FOURmu_card.getline(line,256);
    }
  FOURmu_card.getline(line,256);
  for(int n_line=0; n_line < 3; n_line++)
    {
      FOURmu_card.getline(line,256);
    }
  FOURmu_card.getline(line,256);
  FOURmu_card.getline(line,256);
  FOURmu_card.getline(line,256);
  FOURmu_card.getline(line,256);
  FOURmu_card.getline(line,256);
  //cout << line << endl;
  
  RooRealVar* qqZZ_4mu_norm = new RooRealVar("qqZZ_4mu_norm","qqZZ_norm",0.);
  RooRealVar* ZX_4mu_norm = new RooRealVar("ZX_4mu_norm","ZX_norm",0.);
  
  char* chars_array = strtok(line," ");
  int count = 0;
  while(chars_array)
    {
      if (count == 3) qqZZ_4mu_norm->setVal(atof(chars_array));
      else if (count == 4) ZX_4mu_norm->setVal(atof(chars_array));
      chars_array = strtok(NULL," ");
      count++;
    }
  
  
  TH2* Zjets_4mu_2Dhist = Zjets_4mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //Zjets_4mu_2Dhist->GetXaxis()->SetLimits(220,800);
  Zjets_4mu_2Dhist->SetTitle("Zjets 4mu;m_{4l};D_{gg}");
  TH1D* Zjets_4mu_mass = Zjets_4mu_2Dhist->ProjectionX("Zjets_4mu_mass");
  TH1D* Zjets_4mu_kd = Zjets_4mu_2Dhist->ProjectionY("Zjets_4mu_kd");
  Zjets_4mu_mass->Scale(ZX_4mu_norm->getVal()/Zjets_4mu_mass->Integral());
  Zjets_4mu_kd->Scale(ZX_4mu_norm->getVal()/Zjets_4mu_kd->Integral());
  Zjets_4mu_mass->SetLineColor(1);
  Zjets_4mu_mass->SetLineWidth(2);
  Zjets_4mu_mass->SetFillColor(kGreen-5);
  Zjets_4mu_kd->SetLineColor(1);
  Zjets_4mu_kd->SetLineWidth(2);
  Zjets_4mu_kd->SetFillColor(kGreen-5);

  TH2* qqZZ_4mu_2Dhist = qqZZ_4mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //qqZZ_4mu_2Dhist->GetXaxis()->SetLimits(220,800);
  qqZZ_4mu_2Dhist->SetTitle("qqZZ 4mu;m_{4l};D_{gg}");
  TH1D* qqZZ_4mu_mass = qqZZ_4mu_2Dhist->ProjectionX("qqZZ_4mu_mass");
  TH1D* qqZZ_4mu_kd = qqZZ_4mu_2Dhist->ProjectionY("qqZZ_4mu_kd");
  qqZZ_4mu_mass->Scale(qqZZ_4mu_norm->getVal()/qqZZ_4mu_mass->Integral());
  qqZZ_4mu_kd->Scale(qqZZ_4mu_norm->getVal()/qqZZ_4mu_kd->Integral());
  qqZZ_4mu_mass->SetLineColor(1);
  qqZZ_4mu_mass->SetLineWidth(2);
  qqZZ_4mu_mass->SetFillColor(kAzure-9);
  qqZZ_4mu_kd->SetLineColor(1);
  qqZZ_4mu_kd->SetLineWidth(2);
  qqZZ_4mu_kd->SetFillColor(kAzure-9);

  RooAbsPdf* ggZZ_4mu = FOURmu_w->pdf("ggzz");
  RooFormulaVar* ggZZ_4mu_norm = FOURmu_w->function("ggzz_norm");
  FOURmu_mu->setVal(1.0);
  FOURmu_GGsm->setVal(1.0);

  TH2* ggZZ_4mu_2Dhist = ggZZ_4mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //ggZZ_4mu_2Dhist->GetXaxis()->SetLimits(220,800);
  ggZZ_4mu_2Dhist->SetTitle("ggZZ 4mu;m_{4l};D_{gg}");
  TH1D* ggZZ_4mu_mass = ggZZ_4mu_2Dhist->ProjectionX("ggZZ_4mu_mass");
  TH1D* ggZZ_4mu_kd = ggZZ_4mu_2Dhist->ProjectionY("ggZZ_4mu_kd");
  ggZZ_4mu_mass->Scale(ggZZ_4mu_norm->getVal()/ggZZ_4mu_mass->Integral());
  ggZZ_4mu_kd->Scale(ggZZ_4mu_norm->getVal()/ggZZ_4mu_kd->Integral());
  cout << ggZZ_4mu_norm->getVal() << " " << ggZZ_4mu_mass->Integral() << " " << ggZZ_4mu_kd->Integral() << endl;
  ggZZ_4mu_mass->SetLineColor(kOrange+10);
  ggZZ_4mu_mass->SetLineWidth(2);
  ggZZ_4mu_mass->SetFillColor(0);
  ggZZ_4mu_kd->SetLineColor(kOrange+10);
  ggZZ_4mu_kd->SetLineWidth(2);
  ggZZ_4mu_kd->SetFillColor(0);

  FOURmu_mu->setVal(1.0);
  FOURmu_GGsm->setVal(25.0);

  TH2* BSM_4mu_2Dhist = ggZZ_4mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //BSM_4mu_2Dhist->GetXaxis()->SetLimits(220,800);
  BSM_4mu_2Dhist->SetTitle("BSM 4mu;m_{4l};D_{gg}");
  TH1D* BSM_4mu_mass = BSM_4mu_2Dhist->ProjectionX("BSM_4mu_mass");
  TH1D* BSM_4mu_kd = BSM_4mu_2Dhist->ProjectionY("BSM_4mu_kd");
  BSM_4mu_mass->Scale(ggZZ_4mu_norm->getVal()/BSM_4mu_mass->Integral());
  BSM_4mu_kd->Scale(ggZZ_4mu_norm->getVal()/BSM_4mu_kd->Integral());
  cout << ggZZ_4mu_norm->getVal() << " " << BSM_4mu_mass->Integral() << " " << BSM_4mu_kd->Integral() << endl;
  BSM_4mu_mass->SetLineColor(kOrange+10);
  BSM_4mu_mass->SetLineWidth(2);
  BSM_4mu_mass->SetLineStyle(2);
  BSM_4mu_mass->SetFillColor(0);
  BSM_4mu_kd->SetLineColor(kOrange+10);
  BSM_4mu_kd->SetLineWidth(2);
  BSM_4mu_kd->SetLineStyle(2);
  BSM_4mu_kd->SetFillColor(0);

  RooAbsPdf* VBF_4mu = FOURmu_w->pdf("vbf_offshell");
  RooFormulaVar* VBF_4mu_norm = FOURmu_w->function("vbf_offshell_norm");
  FOURmu_mu->setVal(1.0);
  FOURmu_GGsm->setVal(1.0);

  TH2* VBF_4mu_2Dhist = VBF_4mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //VBF_4mu_2Dhist->GetXaxis()->SetLimits(220,800);
  VBF_4mu_2Dhist->SetTitle("VBF 4mu;m_{4l};D_{gg}");
  TH1D* VBF_4mu_mass = VBF_4mu_2Dhist->ProjectionX("VBF_4mu_mass");
  TH1D* VBF_4mu_kd = VBF_4mu_2Dhist->ProjectionY("VBF_4mu_kd");
  VBF_4mu_mass->Scale(VBF_4mu_norm->getVal()/VBF_4mu_mass->Integral());
  VBF_4mu_kd->Scale(VBF_4mu_norm->getVal()/VBF_4mu_kd->Integral());
  cout << VBF_4mu_norm->getVal() << " " << VBF_4mu_mass->Integral() << " " << VBF_4mu_kd->Integral() << endl;
  VBF_4mu_mass->SetLineColor(kViolet-2);
  VBF_4mu_mass->SetLineWidth(2);
  VBF_4mu_mass->SetFillColor(0);
  VBF_4mu_kd->SetLineColor(kViolet-2);
  VBF_4mu_kd->SetLineWidth(2);
  VBF_4mu_kd->SetFillColor(0);

  FOURmu_mu->setVal(1.0);
  FOURmu_GGsm->setVal(25.0);

  TH2* BSM_VBF_4mu_2Dhist = VBF_4mu->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //BSM_VBF_4mu_2Dhist->GetXaxis()->SetLimits(220,800);
  BSM_VBF_4mu_2Dhist->SetTitle("BSM VBF 4mu;m_{4l};D_{gg}");
  TH1D* BSM_VBF_4mu_mass = BSM_VBF_4mu_2Dhist->ProjectionX("BSM_VBF_4mu_mass");
  TH1D* BSM_VBF_4mu_kd = BSM_VBF_4mu_2Dhist->ProjectionY("BSM_VBF_4mu_kd");
  BSM_VBF_4mu_mass->Scale(VBF_4mu_norm->getVal()/BSM_VBF_4mu_mass->Integral());
  BSM_VBF_4mu_kd->Scale(VBF_4mu_norm->getVal()/BSM_VBF_4mu_kd->Integral());
  cout << VBF_4mu_norm->getVal() << " " << BSM_VBF_4mu_mass->Integral() << " " << BSM_VBF_4mu_kd->Integral() << endl;
  BSM_VBF_4mu_mass->SetLineColor(kViolet-2);
  BSM_VBF_4mu_mass->SetLineWidth(2);
  BSM_VBF_4mu_mass->SetLineStyle(2);
  BSM_VBF_4mu_mass->SetFillColor(0);
  BSM_VBF_4mu_kd->SetLineColor(kViolet-2);
  BSM_VBF_4mu_kd->SetLineWidth(2);
  BSM_VBF_4mu_kd->SetLineStyle(2);
  BSM_VBF_4mu_kd->SetFillColor(0);

  TFile* FOURe_f = new TFile("hzz4l_4eS_8TeV.input.root","OPEN");
  RooWorkspace* FOURe_w = FOURe_f->Get("w");
  RooAbsPdf* Zjets_4e = FOURe_w->pdf("bkg_zjets");
  RooAbsPdf* qqZZ_4e = FOURe_w->pdf("bkg_qqzz");
  RooRealVar* FOURe_mu = FOURe_w->var("CMS_zz4l_mu");
  RooRealVar* FOURe_GGsm = FOURe_w->var("CMS_zz4l_GGsm");
  RooRealVar* FOURe_mass = FOURe_w->var("CMS_zz4l_widthMass");

  RooDataSet* data_4e = FOURe_w->data("data_obs")->Clone("data_4e");
  if(sumChannels)data->append(*data_4e);
  RooPlot* mass_data = FOURe_mass.frame();
  RooPlot* kd_data = TWOeTWOmu_kd.frame();
  data_4e->plotOn(mass_data);
  data_4e->plotOn(kd_data);

  ifstream FOURe_card;
  FOURe_card.open("hzz4l_4eS_8TeV.txt");
  char line[256];
  for (int n_line=0; n_line < 7; n_line++)
    {
      FOURe_card.getline(line,256);
    }
  FOURe_card.getline(line,256);
  for(int n_line=0; n_line < 3; n_line++)
    {
      FOURe_card.getline(line,256);
    }
  FOURe_card.getline(line,256);
  FOURe_card.getline(line,256);
  FOURe_card.getline(line,256);
  FOURe_card.getline(line,256);
  FOURe_card.getline(line,256);
  //cout << line << endl;
  
  RooRealVar* qqZZ_4e_norm = new RooRealVar("qqZZ_4e_norm","qqZZ_norm",0.);
  RooRealVar* ZX_4e_norm = new RooRealVar("ZX_4e_norm","ZX_norm",0.);
  
  char* chars_array = strtok(line," ");
  int count = 0;
  while(chars_array)
    {
      if (count == 3) qqZZ_4e_norm->setVal(atof(chars_array));
      else if (count == 4) ZX_4e_norm->setVal(atof(chars_array));
      chars_array = strtok(NULL," ");
      count++;
    }
  
  
  TH2* Zjets_4e_2Dhist = Zjets_4e->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //Zjets_4e_2Dhist->GetXaxis()->SetLimits(220,800);
  Zjets_4e_2Dhist->SetTitle("Zjets 4e;m_{4l};D_{gg}");
  TH1D* Zjets_4e_mass = Zjets_4e_2Dhist->ProjectionX("Zjets_4e_mass");
  TH1D* Zjets_4e_kd = Zjets_4e_2Dhist->ProjectionY("Zjets_4e_kd");
  Zjets_4e_mass->Scale(ZX_4e_norm->getVal()/Zjets_4e_mass->Integral());
  Zjets_4e_kd->Scale(ZX_4e_norm->getVal()/Zjets_4e_kd->Integral());
  Zjets_4e_mass->SetLineColor(1);
  Zjets_4e_mass->SetLineWidth(2);
  Zjets_4e_mass->SetFillColor(kGreen-5);
  Zjets_4e_kd->SetLineColor(1);
  Zjets_4e_kd->SetLineWidth(2);
  Zjets_4e_kd->SetFillColor(kGreen-5);

  TH2* qqZZ_4e_2Dhist = qqZZ_4e->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //qqZZ_4e_2Dhist->GetXaxis()->SetLimits(220,800);
  qqZZ_4e_2Dhist->SetTitle("qqZZ 4e;m_{4l};D_{gg}");
  TH1D* qqZZ_4e_mass = qqZZ_4e_2Dhist->ProjectionX("qqZZ_4e_mass");
  TH1D* qqZZ_4e_kd = qqZZ_4e_2Dhist->ProjectionY("qqZZ_4e_kd");
  qqZZ_4e_mass->Scale(qqZZ_4e_norm->getVal()/qqZZ_4e_mass->Integral());
  qqZZ_4e_kd->Scale(qqZZ_4e_norm->getVal()/qqZZ_4e_kd->Integral());
  qqZZ_4e_mass->SetLineColor(1);
  qqZZ_4e_mass->SetLineWidth(2);
  qqZZ_4e_mass->SetFillColor(kAzure-9);
  qqZZ_4e_kd->SetLineColor(1);
  qqZZ_4e_kd->SetLineWidth(2);
  qqZZ_4e_kd->SetFillColor(kAzure-9);

  RooAbsPdf* ggZZ_4e = FOURe_w->pdf("ggzz");
  RooFormulaVar* ggZZ_4e_norm = FOURe_w->function("ggzz_norm");
  FOURe_mu->setVal(1.0);
  FOURe_GGsm->setVal(1.0);

  TH2* ggZZ_4e_2Dhist = ggZZ_4e->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //ggZZ_4e_2Dhist->GetXaxis()->SetLimits(220,800);
  ggZZ_4e_2Dhist->SetTitle("ggZZ 4e;m_{4l};D_{gg}");
  TH1D* ggZZ_4e_mass = ggZZ_4e_2Dhist->ProjectionX("ggZZ_4e_mass");
  TH1D* ggZZ_4e_kd = ggZZ_4e_2Dhist->ProjectionY("ggZZ_4e_kd");
  ggZZ_4e_mass->Scale(ggZZ_4e_norm->getVal()/ggZZ_4e_mass->Integral());
  ggZZ_4e_kd->Scale(ggZZ_4e_norm->getVal()/ggZZ_4e_kd->Integral());
  cout << ggZZ_4e_norm->getVal() << " " << ggZZ_4e_mass->Integral() << " " << ggZZ_4e_kd->Integral() << endl;
  ggZZ_4e_mass->SetLineColor(kOrange+10);
  ggZZ_4e_mass->SetLineWidth(2);
  ggZZ_4e_mass->SetFillColor(0);
  ggZZ_4e_kd->SetLineColor(kOrange+10);
  ggZZ_4e_kd->SetLineWidth(2);
  ggZZ_4e_kd->SetFillColor(0);

  FOURe_mu->setVal(1.0);
  FOURe_GGsm->setVal(25.0);

  TH2* BSM_4e_2Dhist = ggZZ_4e->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //BSM_4e_2Dhist->GetXaxis()->SetLimits(220,800);
  BSM_4e_2Dhist->SetTitle("BSM 4e;m_{4l};D_{gg}");
  TH1D* BSM_4e_mass = BSM_4e_2Dhist->ProjectionX("BSM_4e_mass");
  TH1D* BSM_4e_kd = BSM_4e_2Dhist->ProjectionY("BSM_4e_kd");
  BSM_4e_mass->Scale(ggZZ_4e_norm->getVal()/BSM_4e_mass->Integral());
  BSM_4e_kd->Scale(ggZZ_4e_norm->getVal()/BSM_4e_kd->Integral());
  cout << ggZZ_4e_norm->getVal() << " " << BSM_4e_mass->Integral() << " " << BSM_4e_kd->Integral() << endl;
  BSM_4e_mass->SetLineColor(kOrange+10);
  BSM_4e_mass->SetLineWidth(2);
  BSM_4e_mass->SetLineStyle(2);
  BSM_4e_mass->SetFillColor(0);
  BSM_4e_kd->SetLineColor(kOrange+10);
  BSM_4e_kd->SetLineWidth(2);
  BSM_4e_kd->SetLineStyle(2);
  BSM_4e_kd->SetFillColor(0);

  RooAbsPdf* VBF_4e = FOURe_w->pdf("vbf_offshell");
  RooFormulaVar* VBF_4e_norm = FOURe_w->function("vbf_offshell_norm");
  FOURe_mu->setVal(1.0);
  FOURe_GGsm->setVal(1.0);

  TH2* VBF_4e_2Dhist = VBF_4e->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //VBF_4e_2Dhist->GetXaxis()->SetLimits(220,800);
  VBF_4e_2Dhist->SetTitle("VBF 4e;m_{4l};D_{gg}");
  TH1D* VBF_4e_mass = VBF_4e_2Dhist->ProjectionX("VBF_4e_mass");
  TH1D* VBF_4e_kd = VBF_4e_2Dhist->ProjectionY("VBF_4e_kd");
  VBF_4e_mass->Scale(VBF_4e_norm->getVal()/VBF_4e_mass->Integral());
  VBF_4e_kd->Scale(VBF_4e_norm->getVal()/VBF_4e_kd->Integral());
  cout << VBF_4e_norm->getVal() << " " << VBF_4e_mass->Integral() << " " << VBF_4e_kd->Integral() << endl;
  VBF_4e_mass->SetLineColor(kViolet-2);
  VBF_4e_mass->SetLineWidth(2);
  VBF_4e_mass->SetFillColor(0);
  VBF_4e_kd->SetLineColor(kViolet-2);
  VBF_4e_kd->SetLineWidth(2);
  VBF_4e_kd->SetFillColor(0);

  FOURe_mu->setVal(1.0);
  FOURe_GGsm->setVal(25.0);

  TH2* BSM_VBF_4e_2Dhist = VBF_4e->createHistogram("CMS_zz4l_widthMass,CMS_zz4l_widthKD");
  //BSM_VBF_4e_2Dhist->GetXaxis()->SetLimits(220,800);
  BSM_VBF_4e_2Dhist->SetTitle("BSM VBF 4e;m_{4l};D_{gg}");
  TH1D* BSM_VBF_4e_mass = BSM_VBF_4e_2Dhist->ProjectionX("BSM_VBF_4e_mass");
  TH1D* BSM_VBF_4e_kd = BSM_VBF_4e_2Dhist->ProjectionY("BSM_VBF_4e_kd");
  BSM_VBF_4e_mass->Scale(VBF_4e_norm->getVal()/BSM_VBF_4e_mass->Integral());
  BSM_VBF_4e_kd->Scale(VBF_4e_norm->getVal()/BSM_VBF_4e_kd->Integral());
  cout << VBF_4e_norm->getVal() << " " << BSM_VBF_4e_mass->Integral() << " " << BSM_VBF_4e_kd->Integral() << endl;
  BSM_VBF_4e_mass->SetLineColor(kViolet-2);
  BSM_VBF_4e_mass->SetLineWidth(2);
  BSM_VBF_4e_mass->SetLineStyle(2);
  BSM_VBF_4e_mass->SetFillColor(0);
  BSM_VBF_4e_kd->SetLineColor(kViolet-2);
  BSM_VBF_4e_kd->SetLineWidth(2);
  BSM_VBF_4e_kd->SetLineStyle(2);
  BSM_VBF_4e_kd->SetFillColor(0);
  if(sumChannels){
    qqZZ_2e2mu_mass->Add(qqZZ_4mu_mass);
    qqZZ_2e2mu_mass->Add(qqZZ_4e_mass);
    qqZZ_2e2mu_kd->Add(qqZZ_4mu_kd);
    qqZZ_2e2mu_kd->Add(qqZZ_4e_kd);

    ggZZ_2e2mu_mass->Add(ggZZ_4mu_mass);
    ggZZ_2e2mu_mass->Add(ggZZ_4e_mass);
    ggZZ_2e2mu_kd->Add(ggZZ_4mu_kd);
    ggZZ_2e2mu_kd->Add(ggZZ_4e_kd);

    BSM_2e2mu_mass->Add(BSM_4mu_mass);
    BSM_2e2mu_mass->Add(BSM_4e_mass);
    BSM_2e2mu_kd->Add(BSM_4mu_kd);
    BSM_2e2mu_kd->Add(BSM_4e_kd);

    VBF_2e2mu_mass->Add(VBF_4mu_mass);
    VBF_2e2mu_mass->Add(VBF_4e_mass);
    VBF_2e2mu_kd->Add(VBF_4mu_kd);
    VBF_2e2mu_kd->Add(VBF_4e_kd);

    BSM_VBF_2e2mu_mass->Add(BSM_VBF_4mu_mass);
    BSM_VBF_2e2mu_mass->Add(BSM_VBF_4e_mass);
    BSM_VBF_2e2mu_kd->Add(BSM_VBF_4mu_kd);
    BSM_VBF_2e2mu_kd->Add(BSM_VBF_4e_kd);
    
    Zjets_2e2mu_mass->Add(Zjets_4mu_mass);
    Zjets_2e2mu_mass->Add(Zjets_4e_mass);
    Zjets_2e2mu_kd->Add(Zjets_4mu_kd);
    Zjets_2e2mu_kd->Add(Zjets_4e_kd);
  }
  THStack* kd_stack = new THStack("KDStack","KDStack");
  kd_stack->SetTitle(Form(";D_{gg};Events/%1.4f",Zjets_4e_kd->GetXaxis()->GetBinWidth(2)));
  kd_stack->Add(Zjets_4e_kd);
  kd_stack->Add(qqZZ_4e_kd);
  kd_stack->Add(ggZZ_4e_kd);
  kd_stack->Add(VBF_4e_kd);

  THStack* kd_stack2 = new THStack("KDStack2","KDStack2");
  kd_stack2->SetTitle(Form(";D_{gg};Events/%1.4f",Zjets_4e_kd->GetXaxis()->GetBinWidth(2)));
  kd_stack2->Add(Zjets_4e_kd);
  kd_stack2->Add(qqZZ_4e_kd);
  kd_stack2->Add(BSM_4e_kd);
  kd_stack2->Add(BSM_VBF_4e_kd);

  TPaveText *pt = new TPaveText(0.0992462,0.90544,0.899497,0.943005,"NDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  TText *text = pt->AddText(0.01,0.5,"CMS Preliminary");
  //text = pt->AddText(0.2,0.6,Form("#sqrt{s} = 7 TeV, L = %.1f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi7TeV,lumi8TeV));
  text = pt->AddText(0.55,0.6,Form("#sqrt{s} = 8 TeV, L = %.1f fb^{-1}",lumi8TeV)); 

  //  kd_stack2->Draw();
  kd_data->Draw();
  kd_stack2->Draw("same");
  kd_stack->Draw("same");
  kd_data->Draw("same");
  kd_stack2->GetYaxis()->SetRange(0,15);
  gPad->RedrawAxis();
  pt->Draw();

  TLegend *leg2 = new TLegend(.561558,.637306,.894472,.896373);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);

  leg2->AddEntry(Zjets_4e_mass,"Z+X","f");
  leg2->AddEntry(qqZZ_4e_mass,"qqZZ","f");
  leg2->AddEntry(ggZZ_4e_mass,"ggZZ_{SM} (#Gamma = #Gamma_{SM})","f");
  leg2->AddEntry(BSM_4e_mass,"ggZZ_{BSM} (#Gamma = 25 #times #Gamma_{SM})","f");
  leg2->AddEntry(VBF_4e_mass,"VBF_{SM} (#Gamma = #Gamma_{SM})","f");
  leg2->AddEntry(BSM_VBF_4e_mass,"VBF_{BSM} (#Gamma = 25 #times #Gamma_{SM})","f");
  leg2->Draw();

  c1->SaveAs("can_kd_ggsm4e.C");
  c1->SaveAs("can_kd_ggsm4e.root");
  c1->SaveAs("can_kd_ggsm4e.eps");
  c1->SaveAs("can_kd_ggsm4e.png");
  c1->SaveAs("can_kd_ggsm4e.pdf");

  THStack* mass_stack = new THStack("MassStack","MassStack");
  mass_stack->SetTitle(Form(";m_{4l};Events/%2.0fGeV",Zjets_4e_mass->GetXaxis()->GetBinWidth(20)));
  mass_stack->Add(Zjets_4e_mass);
  mass_stack->Add(qqZZ_4e_mass);
  mass_stack->Add(ggZZ_4e_mass);
  mass_stack->Add(VBF_4e_mass);

  THStack* mass_stack2 = new THStack("MassStack2","MassStack2");
  mass_stack2->SetTitle(Form(";m_{4l};Events/%2.0fGeV",Zjets_4e_mass->GetXaxis()->GetBinWidth(20)));
  mass_stack2->Add(Zjets_4e_mass);
  mass_stack2->Add(qqZZ_4e_mass);
  mass_stack2->Add(BSM_4e_mass);
  mass_stack2->Add(BSM_VBF_4e_mass);

  
  mass_stack2->Draw();
  mass_stack->Draw("same");
  mass_data->Draw("same");
  gPad->RedrawAxis();
  pt->Draw();

  TLegend *leg = new TLegend(.418342,.645078,.817839,.84456);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);

  leg->AddEntry(Zjets_4e_mass,"Z+X","f");
  leg->AddEntry(qqZZ_4e_mass,"qqZZ","f");
  leg->AddEntry(ggZZ_4e_mass,"ggZZ_{SM} (#Gamma = #Gamma_{SM})","f");
  leg->AddEntry(BSM_4e_mass,"ggZZ_{BSM} (#Gamma = 25 #times #Gamma_{SM})","f");
  leg->AddEntry(VBF_4e_mass,"VBF_{SM} (#Gamma = #Gamma_{SM})","f");
  leg->AddEntry(BSM_VBF_4e_mass,"VBF_{BSM} (#Gamma = 25 #times #Gamma_{SM})","f"); 

  leg->Draw();

  c1->SaveAs("can_mass_ggsm4e.C");
  c1->SaveAs("can_mass_ggsm4e.root");
  c1->SaveAs("can_mass_ggsm4e.eps");
  c1->SaveAs("can_mass_ggsm4e.png");
  c1->SaveAs("can_mass_ggsm4e.pdf");

  gPad->SetLogy();
  gPad->Update();
  gPad->RedrawAxis();

  c1->SaveAs("can_mass_ggsm_log4e.C");
  c1->SaveAs("can_mass_ggsm_log4e.root");
  c1->SaveAs("can_mass_ggsm_log4e.eps");
  c1->SaveAs("can_mass_ggsm_log4e.png");
  c1->SaveAs("can_mass_ggsm_log4e.pdf");


}
