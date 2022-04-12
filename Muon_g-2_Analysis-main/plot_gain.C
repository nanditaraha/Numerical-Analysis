{

  
  TMultiGraph *mg = new TMultiGraph();
  //units of drop / energy. p0=4.267e-6 -> 100p0(4), 50p0(5.5), 10p0(5),5p0(6.5),p0(6)
  TFile *f1 = TFile::Open("varyDrop1_m100_c500_f16_80_5.root");
  TFile *f2 = TFile::Open("varyDrop5_m100_c500_f16_80_5.root");
  TFile *f3 = TFile::Open("varyDrop10_m100_c500_f16_80_5.root");
  TFile *f4 = TFile::Open("varyDrop15_m100_c500_f16_80_5.root");
  TFile *f5 = TFile::Open("varyDrop25_m100_c500_f16_80_5.root");
  //For 100 muons
  TH1D *h1 = (TH1D*)f1->Get("Laser");
  TH1D *h2 = (TH1D*)f2->Get("Laser");
  TH1D *h3 = (TH1D*)f3->Get("Laser");
  TH1D *h4 = (TH1D*)f4->Get("Laser");
  TH1D *h5 = (TH1D*)f5->Get("Laser");

  double x1_100 = h1->GetBinContent(30);
  double gain_min_100[5]={x1_100, 0.994,0.988,0.955,0.950};
  double  drop_100[5]={1,5,10,15,25};
  
  //For 50 muons
  double gain_min_50[5]={0.999,0.997,0.994,0.973,0.958};
  //doubledrop_50={1,5,10,50,100};

  //For 10 muons
  double gain_min_10[5]={0.9994,0.9991,0.9988,0.9945,0.9889};

  mg->SetTitle("Laser interval 80#mus and step 5#mus ;Drop per 100 MeV for a pulse in units of p_{0} = 4.26 #times 10^{-6};1-G_{min}");

  TGraph *gr1 = new TGraph(5,drop_100,gain_min_100);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  mg->Add(gr1);

  TGraph *gr2 = new TGraph(5,drop_100,gain_min_50);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(20);
  mg->Add(gr2);

  TGraph *gr3 = new TGraph(5,drop_100,gain_min_10);
  gr3->SetMarkerColor(kMagenta);
  gr3->SetMarkerStyle(22);

  
  mg->Add(gr3);
  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  //leg->SetHeader("","C"); 
  leg->AddEntry(gr1,"100 muons","lep");
  leg->AddEntry(gr2,"50 muons","lep");
  leg->AddEntry(gr3,"10 muons","lep");
  
  mg->Draw("apl");
  leg->Draw();
  
}
