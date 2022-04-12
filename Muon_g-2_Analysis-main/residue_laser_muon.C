{
  
  TFile *f1 = TFile::Open("try_640_5.root");//only for 640 as I dont have enough stats. 
  //TFile *f1 = TFile::Open("noWiggle_fac10_g320_m100.root");
  
  TH1D *h = new TH1D("h","",700/5,-2.5,700-2.5);
  TF1 *f = new TF1("f","pol9",0,200);
  Laser->Draw();//Must do this or the laser histo loses all its value
  Muons->Fit("f","RQ",0,200);
  double x[140]={0};
  for (int i=0;i<140;i++){
    if(Muons->GetBinContent(i*5-1)!=0){
      h->SetBinContent(i,Laser->GetBinContent(i*5)-Muons->GetBinContent(i*5-1));
      //h->SetBinError(i,sqrt(Laser->GetBinError(i*5)*Laser->GetBinError(i*5)+Muons->GetBinError(i*5-1)*Muons->GetBinError(i*5-1)));
      h->SetBinError(i,9.3e-6);
      //cout<<"i="<<i<<" Bin value="<<h->GetBinContent(i)<<" Laser Value="<<Laser->GetBinContent(i)<<" Muons="<<Muons->GetBinContent(i)<<"\n";
    }
    //if(Muons->GetBinContent(i)==0)
    else if(Muons->GetBinContent(i*5-1)==0){
      h->SetBinContent(i,Laser->GetBinContent(i*5)-f->Eval(i*5-1));
      h->SetBinError(i,4.3e-6);
      //h->SetBinError(i,Laser->GetBinError(i*5));
    }
    //if(i%5!=0) h->SetBinContent(i,0);
  }
  //  for (int i=200;i<700;i++) h->SetBinContent(i,Muons->GetBinContent(i));

  //for (int i=0;i<700;i++) h->SetBinError(i,Muons->GetBinError(i));
  //h->Draw();
  
  //Laser->Add(Muons,-1);
  h->SetMarkerStyle(kFullStar);
  h->SetMarkerColor(kRed);
  h->SetMarkerSize(1);
  h->GetXaxis()->SetTitle("Time (#mus)");
  h->GetYaxis()->SetTitle("#DeltaG");
  h->GetYaxis()->SetRangeUser(-0.2e-3,0.05e-3);
  h->GetXaxis()->SetRangeUser(0,150);
  h->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  //Laser->Draw();
  h->Draw();
}
