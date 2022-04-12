{
  /*TFile *_file0 = TFile::Open("osc_18_c500_80.root");
  Muons->Draw();
  printf("Bin error %f\n",Muons->GetBinError(30));
  printf("Bin content %f\n",Muons->GetBinContent(30));
  //printf("Bin value %f\n",Muons->GetBinContent(30));
  Muons->GetXaxis()->SetRange(29,30);
  printf("Bin mean content x axis %f\n",Muons->GetMean(1));
  printf("Bin mean content y axis %f\n",Muons->GetMean(2));
  */

  TF1* f2 = new TF1("f2","[0]*(1+[2]*cos(1.439*(1-[4]*1.e-6)*x+([3]+2*pi)))*exp(-x/[1])",0.0,700);
  f2->SetParameters(2e11,64.4, 0.4, 4.5e-6, 0);
  f2->SetNpx(7000);
  printf("Function value %f\n",f2->Eval(0));
}
