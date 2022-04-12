
//Try to simulate and normalise different types of functions of any form with simulated data...

{
  TF1 *f1 = new TF1("f1","2*x*x + x +5",0,10);
  TH1F h2("h2","",10000,0,10);
  for (Int_t i=0;i<10000000;i++)
    {
      r = f1->GetRandom();
      h2.Fill(r,f1->Integral(0,10)*1000/10000000);
      //for (Int_t j=0;j<10;j++)
      //h2.Fill(j,r);
    }
  h2.Draw();
  f1->Draw("same");

}

