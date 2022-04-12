{
    TFile *f = TFile::Open("pulse_laser_1000.root");
    //Laser->SetLineColor(2);
    //Laser->Draw();
    TGraph *Laser = (TGraph*)f->Get("Laser");
    
    TFile *f1 = TFile::Open("sample_only_1000.root");
    //Laser->SetLineColor(kRed);
    //MySample->Draw("same");
    TGraph *Sample = (TGraph*)f->Get("MySample");
    
    double x[1000], y[1000];
    
    for (int i=0; i<1000; i++)
    {
        y[i]=Laser->Eval(i)-MySample->Eval(i);
        x[i]=i;
        
    }
    
    TGraph *gr = new TGraph(1000, x, y);
      //  gr->SetPoint(i,i,Laser->GetY()-Sample->GetY());
        
        gr->Draw("AC");

}
  
  
