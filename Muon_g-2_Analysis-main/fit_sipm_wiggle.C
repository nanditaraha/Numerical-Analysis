/*
Fitting the gain function.
  */

void fit_sipm_wiggle() {
  TRandom3 rd;
  //TFile *f = TFile::Open("wiggle_fac10_g320_m100.root"); //This file includes wiggles in the plot too!
  TFile *f = TFile::Open("noWiggle_fac10_g320_m100.root");
  //TH1D *h1=(TH1D*)f->Get("Laser");// As this file has no lasers - just the basic fit
  TH1D *h1=(TH1D*)f->Get("Muons");// As this file has no lasers - just the basic fit
  h1->SetLineColor(kBlue);
  //h1->GetXaxis()->SetRangeUser(0,700);
   
  for(int i=0; i<700; i++){
   if(h1->GetBinContent(i+1)!=0)
     {
       h1->SetBinError(i+1,0.02/sqrt(2000));
       
       
       h1->SetBinContent((i+1),h1->GetBinContent(i+1)+(0.02/sqrt(2000))*rd.Gaus(0.,1.));
       //h1->SetBinError(i+1,sqrt(0.02*0.02/2000+h1->GetBinError(i+1)*h1->GetBinError(i+1)));
       //Error set as Quadrature of a fixed error due to 2000 pts and 2% rms and error of the plot itself
       
       //h1->SetBinError(i+1,sqrt(0.02/2000));// this works best but all points then have a constant error bar ans so not very convincing...
     }
   else h1->SetBinError(i+1,0.0);
  }
  /*
  TF1 *S_N= new TF1("S_N","24.9828*(exp(-x/64.4)-exp(-x/18.))", 0,700);
  TF1 *S_a= new TF1("S_a","0.000416279*((64.4-18)*exp(-x/64.4)*cos(1.449*x)-exp(-x/18.) + 64.4*18*1.439*exp(-x/64.4)*sin(1.449*x))", 0,700);
  TF1 *S_b= new TF1("S_b","0.000416279*((64.4-18)*exp(-x/64.4)*sin(1.449*x)-exp(-x/18.) + 64.4*18*1.439*exp(-x/64.4)*cos(1.449*x))", 0,700);
  
  TF1 *f2 = new TF1("f2","1-100*[0]*([1]*S_N+[2]*S_a+[3]*S_b)", 0,700);
  f2->SetParameters(2.38344e-05,0.420282,0.0457639,-0.000865762);
  f2->SetNpx(1000);
  */

  TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*((1+[3]*cos(1.449*x))*exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,400);
  //TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*(1+[3]*cos(1.449*x)*exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,400);
  f3->SetParNames("P_{0}","n_{0}","#tau_{r}","A");
  //TF1 *f3 = new TF1("f3","1-100*[0]*([1]*S_N+[2]*S_a+[3]*S_b)", 0,700);
  f3->SetParameters(1.052e-5,100,18,0.006);
  //f3->SetParameters(2.38344e-05,0.420282,0.0457639,-0.000865762);
  f3->SetNpx(1000);
  f3->SetLineColor(kBlue);
    
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  h1->SetTitle("");
  h1->GetYaxis()->SetRangeUser(0.987,1.001);
  h1->GetXaxis()->SetRangeUser(0,250);
  h1->GetXaxis()->SetTitle("Time(#mus)");
  h1->GetYaxis()->SetTitle("#frac{G}{G_{0}}");
  h1->SetStats(1);
  h1->Fit("f3","R");
  
  TLatex latex;
  latex.SetTextSize(0.037);
  latex.SetTextColor(kBlue);
  latex.DrawLatex(121.5,.993,"#frac{G}{G_{0}} = 1 - #frac{p_{0} n_{0} #tau_{r}}{#tau_{#mu}-#tau_{r}}#left([1+A_{g} cos(#omega_{a} t + #phi_{g})]e^{-t/#tau_{#mu}} - e^{-t/#tau_{r}} #right)");
  gPad->SetGridx();
  gPad->SetGridy();
  
}

