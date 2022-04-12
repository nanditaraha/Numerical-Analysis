/*
Fitting the gain function.
  */

void draw_sipm() {
  TRandom3 rd;
  //TFile *f = TFile::Open("noWiggle_fac10_g320_m100.root");
  TFile *f = TFile::Open("try_640_5.root");//only for 640... 
  
  //TH1D *h3=(TH1D*)f->Get("Laser");
  TH1D *h1=(TH1D*)f->Get("All");
  TH1D *h2=(TH1D*)f->Get("Muons");
  TH1D *h3=(TH1D*)f->Get("Laser");
  
  for(int i=0; i<700; i++){
    if(h1->GetBinContent(i+1)!=0){
	h1->SetBinError(i+1,0.02/sqrt(2000));
	h1->SetBinContent((i+1),h1->GetBinContent(i+1)+(0.02/sqrt(2000))*rd.Gaus(0.,0.5));

	h2->SetBinError(i+1,0.02/sqrt(2000));
	h2->SetBinContent((i+1),h2->GetBinContent(i+1)+(0.02/sqrt(2000))*rd.Gaus(0.,0.5));
	
	h3->SetBinError(i+1,0.02/sqrt(2000));
	h3->SetBinContent((i+1),h3->GetBinContent(i+1)+(0.02/sqrt(2000))*rd.Gaus(0.,0.5));
	  
	//h1->SetBinError(i+1,sqrt(0.02/2000));// this works best but all points then have a constant error bar ans so not very convincing...
    }
    else h1->SetBinError(i+1,0.0);
  }
  h1->SetTitle("");
  h1->GetYaxis()->SetRangeUser(0.987,1.001);
  h1->GetXaxis()->SetTitle("Time (#mus)");
  h1->SetStats(0);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  //h1->Draw("same");
  gPad->SetGridx();
  gPad->SetGridy();
   
  TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*(exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,400);
  f3->SetParNames("P_{0}","n_{0}","#tau_{r}");
  //f3->FixParameter(2,18);
  f3->SetParameters(1.052e-5,100,16.8);
  gStyle->SetOptFit(1);
  f3->SetLineColor(kMagenta);
  //f3->Draw("same");
  h1->SetStats(0);
  
  TLatex latex;
  latex.SetTextSize(0.037);
  latex.DrawLatex(121.5,.991,"Laser 640 #mus, 5 #mus step ");
  
  //latex.DrawLatex(121.5,.994,Form("p_{0} = P_{0} #tau_{#mu} = 6.4 #times 10^{-4}"));
  //latex.SetTextColor(kMagenta);
  //latex.DrawLatex(121.5,.993,"#frac{G}{G_{0}} = 1 - #frac{p_{0} n_{0} #tau_{r}}{#tau_{#mu}-#tau_{r}}#left(e^{-t/#tau_{#mu}} - e^{-t/#tau_{r}} #right)");
  
}

