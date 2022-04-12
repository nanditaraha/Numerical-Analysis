double gain(double *y, double *p){
  double x = y[0];
  double n = 0.420484;
  double a = 0.0457639;
  double b = -0.000865762;
  double w = 1.439; //precession freq 
  double nk = n*(exp(-(x/64.4))- exp(-(x/p[2])))/(64.4-p[2]);
  double ak =  a*((64.4-p[2])*(cos(w*x)*exp(-(x/64.4))- exp(-(x/p[2])))+w*64.4*p[2]*(sin(w*x)*exp(-(x/64.4))))/((1+w*64.4*w*64.4)*p[2]*p[2]-2*64.4*p[2]+64.4*64.4);
  double bk = b*((64.4-p[2])*sin(1.439*x)*exp(-(x/64.4)) - w*64.4*p[2]*(cos(w*x)*exp(-(x/64.4)) -exp(-(x/p[2]))))/((1+w*64.4*w*64.4)*p[2]*p[2]-2*64.4*p[2]+64.4*64.4);

  return 1-p[0]*p[1]*64.4*p[2]*(nk + ak + bk);
}

void gain_formula(double tau_r=18)
{
  TF1 *g1 = new TF1("Gain", gain,0,200,3);
  g1->SetParameters(100,1e-5,20);
  g1->SetNpx(1000);
  g1->SetLineColor(kBlack);
  g1->GetXaxis()->SetTitle("Time (#mus)");
  g1->GetYaxis()->SetTitle("#frac{G}{G_{0}}");
  g1->Draw();
  
  TF1 *g = new TF1("g", gain,0,200,3);
  g->SetParameters(100,1e-5,10);
  g->SetNpx(1000);
  g->Draw("same");
  
  TF1 *g2 = new TF1("g2", gain,0,200,3);
  g2->SetParameters(100,1e-5,5);
  g2->SetNpx(1000);
  g2->SetLineColor(kBlue);
  g2->Draw("same");

  TF1 *g3 = new TF1("g3", gain,0,200,3);
  g3->SetParameters(100,1e-5,2);
  g3->SetNpx(1000);
  g3->SetLineColor(kGreen);
  g3->Draw("same");


   TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
   //leg->SetHeader("Recovery times","C"); // option "C" allows to center the header
   leg->AddEntry(g1,"#tau_{r} = 20 #mus","l");
   leg->AddEntry("g","#tau_{r} = 10 #mus","l");
   leg->AddEntry("g2","#tau_{r} = 5 #mus","l");
   leg->AddEntry("g3","#tau_{r} = 2 #mus","l");
   leg->Draw();


  
  
}
