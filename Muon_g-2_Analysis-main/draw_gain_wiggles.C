/*
Fitting the gain function.
  */

void draw_gain_wiggles() {

  
  
  TF1 *S_N1= new TF1("S_N1","64.4*1.*(exp(-x/64.4)-exp(-x/1.))/(64.4-1.)", 0,200);
  TF1 *S_a1= new TF1("S_a1","64.4*1.*((64.4-1.)*exp(-x/64.4)*cos(1.449*x)-exp(-x/1.) + 64.4*1.*1.439*exp(-x/64.4)*sin(1.449*x))/((1+1.449*1.449*64.4*64.4)*1.*1.-2*64.4*1.+64.4*64.4)", 0,200);
  TF1 *S_b1= new TF1("S_b1","64.4*1.*((64.4-1.)*exp(-x/64.4)*sin(1.449*x) - 64.4*1.*1.439*(exp(-x/64.4)*cos(1.449*x)-exp(-x/1.))) /((1+1.449*1.449*64.4*64.4)*1.*1.-2*64.4*1.+64.4*64.4)", 0,200);
  
  TF1 *f21 = new TF1("f21","1-100*[0]*([1]*S_N1+[2]*S_a1+[3]*S_b1)", 0,200);
  f21->SetParameters(1.1*2.38344e-04,0.420282,0.0457639,-0.000865762);
  f21->SetNpx(1000);


TF1 *S_N2= new TF1("S_N2","64.4*5.*(exp(-x/64.4)-exp(-x/5.))/(64.4-5.)", 0,200);
  TF1 *S_a2= new TF1("S_a2","64.4*5.*((64.4-5.)*exp(-x/64.4)*cos(1.449*x)-exp(-x/5.) + 64.4*5.*1.439*exp(-x/64.4)*sin(1.449*x))/((1+1.449*1.449*64.4*64.4)*5.*5.-2*64.4*5.+64.4*64.4)", 0,200);
  TF1 *S_b2= new TF1("S_b2","64.4*5.*((64.4-5.)*exp(-x/64.4)*sin(1.449*x) - 64.4*5.*1.439*(exp(-x/64.4)*cos(1.449*x)-exp(-x/5.))) /((1+1.449*1.449*64.4*64.4)*5.*5.-2*64.4*5.+64.4*64.4)", 0,200);

  
  TF1 *f22 = new TF1("f22","1-100*[0]*([1]*S_N2+[2]*S_a2+[3]*S_b2)", 0,200);
  f22->SetParameters(2.3*2.38344e-05,0.420282,0.0457639,-0.000865762);
  f22->SetNpx(1000);

TF1 *S_N3= new TF1("S_N3","64.4*10.*(exp(-x/64.4)-exp(-x/10.))/(64.4-10.)", 0,200);
  TF1 *S_a3= new TF1("S_a3","64.4*10.*((64.4-10.)*exp(-x/64.4)*cos(1.449*x)-exp(-x/10.) + 64.4*10.*1.439*exp(-x/64.4)*sin(1.449*x))/((1+1.449*1.449*64.4*64.4)*10.*10.-2*64.4*10.+64.4*64.4)", 0,200);
  TF1 *S_b3= new TF1("S_b3","64.4*10.*((64.4-10.)*exp(-x/64.4)*sin(1.449*x) - 64.4*10.*1.439*(exp(-x/64.4)*cos(1.449*x)-exp(-x/10.))) /((1+1.449*1.449*64.4*64.4)*10.*10.-2*64.4*10.+64.4*64.4)", 0,200);

  
  TF1 *f23 = new TF1("f23","1-100*[0]*([1]*S_N3+[2]*S_a3+[3]*S_b3)", 0,200);
  f23->SetParameters(1.2*2.38344e-05,0.420282,0.0457639,-0.000865762);
  f23->SetNpx(1000);

  //Used calculated values with tau_r - check file CalibrationProposalIntro.pdf page 10 for formulae
  TF1 *S_N4= new TF1("S_N4","24.9828*(exp(-x/64.4)-exp(-x/18.))", 0,200);
  TF1 *S_a4= new TF1("S_a4","0.000416279*((64.4-18.)*exp(-x/64.4)*cos(1.449*x)-exp(-x/18.) + 64.4*18*1.439*exp(-x/64.4)*sin(1.449*x))", 0,200);
  TF1 *S_b4= new TF1("S_b4","0.000416279*((64.4-18.)*exp(-x/64.4)*sin(1.449*x)-exp(-x/18.) + 64.4*18*1.439*exp(-x/64.4)*cos(1.449*x))", 0,200);
 
  TF1 *f24 = new TF1("f24","1-100*[0]*([1]*S_N4+[2]*S_a4+[3]*S_b4)", 0,200);
  f24->SetParameters(0.64*2.38344e-05,0.420282,0.0457639,-0.000865762);
  f24->SetNpx(1000);

  f21->SetLineColor(kGreen);
  f22->SetLineColor(kBlue);
  //f23->SetLineColor();
  f23->SetLineColor(kBlack);

  f21->Draw();
  f24->Draw("same");
  f22->Draw("same");
  f23->Draw("same");
  //f24->Draw();
 
  TLatex latex;
  latex.SetTextSize(0.037);
  latex.SetTextColor(kRed);
  latex.DrawLatex(121.5,.993,"#tau_{r} = 18 #mus");
  latex.SetTextColor(kBlack);
  latex.DrawLatex(121.5,.992,"#tau_{r} = 10 #mus");
  
  latex.SetTextColor(kBlue);
  latex.DrawLatex(121.5,.991,"#tau_{r} = 5 #mus");
  
  latex.SetTextColor(kGreen);
  latex.DrawLatex(121.5,.990,"#tau_{r} = 1 #mus");
  
}

