//THIS CODE DOES NOT WORK but has the framework
{

  gStyle->SetOptFit(); // display fit parameters
  
  // GLOBALS
  const double cyclotron_period = 149.; //ns from Cenap Note #385
  
  // Setup fit parameters ~ basic guesses for general wiggle fit
  const double lifetime = 64400.; // ns
  const double Asymmetry = 0.4;
  const double R = 0;
  const double omega_ppm = 0.001439e-06; // over precise version from geant
  const double omega_fit = 0.001439; // 
  const double phase = 1.5707; // pi/2
  const double blindbyppm = 0;
  
  const double A_c = Asymmetry * cos(phase);
  const double A_s = Asymmetry * sin(phase);
  
  TH1D *ph1_T = histo;
  
  const int binsT = ph1_T->GetNbinsX();
  //const int binsT = 100000;
  
  const float fit_start_T = 0;
  const float fit_stop_T = 596000;
  
  
  TH1D *ph1_T = new TH1D("h1_time_T","",binsT,fit_start_T,fit_stop_T);
 
  TH1D * ph1_r_T = new TH1D("h1_time_residual_T","",binsT,fit_start_T,fit_stop_T);
  TH1D * ph1_nr_T = new TH1D("h1_time_norm_residual_T","",binsT,fit_start_T,fit_stop_T);
  

  // Establish Fit functions
  TF1 * pf_T = new TF1("f_T",Form("[0]*exp(-x/[1])*(1.+[3]*cos(  %.15f*(1.-([2]+%f)/1000000.)*x) - [4]*sin(  %.15f*(1.-([2]+%f)/1000000.)*x )) ",omega_fit, blindbyppm, omega_fit, blindbyppm ),fit_start_T,fit_stop_T);
  pf_T->SetParName(0,"N_{0}");
  pf_T->SetParName(1,"#tau");
  pf_T->SetParName(2,"R");
  pf_T->SetParName(3,"A_c");
  pf_T->SetParName(4,"A_s");
  pf_T->SetParameter(0,ph1_T->GetMaximum()); 
  pf_T->SetParameter(1,lifetime);
  pf_T->SetParameter(2,R);
  pf_T->SetParameter(3,A_c);
  pf_T->SetParameter(4,A_s);
  pf_T->SetNpx(1000);pf_T->SetLineColor(2);
  
  // Do some Fits
  TCanvas* c0 = new TCanvas("c0","",0,0,900,600);
  TH1D *ph1_Tcopy = (TH1D*)ph1_T->Clone();
  ph1_Tcopy->Fit("f_T","R");
  ph1_Tcopy->SetAxisRange(0,200000,"x");
  ph1_Tcopy->Draw();

  double fit_A_c = pf_T->GetParameter(3);
  double fit_A_s = pf_T->GetParameter(4);
  double fit_asymmetry = sqrt(fit_A_c*fit_A_c + fit_A_s*fit_A_s); // always positive
  double fit_phase = atan(fit_A_s / fit_A_c);
  if (fit_A_c < 0) fit_phase += 3.14159; // add pi to reach other quadrant

  cout << "asymmetry = " << fit_asymmetry << endl;
  cout << "phase = " << fit_phase << endl;
  
  
  TCanvas* c1 = new TCanvas("c1","",0,0,900,600);
  
  // Compute Residuals
  for ( int i=1;i<binsT;i++){
    ph1_r_T->SetBinContent(i,ph1_T->GetBinContent(i)-pf_T->Eval(ph1_r_T->GetBinCenter(i)));
  }
  
  // Compute Studentized Residuals
  for ( int i=1;i<binsT;i++){
    if(ph1_r_T->GetBinContent(i) / ph1_T->GetBinError(i)>-100. &&ph1_r_T->GetBinContent(i) / ph1_T->GetBinError(i)<100. ){    ph1_nr_T->SetBinContent(i, ph1_r_T->GetBinContent(i) / ph1_T->GetBinError(i));}
  }


  // Establish Variable Fit Function
  TF1 * pf_Tv = new TF1("f_Tv",Form("[0]*exp(-x/[1])*(1.+[3]*cos(  %.15f*(1.-([2]+%f)/1000000.)*x) - [4]*sin(  %.15f*(1.-([2]+%f)/1000000.)*x )) ",omega_fit, blindbyppm, omega_fit, blindbyppm ));
  pf_Tv->SetParName(0,"N_{0}");
  pf_Tv->SetParName(1,"#tau");
  pf_Tv->SetParName(2,"R");
  pf_Tv->SetParName(3,"A_c");
  pf_Tv->SetParName(4,"A_s");
  pf_Tv->SetNpx(1000);pf_Tv->SetLineColor(2);

  // Establish Histogram to Hold Fit Varieties
  const int steps = 100;

  const float search_around_start = 50000;
  const float start_step_size = search_around_start/steps*2.;
  const float best_stop_time = 650000;

const int stepsInns = 440;
const float start_offset = 0;	  
  
  TH1F * ph1_start_T =   new TH1F("h1_start_T",  ";fit start time (ns);#chi^{2}/ndf",steps,start_offset,start_offset+steps*stepsInns);
  TH1F * ph1_start_tau = new TH1F("h1_start_tau",";fit start time (ns);#tau"        ,steps,start_offset,start_offset+steps*stepsInns);
  TH1F * ph1_start_rrr = new TH1F("h1_start_rrr",";fit start time (ns);R"           ,steps,start_offset,start_offset+steps*stepsInns);
  TH1F * ph1_start_Ac =  new TH1F("h1_start_Ac", ";fit start time (ns);A_c"         ,steps,start_offset,start_offset+steps*stepsInns);
  TH1F * ph1_start_As =  new TH1F("h1_start_As", ";fit start time (ns);A_s"         ,steps,start_offset,start_offset+steps*stepsInns);

  for (int i=1;i<=steps;i++){

    pf_Tv->SetParameter(0,ph1_T->GetBinContent(ph1_T->FindBin(search_around_start-steps*start_step_size/2.+i*start_step_size)));
    pf_Tv->SetParameter(1,lifetime);
    pf_Tv->SetParameter(2,R);
    pf_Tv->SetParameter(3,A_c);
    pf_Tv->SetParameter(4,A_s);
    ph1_T->Fit("f_Tv","Q","",start_offset+i*stepsInns,best_stop_time);
    cout << i << " " << start_offset+i*stepsInns << "  chi2/ndf = " << pf_Tv->GetChisquare()/pf_Tv->GetNDF()  << "  vs  1 +/- " << sqrt(2./pf_Tv->GetNDF()) << endl;

    double fit_A_c = pf_Tv->GetParameter(3);
    double fit_A_s = pf_Tv->GetParameter(4);
    double fit_asymmetry = sqrt(fit_A_c*fit_A_c + fit_A_s*fit_A_s); // always positive
    double fit_phase = atan(fit_A_s / fit_A_c);
    if (fit_A_c < 0) fit_phase += 3.14159; // add pi to reach other quadrant

    ph1_start_T->SetBinContent(i,pf_Tv->GetChisquare()/pf_Tv->GetNDF());    ph1_start_T->SetBinError(i,sqrt(2./pf_Tv->GetNDF()));
    ph1_start_tau->SetBinContent(i,pf_Tv->GetParameter(1));                    ph1_start_tau->SetBinError(i,pf_Tv->GetParError(1));
    ph1_start_rrr->SetBinContent(i,pf_Tv->GetParameter(2));                    ph1_start_rrr->SetBinError(i,pf_Tv->GetParError(2));
    ph1_start_Ac->SetBinContent(i,pf_Tv->GetParameter(3));                    ph1_start_Ac->SetBinError(i,pf_Tv->GetParError(3));
    ph1_start_As->SetBinContent(i,pf_Tv->GetParameter(4));                    ph1_start_As->SetBinError(i,pf_Tv->GetParError(4));

  }

/*
  const float stop_step_size = cyclotron_period*30.;
  const float search_around_stop = 325000;
  const float best_start_time = 0;

  TH1F * ph1_stop_T = new TH1F("h1_stop_T",";fit stop time (ns);#chi^{2}/ndf",steps,search_around_stop-steps*stop_step_size/2.,search_around_stop+steps*stop_step_size/2.);
  TH1F * ph1_stop_tau = new TH1F("h1_stop_tau",";fit stop time (ns);#tau",steps,search_around_stop-steps*stop_step_size/2.,search_around_stop+steps*stop_step_size/2.);
  TH1F * ph1_stop_asy = new TH1F("h1_stop_asy",";fit stop time (ns);A"   ,steps,search_around_stop-steps*stop_step_size/2.,search_around_stop+steps*stop_step_size/2.);
  TH1F * ph1_stop_rrr = new TH1F("h1_stop_rrr",";fit stop time (ns);R"   ,steps,search_around_stop-steps*stop_step_size/2.,search_around_stop+steps*stop_step_size/2.);
  TH1F * ph1_stop_phi = new TH1F("h1_stop_phi",";fit stop time (ns);#phi",steps,search_around_stop-steps*stop_step_size/2.,search_around_stop+steps*stop_step_size/2.);

  for (int i=1;i<=steps;i++){
    pf_Tv->SetParameter(0,ph1_T->GetBinContent(ph1_T->FindBin(best_start_time)));
    pf_Tv->SetParameter(1,lifetime);
    pf_Tv->SetParameter(2,Asymmetry);
    pf_Tv->SetParameter(3,R);
    pf_Tv->SetParameter(4,phase);
    ph1_T->Fit("f_Tv","Q","",best_start_time,search_around_stop-steps*stop_step_size/2.+i*stop_step_size);

    cout << i << " " << search_around_stop-steps*stop_step_size/2.+i*stop_step_size << "  chi2/ndf = " << pf_Tv->GetChisquare()/pf_Tv->GetNDF()  << "  vs  1 +/- " << sqrt(2./pf_Tv->GetNDF()) << endl;
    ph1_stop_T->SetBinContent(i,pf_Tv->GetChisquare()/pf_Tv->GetNDF());     ph1_stop_T->SetBinError(i,sqrt(2./pf_Tv->GetNDF()));
    ph1_stop_tau->SetBinContent(i,pf_Tv->GetParameter(1));                    ph1_stop_tau->SetBinError(i,pf_Tv->GetParError(1));
    ph1_stop_asy->SetBinContent(i,pf_Tv->GetParameter(2));                    ph1_stop_asy->SetBinError(i,pf_Tv->GetParError(2));
    ph1_stop_rrr->SetBinContent(i,pf_Tv->GetParameter(3));                    ph1_stop_rrr->SetBinError(i,pf_Tv->GetParError(3));
    ph1_stop_phi->SetBinContent(i,pf_Tv->GetParameter(4));                    ph1_stop_phi->SetBinError(i,pf_Tv->GetParError(4));
  }
*/

  TCanvas * pcan03 = new TCanvas(); pcan03->Divide(1,2);
  pcan03->cd(1); ph1_start_T->Draw();
  pcan03->cd(2); ph1_T->Draw();
//  pcan03->cd(2); ph1_stop_T->Draw();

/*
TCanvas * pcan04 = new TCanvas(); pcan04->Divide(1,4);
  pcan04->cd(1); ph1_stop_tau->Draw();
  pcan04->cd(2); ph1_stop_asy->Draw();
  pcan04->cd(3); ph1_stop_rrr->Draw();
  pcan04->cd(4); ph1_stop_phi->Draw();
*/

  TCanvas * pcan05 = new TCanvas("pcan05","",0,0,900,900); pcan05->Divide(1,4);
  pcan05->cd(1); ph1_start_tau->Draw();
  pcan05->cd(2); ph1_start_rrr->Draw();
  pcan05->cd(3); ph1_start_Ac->Draw();
  pcan05->cd(4); ph1_start_As->Draw();

}
