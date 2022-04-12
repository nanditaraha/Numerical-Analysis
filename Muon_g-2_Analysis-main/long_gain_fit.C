/* Few changes
Fill time = 700 us..
Took one fill looped over 20 for 750 data pts. keeping an event fixed.
Repeat the above for 2000 events... 
N0 is for a single fill for our function, but fit function will have 20*2000 a fill 
*/
#include <TH1D.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include "TApplication.h"
#include <TRandom3.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <fstream>

#include<math.h>
#include<TAttParticle.h>


void wiggle(int events=2000, int fills=16) {
  //Variables
  Double_t pi = 3.14159265358979323846;
  Double_t Emax = 3.1;
  Int_t bin = 100;
  // Int_t bin = 1000;
  Double_t eps = 0.001;//Gain change of SiPM for one fill ;
  //Double_t eps = 0.05;//for one hour;
  Double_t eps_pmt = 0.05/(events*fills);//for one hour;
  //for one hour we take eps as 5%
  //Double_t eps_pmt = 0.05*700.0/3.6e8; //PMTs vary by 5% in an hour (laser monitoring PMT's)
  Double_t prob;
  Double_t y;
  Double_t pert_y;
  Double_t mean;
  Double_t E;
  Double_t var;
  //Double_t Nt;
  Double_t fNt;
  Double_t gain = 0.;
  double gain_pmt = 0;
  double gain_ratio = 0;
  Double_t gain_array[3000][5000]={0};
  Double_t turn = 149; //considering 1 turn
  //   Double_t turn = 1; //considering every ns 
  //Double_t Norm = 2.e11; //Number in Total data set 
  //Double_t Norm = 5.e8; //N0 for total data set
  //Double_t Norm = 1.15e4; //single fill
  Double_t Norm = 1.15e4; //single fill*2000  - our req for testing
  Double_t tau = 64400.;
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  //Double_t fphase = pi/2;
  Double_t fphase = 4.5e-6;
  Double_t time = 4700; //700 mus if turn 149 ns
  //Double_t time = 4000; //600 mus if turn 149 ns
  //   Double_t time = 1000; //150 mus if turn 149 ns
  //   Double_t Ntot = 0;
  TRandom3 rd;
  //parameters for the phase function
  double p0 =    -0.255134;
  double p1 =      65.3034;
  double p2 =     -705.492;
  double p3 =      5267.21;
  double p4 =     -23986.5;
  double p5 =      68348.1;
  double p6 =      -121761;
  double p7 =       131393;
  double p8 =       -78343;
  double p9 =      19774.1;



  //N function
  //TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set
  TF1* Nfunc = new TF1("Nfunc","[0]*[1]*1.15e4*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);//single fill*20*events
  //   Nfunc->SetParNames("Norm");
  //   Nfunc->SetParameters(Norm);
  Nfunc->SetNpx(10000000);
  Nfunc->SetLineColor(1);
  Nfunc->SetParameters(events,fills);
  //  Nfunc->SetParameter(0,events);
  
  TH1D* hn = new TH1D("hn","N(y)",bin,0.,1.); 
  hn->Eval(Nfunc);
  
  //phase function
  TF1* phase = new TF1("phase","[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5 + [6]*x^6 + [7]*x^7 + [8]*x^8 + [9]*x^9",0.,1.);
  phase->SetParNames("p0","p1","p2","p3","p4","p5","p6","p7","p8","p9");
  phase->SetParameters(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);
  phase->SetNpx(10000000);
  phase->SetLineColor(2);
  TH1D* hp = new TH1D("hp","Phi(y)",bin,0.,1.);
  hp->Eval(phase);
  
  //Asymmetri function
  TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.);
  Asy->SetNpx(10000000);
  Asy->SetLineColor(3);
  ha->Eval(Asy);
  
  //Fit Function
  TF1* fit_func1 = new TF1("fit_func1","[0]*(1+[2]*cos(0.001439*(1-[4]*1.e-6)*x+([3]+2*pi)))*exp(-x/[1])",0.0,time*turn);
  fit_func1->SetParNames("Norm","tau","Asymmetry","Phase","R(ppm)");
  fit_func1->SetParameters(Norm*events*fills,tau, fAsy, fphase, R);
  //   fit_func1->SetParLimits(fphase,0., 2*pi);
  fit_func1->SetNpx(time);
  fit_func1->SetLineColor(2);
  
  //wiggle plot
  TF1* wiggle = new TF1("wiggle","[0]*(1+[1]*cos(0.001439*x+[2]))*exp(-x/64400)",0.,time*turn);
  wiggle->SetParNames("N","A","phi");
  wiggle->SetParameters(Norm*events*fills, fAsy, fphase);
  //wiggle->SetNpx(10000000);
  
  TH1D* h_wiggletot = new TH1D("h_wiggletot","wiggle(t)",time, 0., time*turn);
    
  //gain function
  TH1D* h_gain = new TH1D("h_gain", "Gain(t)", time, 0., time*turn);
  TH2D* h_gain_ratio = new TH2D("h_gain_ratio", "Gain(t) ratio", time, 0., time*turn, events, 0., events);
  TH1D* h_gain_pmt = new TH1D("h_gain_pmt", "Gain(t) PMT", events, 0., events);
    
  //Loop for evaluating SiPM gains for each /all fills which is same
  for(int it=0; it<time;it++){
    double laser_offset_t = 5e3/turn;//corresponds to 5 us in bin size for each laser firing point
    gain = 1. + eps*TMath::Exp((-(it- laser_offset_t)*turn)/tau); //exponential gain fluctutation
    h_gain->SetBinContent(it,gain);
  }  

  //Loop for events or injector cycle                                                                                                                  
  for(int i=0;i<events;i++){
    gain_pmt=1+eps_pmt*i;  
    for(int it=0; it<time;it++){
      Double_t Ntot=0;
      Double_t Nt=0;
      Double_t del = rd.Gaus(0.,1.);
      gain_ratio=h_gain->GetBinContent(it)/gain_pmt; 

      for(int i=60;i<bin;i++){
	Double_t Nytot=0;
	for(int j=0; j<bin; j++){
	  y = float(j)/bin;
	  y = y/(gain_ratio);                                                                                                                                        
	  mean = y+(1./(2*bin));
	  E = mean*Emax;
	  var = 0.05*sqrt(E);
	  //var = 1e-6;                                                                                                                                          
	  Double_t probt=0;
	  Double_t Na = hn->GetBinContent(j);
	  Double_t A = ha->GetBinContent(j);
	  Double_t phi = hp->GetBinContent(j);
	  prob = (-TMath::Erfc((((float(i+1))/bin)-mean)/var)+TMath::Erfc(((float(i)/bin)-mean)/var))/(2.);
	  probt= probt+prob;
	  //p_matrix->Fill(i,j, prob);
	  //h_prob->SetBinContent(j,prob);
	  Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
	  Nytot = Nytot+Ny;
	} //end for j                                                                                                                                            
	Ntot = Ntot+Nytot;
      }//end for i                                                                                                                                               
      Nt = Ntot*TMath::Exp(-(it*turn)/tau);
      Nt = (Nt + del*sqrt(Nt));
      h_wiggletot->SetBinContent(it,Nt);
      //h_gain->SetBinContent(it, gain);
      
      
    }//end for it
  }





  TCanvas* c2 = new TCanvas("c2");
  //c2->Divide(1,2);
  h_wiggletot->Draw();
  //c2->cd(1);
  //h_gain->Draw();
  //c2->cd(2);
  //h_gain_pmt->Draw();
  //TCanvas* c1 = new TCanvas("c1");
  //h_gain_ratio->Draw("colz");
  //fit_func1->Draw("same");
  //fit_func1->Draw();
  h_wiggletot->Fit("fit_func1","R");
  h_wiggletot->SetXTitle("t[ns]");
  h_wiggletot->SetYTitle("Entries");
  gStyle->SetOptFit(1);

  TFile f("exp_wiggle.root","new");
  h_wiggletot->Write();
  
}


void stable_wiggle(int events=2700, int fills=16)
{

  //Variables                                                                                                                                          
  Double_t pi = 3.14159265358979323846;
  Double_t Emax = 3.1;
  Int_t bin = 100;
    Double_t prob;
  Double_t y;
  Double_t pert_y;
  Double_t mean;
  Double_t E;
  Double_t var;
  //Double_t Nt;                                                                                                                                       
  Double_t fNt;
  Double_t turn = 149; //considering 1 turn                                                                                                            
  Double_t Norm = 1.15e4; //single fill*2000  - our req for testing                                                                                    
  Double_t tau = 64400.;
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  Double_t fphase = 4.5e-6;
  Double_t time = 4700; //700 mus if turn 149 ns                                                                                                       
  TRandom3 rd;
  //parameters for the phase function                                                                                                                  
  double p0 =    -0.255134;
  double p1 =      65.3034;
  double p2 =     -705.492;
  double p3 =      5267.21;
  double p4 =     -23986.5;
  double p5 =      68348.1;
  double p6 =      -121761;
  double p7 =       131393;
  double p8 =       -78343;
  double p9 =      19774.1;

  //N function                                                                                                                                         
  //TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set                                         
  TF1* Nfunc = new TF1("Nfunc","[0]*[1]*1.15e4*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);//single fill*20*events                       
  Nfunc->SetNpx(10000000);
  Nfunc->SetLineColor(1);
  Nfunc->SetParameters(events,fills);

  TH1D* hn = new TH1D("hn","N(y)",bin,0.,1.);
  hn->Eval(Nfunc);

  //phase function                                                                                                                                     
  TF1* phase = new TF1("phase","[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5 + [6]*x^6 + [7]*x^7 + [8]*x^8 + [9]*x^9",0.,1.);
  phase->SetParNames("p0","p1","p2","p3","p4","p5","p6","p7","p8","p9");
  phase->SetParameters(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);
  phase->SetNpx(10000000);
  phase->SetLineColor(2);
  TH1D* hp = new TH1D("hp","Phi(y)",bin,0.,1.);
  hp->Eval(phase);

  //Asymmetri function                                                                                                                                 
  TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.);
  Asy->SetNpx(10000000);
  Asy->SetLineColor(3);
  ha->Eval(Asy);

  //Fit Function  
  TF1* fit_func1 = new TF1("fit_func1","[0]*(1+[2]*cos(0.001439*(1-[4]*1.e-6)*x+([3]+2*pi)))*exp(-x/[1])",0.0,time*turn);
  fit_func1->SetParNames("Norm","tau","Asymmetry","Phase","R(ppm)");
  fit_func1->SetParameters(Norm*events*fills,tau, fAsy, fphase, R);
  fit_func1->SetNpx(time);
  fit_func1->SetLineColor(2);

  //wiggle plot  
  TF1* wiggle = new TF1("wiggle","[0]*(1+[1]*cos(0.001439*x+[2]))*exp(-x/64400)",0.,time*turn);
  wiggle->SetParNames("N","A","phi");
  wiggle->SetParameters(Norm*events*fills, fAsy, fphase);
  //wiggle->SetNpx(10000000);                               

  TH1D* h_wiggletot_stable = new TH1D("h_wiggletot_stable","wiggle(t)",time, 0., time*turn);
  for(int it=0; it<time;it++){
    Double_t Ntot=0;
    Double_t Nt=0;
    Double_t del = rd.Gaus(0.,1.);
    
    for(int i=60;i<bin;i++){
      Double_t Nytot=0;
      for(int j=0; j<bin; j++){
	y = float(j)/bin;
	mean = y+(1./(2*bin));
	E = mean*Emax;
	var = 0.05*sqrt(E);
	
	Double_t probt=0;
	Double_t Na = hn->GetBinContent(j);
	Double_t A = ha->GetBinContent(j);
	Double_t phi = hp->GetBinContent(j);
	prob = (-TMath::Erfc((((float(i+1))/bin)-mean)/var)+TMath::Erfc(((float(i)/bin)-mean)/var))/(2.);
	probt= probt+prob;
	Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
	Nytot = Nytot+Ny;
      } //end for j

      Ntot = Ntot+Nytot;
    }//end for i           
    Nt = Ntot*TMath::Exp(-(it*turn)/tau);
    Nt = (Nt + del*sqrt(Nt));
    h_wiggletot_stable->SetBinContent(it,Nt);
  }

  TCanvas* c2 = new TCanvas("c2");
  h_wiggletot_stable->Draw();
  h_wiggletot_stable->Fit("fit_func1","R");
  h_wiggletot_stable->SetXTitle("t[ns]");
  h_wiggletot_stable->SetYTitle("Entries");
  gStyle->SetOptFit(1);

  TFile f("stable_wiggle.root","new");
  h_wiggletot_stable->Write();
}
