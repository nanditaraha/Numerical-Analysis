/**** Original Author Antonio  modified by Nandita
********/
//Call startScan() to run scan and wiggle to just plot one histogram

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

//Global Variables and arrays for fit scan:
Double_t R_w[500]={0};
Double_t R_err[500]={0};//statistical error
Double_t R_errs[500]={0};//allowed set subset error
Double_t errx[500]={0};

Double_t T[500]={0};
Double_t T_err[500]={0};//statistical error 
Double_t T_errs[500]={0};//allowed set subset error 

Double_t phi[500]={0};
Double_t phi_err[500]={0};//statistical error   
Double_t phi_errs[500]={0};//allowed set subset error  

Double_t A[500]={0};
Double_t A_err[500]={0};//statistical error 
Double_t A_errs[500]={0};//allowed set subset error

Double_t chi[500]={0};
Double_t chi_err[500]={0};//statistical error 
Double_t chi_errs[500]={0};//allowed set subset error 

Double_t eps[7] = {-0.003,-0.002,-0.001,0.,0.001,0.002,0.003};

//Fit func def...
double fit_func(double *x, double *par){

  return par[0]*(1+par[2]*cos(0.001439*(1-par[4]*1.e-6)*x[0]+(par[3]+2*TMath::Pi())))*exp(-x[0]/par[1]);
}


void wiggle(int n=7, double start=0.) {
  //Variables
  //Double_t pi = 3.14159265358979323846;
  Double_t Emax = 3.1;
  Int_t bin = 100;
  //Int_t bin = 1000;
  //Double_t eps = 0.001;
  
  Double_t prob, prob1, prob2, prob3;
  Double_t y;
  Double_t pert_y;
  Double_t mean;
  Double_t E;
  Double_t resolution;
  Double_t Nt;
  Double_t fNt;
  Double_t gain = 0.;
  Double_t turn = 149; //considering 1 turn
  //   Double_t turn = 1; //considering every ns 
  Double_t Norm = 2.e11; //Total data set
  //   Double_t Norm = 1.15e4; //single fill
  Double_t tau = 64400.;
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  //Double_t fphase = TMath::Pi()/2;
  Double_t fphase = 4.5e-5;
  //   Double_t time = 4700; //700 mus if turn 149 ns
  Double_t time = 4000; //600 mus if turn 149 ns
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

  //Parameters for asymmetry
  
  double ap0 =    -0.138813;
  double ap1 =      -0.3202;
  double ap2 =     -2.33637;
  double ap3 =      55.7403;
  double ap4 =     -299.495;
  double ap5 =      833.345;
  double ap6 =     -1330.95;
  double ap7 =      1230.64;
  double ap8 =     -612.486;
  double ap9 =      126.984;
  

  //N function
  TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set
  //   TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);//single fill
  //   Nfunc->SetParNames("Norm");
  //   Nfunc->SetParameters(Norm);
  Nfunc->SetNpx(10000000);
  Nfunc->SetLineColor(1);
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

  //Asymmetry function
  //TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.); 
  //To be used when I use assymetry parameters
    
  TF1* Asy = new TF1("Asy","[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5 + [6]*x^6 + [7]*x^7 + [8]*x^8 + [9]*x^9",0.,1.);
  Asy->SetParNames("ap0","ap1","ap2","ap3","ap4","ap5","ap6","ap7","ap8","ap9");
  Asy->SetParameters(ap0,ap1,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap9);
  
  Asy->SetNpx(10000000);
  Asy->SetLineColor(3);
  ha->Eval(Asy);
  
  //Fit Function
  TF1* fit_func1 = new TF1("fit_func1",fit_func,start,(time)*turn,5);
  fit_func1->SetParNames("Norm","#tau","Asymmetry","#phi","R(ppm)");
  fit_func1->SetParameters(Norm,tau, fAsy, fphase, R);
  //   fit_func1->SetParLimits(fphase,0., 2*pi);
  //fit_func1->FixParameter(3,fphase);
  fit_func1->SetNpx(time);
  fit_func1->SetLineColor(2);
  //probability matrix      
  TH2D* p_matrix = new TH2D("p_matrix","Probability matrix",bin,0.,bin,bin,0.,bin);

  //wiggle plot
  TF1* wiggle = new TF1("wiggle","[0]*(1+[1]*cos(0.001439*x+[2]))*exp(-x/64400)",0.,time*turn);
  wiggle->SetParNames("N","A","phi");
  wiggle->SetParameters(Norm, fAsy, fphase);
  wiggle->SetNpx(10000000);

  TH1D* h_wiggletot = new TH1D("Linear pert","Wiggle(t)",time, 0., time*turn);
  //   TH1D* h_wiggletot_pert = new TH1D("h_wiggletot_pert","perturbed wiggle(t)",time, 0., time*turn);  
  //   TH1D* h_wiggletot_rate = new TH1D("h_wiggletot_rate","Rate vs time",time, 0., time*turn);
  //   TGraph* wiggletot_freq = new TGraph();
  
  //gain function
  TH1D* h_gain = new TH1D("h_gain", "Gain(t)", time, 0., time*turn);
  TF1* ATF_gain_func = new TF1("ATF_gain_func","1-[0]*exp(-x/[1])",0., time*turn);
  ATF_gain_func->SetParNames("C","tau");
  ATF_gain_func->SetParameters(eps[n], tau);
  
  
  TH1D* h_prob = new TH1D("h_prob", "", 100, 0., 1);
  //=======================================================================  
  //ROOT FILE
  //   TFile *wigglefile =new TFile ("wiggle.root","recreate");
  //   TTree t1("Output","roottupla");
  //   t1.Branch("Nt",&Nt,"Nt/D");
  //   t1.Branch("fNt",&fNt,"fNt/D");
  //   TObjArray glist(0);
  //========================================================================
  // non perturbed histogram
  for(int it=0; it<time;it++){
    Double_t Ntot=0;
    Double_t Nt=0;
    Double_t del = rd.Gaus(0.,1.);
    //To insert a perturbation uncomment the line with the function interested in and uncomment line y = y/(gain) (line 175)
    //gain = 1. + (eps*(it*turn))/(time*turn); // linear gain fluctuation
    //gain = 1. + (eps[n]*(time*turn-it*turn)/(time*turn));//NOT USING THIS for y = y/gain; Use for Dave code
    //gain = 1. + (eps[n]*TMath::Exp((-it*turn)/(tau))); //exponential gain fluctutation
    //gain = 1. + (eps[n]*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine Gain fluctuation
    gain = 1. + (TMath::Exp((-it*turn)/(tau)))*(eps[n]*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine + exp Gain fluctuation
    //   gain = 1. + (-3*eps*(TMath::Exp(-(it*turn)/tau) - TMath::Exp(-(it*turn)/3000.))); //ATF gain function
    //   cout<<"Gain = "<<gain<<" Time = "<<it*turn<<endl;
    for(int i=60;i<bin;i++){// 60 is threshold energy for 100 bins
      Double_t Nytot=0;
      for(int j=0; j<bin; j++){
	y = float(j)/bin;	
	y = y*gain;
	//y=y/gain; // Use this with aa code 
        mean = y+(1./(2*bin));
	E = y*Emax; 
        //E = mean*Emax;// For aa code
        resolution = 0.05*sqrt(E);//mine and aa code
        Double_t probt=0;
        Double_t Na = hn->GetBinContent(j);
        Double_t A = ha->GetBinContent(j);
        Double_t phi = hp->GetBinContent(j);
	//prob = (-TMath::Erfc((((float(i+1))/bin)-mean)/resolution)+TMath::Erfc(((float(i)/bin)-mean)/resolution))/(2.);
	//prob = (TMath::Erf((((float(i+1))/bin)-mean)/resolution)-TMath::Erf(((float(i)/bin)-mean)/resolution))/(2.);
	prob = (TMath::Erf((((float(i+1))/bin)-y)/resolution)-TMath::Erf(((float(i)/bin)-y)/resolution))/2.;
	//prob = (TMath::Erf(((float(i)/bin)-mean)/resolution)+1)/2.; // does not work
	probt= probt+prob;
	//p_matrix->Fill(i,j, prob); // Leaving the prob matrix for the time being
	h_prob->SetBinContent(j,prob);
	Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
	Nytot = Nytot+Ny;
      } //end for j
      Ntot = Ntot+Nytot;
    }//end for i
    Nt = Ntot*TMath::Exp(-(it*turn)/tau);
    Nt = (Nt + del*sqrt(Nt));
    h_wiggletot->SetBinContent(it,Nt);
    h_gain->SetBinContent(it, gain);   
    
  }//end for 
  
  h_wiggletot->Fit("fit_func1","R"); 

  R_w[n]=fit_func1->GetParameter(4);
  R_err[n]=fit_func1->GetParError(4);
  

  //  R_errs[n]=sqrt(fit_func1->GetParError(4); //leave random walk
  h_wiggletot->SetXTitle("t (ns)");  
  //printf("R %f\n",Res[n]);
  gStyle->SetOptFit(1);  

}



void plotWiggle(int n=7, double start=0.){
  //double x[7]={0};
  for(int i=0; i<n; i++){
    wiggle(i,start);
    printf("delta_w = %f, ",R_w[i]-(-0.027004));
    
  } 
}
