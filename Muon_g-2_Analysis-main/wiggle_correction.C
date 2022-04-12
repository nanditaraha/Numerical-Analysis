//THIS IS MY FINAL CHANGED VERSION OF ORIGINAL CODE GIVEN BY ANTONIO - this includes correction with lasers..
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

//TFile *g1 =new TFile ("gain.root","recreate");
//TH1F *h1 = (TH1F*)g1->Get("Gain");

void wiggle_correction() {
  //Variables
  Double_t pi = 3.14159265358979323846;
  Double_t Emax = 3.1;
  Int_t bin = 100;
  // Int_t bin = 1000;
  Double_t eps = 0.001;
  Double_t prob;
  Double_t y;
  Double_t pert_y;
  Double_t mean;
  Double_t E;
  Double_t var;
  //Double_t Nt;
  Double_t fNt;
  Double_t gain_correction=0;
  Double_t gain = 0.;
  Double_t turn = 149; //considering 1 turn time in ns
  //   Double_t turn = 1; //considering every ns 
  Double_t Norm = 2.e11; //Total data set
  //   Double_t Norm = 1.15e4; //single fill
  Double_t tau = 64400.;
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  Double_t fphase = 4.5e-6;
  //Double_t fphase = pi/2;
  //Double_t time = 4700; //700 mus if turn 149 ns
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

  //Lab frame number dist and decay asymmetry are from equrion 3.19 of TDR -page 82 chapter 3
  //N function
  TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set
  //TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);//single fill
  //Nfunc->SetParNames("Norm");
  //Nfunc->SetParameters(Norm);
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

  //Asymmetri function
  TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.);
  Asy->SetNpx(10000000);
  Asy->SetLineColor(3);
  ha->Eval(Asy);
  
  //Fit Function
  TF1* fit_func1 = new TF1("fit_func1","[0]*(1+[2]*cos(0.001439*(1-[4]*1.e-6)*x+([3]+2*pi)))*exp(-x/[1])",0.0,time*turn);
  fit_func1->SetParNames("N","#tau","A","#phi","R (ppm)");
  fit_func1->SetParameters(Norm,tau, fAsy, fphase, R);
  //fit_func1->SetParLimits(fphase,0., 2*pi);
  fit_func1->SetNpx(time);
  fit_func1->SetLineColor(2);
    
  //probability matrix
  TH2D* p_matrix = new TH2D("p_matrix","Probability matrix",bin,0.,bin,bin,0.,bin);

  //wiggle function
  //TH1D* h_wiggle = new TH1D("h_wiggle","wiggle normal",4700, 0., time*turn);
  //TH1D* h_wiggley = new TH1D("h_wiggley","wiggle(y)",4700, 0., time*turn);
  TH1D* h_wiggletot = new TH1D("h_wiggletot","wiggle(t)",time, 0., time*turn);
  //TH1D* h_wiggletot_pert = new TH1D("h_wiggletot_pert","perturbed wiggle(t)",time, 0., time*turn);  
  //TH1D* h_wiggletot_rate = new TH1D("h_wiggletot_rate","Rate vs time",time, 0., time*turn);
  //grafico in frequenza
  //TGraph* wiggletot_freq = new TGraph();
    
  //gain function
  TH1D* h_gain = new TH1D("h_gain", "Gain(t)", time, 0., time*turn);
  TF1* ATF_gain_func = new TF1("ATF_gain_func","1-[0]*exp(-x/[1])",0., time*turn);
  ATF_gain_func->SetParNames("C","tau");
  ATF_gain_func->SetParameters(eps, tau);
  
  TH1D* h_prob = new TH1D("h_prob", "", bin, 0., 1);
  //=======================================================================  
  //ROOT FILE
  //TFile *wigglefile =new TFile ("wiggle.root","recreate");
  //TTree t1("Output","roottupla");
  //t1.Branch("Nt",&Nt,"Nt/D");
  //t1.Branch("fNt",&fNt,"fNt/D");
  //TObjArray glist(0);
  //=======================================================================
  //non perturbed histogram  


  //TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*(exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,700.0);
  //f3->FixParameter(0,1.057e-5);
  //f3->FixParameter(1,100.5-6.9);//min value
  //f3->FixParameter(2,16.5);//Fixed from fit results...

  //f3->SetParNames("p_{0}","n_{0}","#tau_{r}");

 
/*
  TFile *g1 = TFile::Open("noWiggle_constE_fac10_g80_m100.root");
  // TH1F *h1 = (TH1F*)g1->Get("Laser");
  TH1F *h1 = (TH1F*)g1->Get("All");
  
  for(int i=0; i<700;i++)
    {
      if(h1->GetBinContent(i+1)!=0)
	//h1->SetBinError(i+1,0.02/sqrt(2000));
	{
	  
	  h1->SetBinContent((i+1),h1->GetBinContent(i+1)+(0.02/sqrt(2000))*rd.Gaus(0.,1.));
	  h1->SetBinError(i+1,0.02/sqrt(2000));
	  //h1->SetBinError(i+1,2e-4);
	}
      else h1->SetBinError(i+1,0.0);
      }
*/
  //This file has the same gain for 320 us which we plan to use
  //Without wiggle below
  
 
  //f2->Draw();
  //This function below is with muons only - G - case of 320 us
  
  TF1 *S_N= new TF1("S_N","24.9828*(exp(-(x+1.2)/64.4)-exp(-(x+1.2)/18.))/64.4", 0,700);
  TF1 *S_a= new TF1("S_a","0.000416279*((64.4-18)*exp(-(x+1.2)/64.4)*cos(1.449*(x+1.2))-exp(-(x+1.2)/18.) + 64.4*18*1.439*exp(-(x+1.2)/64.4)*sin(1.449*(x+1.2)))/64.4", 0,700);
  TF1 *S_b= new TF1("S_b","0.000416279*((64.4-18)*exp(-(x+1.2)/64.4)*sin(1.449*(x+1.2))-exp(-(x+1.2)/18.) + 64.4*18*1.439*exp(-(x+1.2)/64.4)*cos(1.449*(x+1.2)))/64.4", 0,700);

 TF1 *f2 = new TF1("f2","1-100*[0]*([1]*S_N+[2]*S_a+[3]*S_b)", 0,700);
 f2->FixParameter(0,4.26658e-04);
 f2->FixParameter(1,1.44885);
 f2->FixParameter(2,-6.65059e-02 );
 f2->FixParameter(3,-5.72435e-03);
  
  /*
  //This function below is lasers only for correction - G'
  TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*(exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,700.0);
  f3->FixParameter(0,1.01264e-05);
  f3->FixParameter(1,9.62583e+01);
  f3->FixParameter(2,1.78445e+01);//Fixed from fit results...
  f3->SetNpx(1000);
  */
  //f3->SetParNames("p_{0}","n_{0}","#tau_{r}");

  //h1->Fit("f3");
  
  /////////////////////////////////////////////////////////

  
  /////////////////////////////////////////////////////////

  
  //h1->Draw("same");
  //cout<<"Found file\n";
  //printf("Fit tau_r = %f\n",f3->GetParameter(2));
  for(int it=0; it<time;it++){
    Double_t Ntot=0;
    Double_t Nt=0;
    Double_t del = rd.Gaus(0.,1.);
    //Stable 
    //gain = 1;
    //To insert a perturbation uncomment the line with the function interested in and uncomment line y=y/(gain) (line 162)
    //   gain = 1 + 1*eps;
    //  gain = 1. + (-3*eps*(it*turn))/(time*turn); // linear gain fluctuation
    //gain = 0.001;
    //gain = 1. + (eps*(it*turn))/(time*turn); // linear gain fluctuation
    //gain = 1 + (eps*(time*turn-it*turn)/(time*turn));//Linear gain Nandita used
    //gain = 1. + (eps*TMath::Exp((-it*turn)/(tau))); //exponential gain fluctutation
    // gain = 1. + (-3*eps*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine Gain fluctuation
    //gain = 1 -3*eps*(TMath::Exp(-it*turn/tau) - TMath::Exp(-it*turn/18000.)); //ATF gain function
    //gain = h1->GetBinContent(it);

    gain = f2->Eval(it*turn/1000.);//Theoretical gain 
    //cout<<"Gain = "<<gain<<" Time = "<<it*turn<<endl;
    
    //This is only in case we apply a laser correction
    //gain_correction =  f3->Eval(it*turn/1000.);
    
    gain_correction = 1;    

    //if(it*turn<700000.)
    //printf("Gain = %f, Gain correction = %f, Time = %f\n",gain, gain_correction,it*turn);
    for(int i=60;i<bin;i++){
      Double_t Nytot=0;
      for(int j=0; j<bin; j++){
	y = float(j)/bin;
	//y = y+(gain);//add only for a constant offset
	//y=y/gain;
	
	y = y*gain_correction/(gain);
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
	p_matrix->Fill(i,j, prob);
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
    
  }//end for it
  
  //========================================================================================================================================================
  

  //DONE ONLY AFTER EVRTHG IS OK
   
  TCanvas* c2 = new TCanvas("c2");
  //h_wiggletot->Draw();
  fit_func1->Draw();
  fit_func1->FixParameter(1,tau);
  h_wiggletot->Fit("fit_func1","R","",0.,700000.);
  h_wiggletot->SetXTitle("t[ns]");
  h_wiggletot->SetYTitle("Entries");
  gStyle->SetOptFit(1);
  
  //Just to print data;
  //for(int i=0;i<4000;i++)
  //  printf("%.2f,\t",h_wiggletot->GetBinContent(i));
  
}
