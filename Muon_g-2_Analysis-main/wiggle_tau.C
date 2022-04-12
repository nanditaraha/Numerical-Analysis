//This finds effect of a gian changing with varying epsilon. Range of tau is 0 to 2*tau (2*64.4 us)
//Change to histograms to have the names.... #include <TH1D.h>
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
//#include<TAttParticle.h>
//Global Variables and arrays for fit scan with changing tau:
Double_t R_w[500]={0};
Double_t R_err[500]={0};//statistical error

Double_t A[500]={0};
Double_t A_err[500]={0};//statistical error

Double_t t[500]={0};
Double_t t_err[500]={0};//statistical error

Double_t chi[500]={0};
Double_t chi_err[500]={0};//statistical error

Double_t phi[500]={0};
Double_t phi_err[500]={0};//statistical error

Double_t errx[500]={0};

Double_t tau = 64400.;

void wiggle(int n=50, Double_t start=tau) {
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
  Double_t gain = 0.;
  Double_t turn = 149; //considering 1 turn
  //   Double_t turn = 1; //considering every ns 
  Double_t Norm = 2.e11; //Total data set
  //   Double_t Norm = 1.15e4; //single fill
  //Double_t tau = 64400.;
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  Double_t fphase = 4.5e-6;
  //Double_t fphase = pi/2;
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

//Asymmetri function
  TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.);
  Asy->SetNpx(10000000);
  Asy->SetLineColor(3);
  ha->Eval(Asy);
  
//Fit Function
  //TF1* fit_func1 = new TF1("fit_func1","[0]*(1+[2]*cos(0.001439*(1-[4]*1.e-6)*x+([3]+2*pi)))*exp(-x/[1])",30000.0,time*turn);
    //start fitting from 30 us
    TF1* fit_func1 = new TF1("fit_func1","[0]*(1+[2]*cos(0.001439*(1-[4]*1.e-6)*x+([3]+2*pi)))*exp(-x/[1])",0.0,time*turn);
  fit_func1->SetParNames("Norm","tau","Asymmetry","Phase","R(ppm)");
  fit_func1->SetParameters(Norm,tau, fAsy, fphase, R);
//   fit_func1->SetParLimits(fphase,0., 2*pi);
  fit_func1->SetNpx(time);
  fit_func1->SetLineColor(2);
    
//wiggle plot
  TF1* wiggle = new TF1("wiggle","[0]*(1+[1]*cos(0.001439*x+[2]))*exp(-x/64400)",0.,time*turn);
  wiggle->SetParNames("N","A","phi");
  wiggle->SetParameters(Norm, fAsy, fphase);
  wiggle->SetNpx(10000000);

//probability matrix
  //TH2D* p_matrix = new TH2D("p_matrix","Probability matrix",bin,0.,bin,bin,0.,bin);

//wiggle function
//   TH1D* h_wiggle = new TH1D("h_wiggle","wiggle normal",4700, 0., time*turn);
//   TH1D* h_wiggley = new TH1D("h_wiggley","wiggle(y)",4700, 0., time*turn);
  TH1D* h_wiggletot = new TH1D("h_wiggletot","wiggle(t)",time, 0., time*turn);
//   TH1D* h_wiggletot_pert = new TH1D("h_wiggletot_pert","perturbed wiggle(t)",time, 0., time*turn);  
//   TH1D* h_wiggletot_rate = new TH1D("h_wiggletot_rate","Rate vs time",time, 0., time*turn);
//grafico in frequenza
//   TGraph* wiggletot_freq = new TGraph();
    
//gain function
    /*
  TH1D* h_gain = new TH1D("h_gain", "Gain(t)", time, 0., time*turn);
  TF1* ATF_gain_func = new TF1("ATF_gain_func","1-[0]*exp(-x/[1])",0., time*turn);
  ATF_gain_func->SetParNames("C","tau");
  ATF_gain_func->SetParameters(eps, tau);
  
  
  TH1D* h_prob = new TH1D("h_prob", "", 100, 0., 1);
//=======================================================================  
//ROOT FILE
//   TFile *wigglefile =new TFile ("wiggle.root","recreate");
//   TTree t1("Output","roottupla");
//   t1.Branch("Nt",&Nt,"Nt/D");
//   t1.Branch("fNt",&fNt,"fNt/D");
//   TObjArray glist(0);
//========================================================================================================================================================
     */
// non perturbed histogram
for(int it=0; it<time;it++){
  Double_t Ntot=0;
  Double_t Nt=0;
  Double_t del = rd.Gaus(0.,1.);
//To insert a perturbation uncomment the line with the function interested in and uncomment line y=y/(gain) (line 162)
//   gain = 1 + 1*eps;
//  gain = 1. + (-3*eps*(it*turn))/(time*turn); // linear gain fluctuation
//   gain = 1 + (eps*(time*turn-it*turn)/(time*turn));
  gain = 1. + (eps*TMath::Exp((-it*turn)/(start))); //exponential gain fluctutation
//   gain = 1. + (-3*eps*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine Gain fluctuation
//   gain = 1 + (-3*eps*(TMath::Exp(-(it*turn)/tau) - TMath::Exp(-(it*turn)/3000.))); //ATF gain function
//   cout<<"Gain = "<<gain<<" Time = "<<it*turn<<endl;
      for(int i=60;i<bin;i++){
        Double_t Nytot=0;
        for(int j=0; j<bin; j++){
          y = float(j)/bin;
 	  y = y/(gain);
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
     

}
   TCanvas* c2 = new TCanvas("c2");
   h_wiggletot->Draw();
   h_wiggletot->Fit("fit_func1","R");
   h_wiggletot->SetXTitle("t[ns]");
   h_wiggletot->SetYTitle("Entries");
   gStyle->SetOptFit(1);
    
    
    chi[n]=fit_func1->GetChisquare()/fit_func1->GetNDF();
    chi_err[n]=1.0/sqrt(fit_func1->GetNDF());
    
    phi[n]=fit_func1->GetParameter(3);
    phi_err[n]=fit_func1->GetParError(3);
    
    R_w[n]=fit_func1->GetParameter(4);
    R_err[n]=fit_func1->GetParError(4);
    
    A[n]=fit_func1->GetParameter(2);
    A_err[n]=fit_func1->GetParError(2);
    
    t[n]=fit_func1->GetParameter(1);
    t_err[n]=fit_func1->GetParError(1);
}//end of wiggle


void startTauScan(int n=100,int k=50){
    double x[500]={0};
    for(int i=0; i<=n; i++){
        //x[i]=i*tau/50.; //as bin size is 149
        //x[i]=i*tau/10.; //large asymptotic values n>>1 - take n=500 here range is 50*tau
        //x[i]=i*tau; //multiples of tau now
        x[i]=i*tau/k; //time range from 0 to tau/2 - fraction of tau..
        wiggle(i,x[i]);
        
    }
    

    //TH1D *hist_phi = new TH1D("Phase","Phase",n,0,2*tau);//range of 2*tau
    //TH1D *hist_phi = new TH1D("Phase","Phase",n,0,50*tau);//range of 50*tau
    TH1D *hist_phi = new TH1D("Phase","Phase",n,0,n*tau/k);//range of tau/2
    hist_phi->GetXaxis()->SetTitle("Gain time constant(ns)");
    hist_phi->GetXaxis()->CenterTitle();
    hist_phi->GetYaxis()->SetTitle("#phi");
    hist_phi->GetYaxis()->CenterTitle();
   
    hist_phi->SetTitle("");
    for(int i=0;i<n;i++){
        hist_phi->SetBinContent(i,phi[i]);
        hist_phi->SetBinError(i,phi_err[i]);
    }
 
    
    //TGraphErrors *range_chi = new TGraphErrors(n,x,chi,errx,chi_err);
    //TH1D *hist_chi = new TH1D("Chisq","Chisq",n,0,2*tau);
    //TH1D *hist_chi = new TH1D("Chisq","Chisq",n,0,50*tau);
    TH1D *hist_chi = new TH1D("Chisq","Chisq",n,0,n*tau/k);
    hist_chi->GetXaxis()->SetTitle("Gain time constant (ns)");
    hist_chi->GetXaxis()->CenterTitle();
    hist_chi->GetYaxis()->SetTitle("#chi^{2}/NDF");
    hist_chi->GetYaxis()->CenterTitle();
   
    hist_chi->SetTitle("");
     
    for(int i=0;i<n;i++){
        hist_chi->SetBinContent(i,chi[i]);
        hist_chi->SetBinError(i,chi_err[i]);
    }
    //c->cd(3);
    //range_chi->Draw("alp");
    
    
    //TGraphErrors *range_A = new TGraphErrors(n,x,A,errx, A_err);
    //TH1D *hist_A = new TH1D("Asy","Asy",n,0,2*tau);
    //TH1D *hist_A = new TH1D("Asy","Asy",n,0,50*tau);
    TH1D *hist_A = new TH1D("Asy","Asy",n,0,n*tau/k);
    hist_A->GetXaxis()->SetTitle("Gain time constant (ns)");
    hist_A->GetXaxis()->CenterTitle();
    //hist_A->GetYaxis()->SetLabelSize(0.07);
    //hist_A->GetYaxis()->SetTitleSize(0.06);
    hist_A->GetYaxis()->SetTitle("Asymmetry");
    hist_A->GetYaxis()->CenterTitle();
    //hist_A->SetMarkerStyle(22);
    //hist_A->SetMarkerSize(0.3);
    //hist_A->SetMarkerColor(kBlue);
    hist_A->SetTitle("");
    for(int i=0;i<n;i++){
        hist_A->SetBinContent(i,A[i]);
        hist_A->SetBinError(i,A_err[i]);
    }
    //c1->cd(1);
    //range_A->Draw("alp");
    
    //TGraphErrors *range_t = new TGraphErrors(n,x,t,errx,t_err);
    //TH1D *hist_t = new TH1D("tau","tau",n,0,2*tau);
    //TH1D *hist_t = new TH1D("tau","tau",n,0,50*tau);
    TH1D *hist_t = new TH1D("tau","tau",n,0,n*tau/k);
    hist_t->GetXaxis()->SetTitle("Gain time constant (ns)");
    hist_t->GetXaxis()->CenterTitle();
    hist_t->GetYaxis()->SetTitle("#tau");
    hist_t->GetYaxis()->CenterTitle();
   
    hist_t->SetTitle("");
    for(int i=0;i<n;i++){
        hist_t->SetBinContent(i,t[i]);
        hist_t->SetBinError(i,t_err[i]);
    }
    //c1->cd(2);
    //range_t->Draw("alp");
    
    //TH1D *hist_R = new TH1D("delta_w","delta_w",n,0,2*tau);
    //TH1D *hist_R = new TH1D("delta_w","delta_w",n,0,50*tau);
    //TH1D *hist_R = new TH1D("delta_w","delta_w",n,0,n*tau/k);
    TH1D *hist_R = new TH1D("delta_w","delta_w",n,0,n/k);
    hist_R->GetXaxis()->SetTitle("Gain time constant (ns)");

    hist_R->GetXaxis()->CenterTitle();
    //hist_R->GetYaxis()->SetLabelSize(0.07);
    //hist_R->GetYaxis()->SetTitleSize(0.06);
    hist_R->GetYaxis()->SetTitle("R");
    
    hist_R->GetYaxis()->CenterTitle();
    hist_R->SetTitle("");
    //hist_R->Fill(R_w);
    for(int i=0;i<n;i++){
        hist_R->SetBinContent(i,R_w[i]);
        hist_R->SetBinError(i,R_err[i]);
    }
    
    //TFile *f = new TFile("tau10_scan_30.root", "recreate");//0 to 10*tau
    //TFile *f = new TFile("tau_scan_30.root", "recreate");//0 to tau
    //TFile *f = new TFile("high_asymp_scan_30.root", "recreate");//0 to 50*tau
    //TFile *f = new TFile("low_asymp_scan_30.root", "recreate");//0 to 300 us
    TFile *f = new TFile("multiple10_scan.root", "recreate");//n from 1*tau to 10*tau
    //TFile *f = new TFile("test_tau_scan.root", "recreate");//just a test
    hist_t->Write();
    hist_A->Write();
    hist_chi->Write();
    hist_phi->Write();
    //range_R->Write();
    hist_R->Write();
    
}//end of startTauScan;

