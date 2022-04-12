/**** Original Author Antonio  modified by Nandita
Plan: plot gains for different pert types
See if it changes from no pert values by how much & why does it not effect Na, phi and A
Next see how y is affected by gains and take care of ideal detectors..
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
#include "TLegend.h"
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
Int_t bin = 100;
Double_t tau = 64400.;
Double_t turn = 149; //considering 1 turn

Double_t prob_ideal[100][100]={0};
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

//Fit func def...
double fit_func(double *x, double *par){

  return par[0]*(1+par[2]*cos(0.001439*(1-par[4]*1.e-6)*x[0]+(par[3]+2*TMath::Pi())))*exp(-x[0]/par[1]);
}

void fill_ideal(){
  for(int i=0;i<bin;i++)
    for(int j=0;j<bin;j++){ 
      if(i==j)
	prob_ideal[i][j]=1;
      else 
	prob_ideal[i][j]=0;
    }
     
}


//Other function def - like N, A and phi
/*
double N_func(double *x, double *par){
  return 5.e9*(1./3.)*(x[0]*par[0]-1.)*(4.*x[0]*par[0]*x[0]*par[0]-5*x[0]*par[0]-5.)*(3./5.)*0.008536// Total data set  
}
*/
void wiggle(int n=50, double start=0.) {
  //Variables
  //Double_t pi = 3.14159265358979323846;
  Double_t Emax = 3.1;
  //Int_t bin = 1000;
  Double_t eps = 0.001;
  Double_t prob, prob1, prob2, prob3;
  Double_t y=1, y1=1,y2=1;
  Double_t pert_y;
  Double_t mean, mean1, mean2;
  Double_t E, E1, E2;
  Double_t resolution, resolution1, resolution2;
  //Double_t Nt, Nt1, Nt2;
  Double_t fNt;
  //Double_t gain = 0.;
  Double_t gain = 1., gain1 = 1.0, gain2=1.0 ;
  //Double_t turn = 149; //considering 1 turn - defined globally
  //   Double_t turn = 1; //considering every ns 
  Double_t Norm = 2.e11; //Total data set
  //   Double_t Norm = 1.15e4; //single fill
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  //Double_t fphase = TMath::Pi()/2;
  Double_t fphase = 4.5e-5;
  //   Double_t time = 4700; //700 mus if turn 149 ns
  Double_t time = 4000; //600 mus if turn 149 ns
  //   Double_t time = 1000; //150 mus if turn 149 ns
  //   Double_t Ntot = 0;
  TRandom3 rd;
  fill_ideal(); //fill the ideal case probability matrix w/o detector resolution
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
  //TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*([0]*x-1.)*(4.*(x*[0])^2-5*x*[0]-5.)*(3./5.)*0.008536",0,1);// Total data set
  TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set  
  //   TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);//single fill
  //   Nfunc->SetParNames("Norm");
  //   Nfunc->SetParameters(Norm);
  Nfunc->SetNpx(10000000);
  //TF1* N_func = new TF1("N_func",N_func,start,(time)*turn,5);
  TH1D* hn = new TH1D("hn","N(y)",bin,0.,1.); 
  hn->Eval(Nfunc);
  TH2D* hn0 = new TH2D("hn0","N(y,t)",bin,0.,1., time, 0., time*turn);
  TH2D* hn1 = new TH2D("hn1","N(y,t)",bin,0.,1., time, 0., time*turn);
  TH2D* hn2 = new TH2D("hn2","N(y,t)",bin,0.,1., time, 0., time*turn);
  //hn2->Eval(Nfunc);
  
  //phase function
  TF1* phase = new TF1("phase","[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5 + [6]*x^6 + [7]*x^7 + [8]*x^8 + [9]*x^9",0.,1.);
  phase->SetParNames("p0","p1","p2","p3","p4","p5","p6","p7","p8","p9");
  phase->SetParameters(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);
  phase->SetNpx(10000000);
  phase->SetLineColor(2);
  TH1D* hp = new TH1D("hp","Phi(y)",bin,0.,1.);
  TH2D* hp0 = new TH2D("hp0","#phi(y,t)",bin,0.,1., time, 0., time*turn);
  TH2D* hp1 = new TH2D("hp1","#phi(y,t)",bin,0.,1., time, 0., time*turn);
  TH2D* hp2 = new TH2D("hp2","#phi(y,t)",bin,0.,1., time, 0., time*turn);
  hp->Eval(phase);

  //Asymmetry function
  //TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.); 
  TH2D* ha0 = new TH2D("ha0","A(y,t)",bin,0.,1., time, 0., time*turn);
  TH2D* ha1 = new TH2D("ha1","A(y,t)",bin,0.,1., time, 0., time*turn);
  TH2D* ha2 = new TH2D("ha2","A(y,t)",bin,0.,1., time, 0., time*turn);
  //To be used when I use asymmetry parameters
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
  //fit_func1->FixParameter(2,0.4);
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

  TH1D* h_wiggletot = new TH1D("Linear-[1+ #epsilon(T-t)/T]","Linear-[1+ #epsilon(T-t)/T]",time, 0., time*turn);
  TH1D* h_wiggletot_exp = new TH1D("Exp-[1+  #epsilone^{-(t/#tau)}]","Exp-[1+  #epsilone^{-(t/#tau)}]",time, 0., time*turn);  
  TH1D* h_wiggletot_cos = new TH1D("Phase-[1+  #epsilonCos(#omegat+#phi)]","Phase-[1+  #epsilonCos(#omegat+#phi)]",time, 0., time*turn);
  //   TGraph* wiggletot_freq = new TGraph();
  
  //gain function
  TH1D* h_gain = new TH1D("h_gain", "New Linear-[1+ #epsilon(T-t)/T]", time, 0., time*turn);
  TH1D* h_gain1 = new TH1D("h_gain1", "New exp-[1+  #epsilone^{-(t/#tau)}]", time, 0., time*turn);
  TH1D* h_gain2 = new TH1D("h_gain2", "New phase-[1+  #epsilonCos(#omegat+#phi)]", time, 0., time*turn);
  TF1* ATF_gain_func = new TF1("ATF_gain_func","1-[0]*exp(-x/[1])",0., time*turn);
  ATF_gain_func->SetParNames("C","tau");
  ATF_gain_func->SetParameters(eps, tau);
  
  //y function
  TH1D* h_y = new TH1D("h_y", "Gain(t) - stable", bin, 0., 1);
  TH1D* h_y1 = new TH1D("h_y1", "Gain(t) - exp", bin, 0., 1);
  TH1D* h_y2 = new TH1D("h_y2", "Gain(t) - cos", bin, 0., 1);

  
  TH1D* h_prob = new TH1D("h_prob", "", bin, 0., 1);
  //=======================================================================  
  //ROOT FILE
  //   TFile *wigglefile =new TFile ("wiggle.root","recreate");
  //   TTree t1("Output","roottupla");
  //   t1.Branch("Nt",&Nt,"Nt/D");
  //   t1.Branch("fNt",&fNt,"fNt/D");
  //   TObjArray glist(0);
  //========================================================================
  //Gain def. and filling N, A and phi with gain and y...
  for(int it=0; it<time;it++){
    //To insert a perturbation uncomment the line with the function interested in and uncomment line y = y/(gain) (line 175)
    //gain = 1. + (eps*(it*turn))/(time*turn); // linear gain fluctuation
    gain = 1. + (eps*(time*turn-it*turn)/(time*turn));//NOT USING THIS for y = y/gain; Use for Dave code
    gain1 = 1. + (eps*TMath::Exp((-it*turn)/(tau))); //exponential gain fluctutation
    //gain = 1. + (eps*TMath::Exp((-it*turn)/(tau))); //exponential gain fluctutation varying tau
    gain2 = 1. + (eps*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine Gain fluctuation
    //gain = 1. + (TMath::Exp((-it*turn)/(tau)))*(eps*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine + exp Gain fluctuation
    //   gain = 1. + (-3*eps*(TMath::Exp(-(it*turn)/tau) - TMath::Exp(-(it*turn)/3000.))); //ATF gain function
    //cout<<"Gain1 = "<<gain1<<"\tGain2 = "<<gain2<<"\tTime = "<<it*turn<<endl;
    for(int i=0;i<bin;i++){
      //Uncomment only for linear case not stable
      hn0->SetBinContent(i,it,hn->GetBinContent(i)*gain);
      hp0->SetBinContent(i,it,hp->GetBinContent(i)*gain);
      ha0->SetBinContent(i,it,ha->GetBinContent(i)*gain);

      hn1->SetBinContent(i,it,hn->GetBinContent(i)*gain1);
      hn2->SetBinContent(i,it,hn->GetBinContent(i)*gain2);
      hp1->SetBinContent(i,it,hp->GetBinContent(i)*gain1);
      hp2->SetBinContent(i,it,hp->GetBinContent(i)*gain2);
      ha1->SetBinContent(i,it,ha->GetBinContent(i)*gain1);
      ha2->SetBinContent(i,it,ha->GetBinContent(i)*gain2);

      //if(it<3)
	//printf("N(y)*gain at y=%d %f\n",i,hn->GetBinContent(i)*gain1);
    }
  }//loop for hn1 etc. def. end


  for(int it=0; it<time;it++){
    Double_t Ntot=0,Ntot1=0,Ntot2=0;
    Double_t Nt=0,Nt1=0,Nt2=0;
    Double_t del = rd.Gaus(0.,1.);
    //To insert a perturbation uncomment the line with the function interested in and uncomment line y = y/(gain) (line 175)     
    //gain = 1. + (eps*(it*turn))/(time*turn); // linear gain fluctuation                                                       
    gain = 1. + (eps*(time*turn-it*turn)/(time*turn));//NOT USING THIS for y = y/gain; Use for Dave code                       
    gain1 = 1. + (eps*TMath::Exp((-it*turn)/(tau))); //exponential gain fluctutation                                                
    //gain = 1. + (TMath::Exp(eps(-it*turn)/(tau))); //exponential gain fluctutation varying tau                               
    gain2 = 1. + (eps*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine Gain fluctuation                                   
    //gain = 1. + (TMath::Exp((-it*turn)/(tau)))*(eps*TMath::Cos(0.001439*it*turn + fphase/2)); //cosine + exp Gain fluctuation  
    //   gain = 1. + (-3*eps*(TMath::Exp(-(it*turn)/tau) - TMath::Exp(-(it*turn)/3000.))); //ATF gain function (Aaron's fn found expt.)                           
    //cout<<"Gain = "<<gain1<<" Time = "<<it*turn<<endl;  
    for(int i=0.6*bin;i<bin;i++){// 60 is threshold energy for 100 bins
      Double_t Nytot=0, Nytot1=0,Nytot2=0;
      for(int j=0; j<bin; j++){
	y = float(j)/bin;
	y = y*gain;
	y1 = y*gain1;
	y2 = y*gain2;
	//y=y/gain; // Use this with aa code 
        mean = y+(1./(2*bin));
	mean1 = y1+(1./(2*bin));
	mean2 = y2+(1./(2*bin));

	//E = y*Emax; 
	//E1 = y1*Emax; 
	//E2 = y2*Emax; 
        
	E = mean*Emax;// For aa code
	E1 = mean1*Emax;
	E2 = mean2*Emax;
        resolution = 0.05*sqrt(E);//mine and aa code
	resolution1 = 0.05*sqrt(E1);
	resolution2 = 0.05*sqrt(E2);

        Double_t probt=0;
	Double_t Na=0, A=0, phi=0, Na1=0, A1=0, phi1=0, Na2=0, A2=0, phi2=0;
	//ideal case

	prob = prob_ideal[i][j];
	prob1 = prob_ideal[i][j];
	prob2 = prob_ideal[i][j];

	//prob = (-TMath::Erfc((((float(i+1))/bin)-mean)/resolution)+TMath::Erfc(((float(i)/bin)-mean)/resolution))/(2.);
	//prob1 = (-TMath::Erfc((((float(i+1))/bin)-mean1)/resolution1)+TMath::Erfc(((float(i)/bin)-mean1)/resolution1))/(2.);
	//prob2 = (-TMath::Erfc((((float(i+1))/bin)-mean2)/resolution2)+TMath::Erfc(((float(i)/bin)-mean2)/resolution2))/(2.);// USING THIS ONE...
	
	//prob = (TMath::Erf((((float(i+1))/bin)-y)/resolution)-TMath::Erf(((float(i)/bin)-y)/resolution))/2.;// I use this
	//prob1 = (TMath::Erf((((float(i+1))/bin)-y1)/resolution1)-TMath::Erf(((float(i)/bin)-y1)/resolution1))/2.;// I use this
	//prob2 = (TMath::Erf((((float(i+1))/bin)-y2)/resolution2)-TMath::Erf(((float(i)/bin)-y2)/resolution2))/2.;// I use this
	
	//stable case below:
	/*
	Na = hn->GetBinContent(j);//Evalauted from 0 time 
	A = ha->GetBinContent(j);  
	phi = hp->GetBinContent(j);                                                                                      
	
        Na1 = hn->GetBinContent(j);//bin counted from 1 and not zero while setting bin content                                    
        A1 = ha->GetBinContent(j);                                                                                                
        phi1 = hp->GetBinContent(j);                                                                                              
                                                                                                                                      
        Na2 = hn2->GetBinContent(j);                                                                                               
        A2 = ha2->GetBinContent(j);                                                                                                
        phi2 = hp2->GetBinContent(j);
	*/
	
	//Linear case below:
	Na = hn0->GetBinContent(j,it);//Evalauted from 0 time
	A = ha0->GetBinContent(j,it);
        phi = hp0->GetBinContent(j,it);
	
	//These are for exp and cos cases.
	Na1 = hn1->GetBinContent(j,it);//bin counted from 1 and not zero while setting bin content 
	A1 = ha1->GetBinContent(j,it); 
	phi1 = hp1->GetBinContent(j,it);                                                                               
	
	Na2 = hn2->GetBinContent(j,it);                                                                            
	A2 = ha2->GetBinContent(j,it);                                                                              
	phi2 = hp2->GetBinContent(j,it);
	//if(it==15) printf("Na=%f, Na1=%f, Na2=%f\n", Na, Na1, Na2);
	
	probt= probt+prob;
	p_matrix->Fill(i,j, prob); // Leaving the prob matrix for the time being
	h_prob->SetBinContent(j,prob);

	Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
	Double_t Ny1 = Na1*(1+A1*TMath::Cos(0.001439*it*turn + phi1*1.e-3))*prob1;
	Double_t Ny2 = Na2*(1+A2*TMath::Cos(0.001439*it*turn + phi2*1.e-3))*prob2;

	Nytot = Nytot+Ny;
	Nytot1 = Nytot1+Ny1;
	Nytot2 = Nytot2+Ny2;
      } //end for j
      Ntot = Ntot+Nytot;
      Ntot1 = Ntot1+Nytot1;
      Ntot2 = Ntot2+Nytot2;
    }//end for i
    Nt = Ntot*TMath::Exp(-(it*turn)/tau);
    Nt1 = Ntot1*TMath::Exp(-(it*turn)/tau);
    Nt2 = Ntot2*TMath::Exp(-(it*turn)/tau);
    
    Nt = (Nt + del*sqrt(Nt));
    Nt1 = (Nt1 + del*sqrt(Nt1));
    Nt2 = (Nt2 + del*sqrt(Nt2));

    //if (Nt!=0)
    h_wiggletot->SetBinContent(it,Nt);
    h_wiggletot_exp->SetBinContent(it,Nt1);
    h_wiggletot_cos->SetBinContent(it,Nt2);
    //printf("Time %d: exp = %f, cos=%f\n",it,Nt1,Nt2);
    
    h_gain->SetBinContent(it, gain);   
    h_gain1->SetBinContent(it, gain1);   
    h_gain2->SetBinContent(it, gain2);   
    //cout<<"Prob = "<<prob<<" Time = "<<it*turn<<endl;
  }//end for 
  //----------------------------------------------------------------------------------------------------------------------------------

  //==================================================================================================================================
  //Removed the matrix and gain plots for now....
  
  TCanvas* c1 = new TCanvas("c1");
  /*
  c1->Divide(1,2);
  c1->cd(1);
  //p_matrix->Draw("LEGO");
  // h_prob->Draw();
  h_y1->SetLineColor(kRed);
  h_y2->SetLineColor(kGreen);
  h_y->Draw();
  h_y1->Draw("same");
  h_y2->Draw("same");
  c1->cd(2);
  //    h_wiggletot_rate->Draw(); 
  //    h_wiggletot_rate->SetXTitle("t[ns]");  
  //    h_wiggletot_rate->SetYTitle("Rate [Hz]");   
  */
  h_gain->SetStats(0);
  h_gain->Draw();
  h_gain1->SetLineColor(kRed);
  h_gain2->SetLineColor(kGreen);
  h_gain1->Draw("same");
  h_gain2->Draw("same");
  h_gain->SetXTitle("Time (ns)");
  c1->BuildLegend();

/*
  TCanvas* c3 = new TCanvas("c3");
  c3->Divide(1,3);
  c3->cd(1);
  hn1->SetLineColor(kRed);
  hn2->SetLineColor(kGreen);
  hn->Draw();
  //hn1->Draw();
  //hn2->Draw("same");
  
  Nfunc->Draw("same");
  c3->cd(2);
  hp->Draw();
  phase->Draw("same");
  c3->cd(3);
  ha->Draw();
  Asy->Draw("same");
  */

  TCanvas* c = new TCanvas("c");
  h_wiggletot->SetStats(0);
  h_wiggletot->SetLineWidth(2);
  h_wiggletot->Draw(); 
  h_wiggletot_exp->SetLineColor(kRed);
  h_wiggletot_cos->SetLineColor(kGreen);
  h_wiggletot_cos->Draw("same"); 
  //h_wiggletot_cos->Draw(); 
  h_wiggletot_exp->SetLineStyle(2);
  h_wiggletot_exp->Draw("same"); 
  //h_wiggletot_exp->Draw(); 
  //fit_func1->SetLineColor(kBlue);
  //h_wiggletot->Fit("fit_func1","R"); 
  //h_wiggletot_cos->Fit("fit_func1","R"); 
  //h_wiggletot_exp->Fit("fit_func1","R"); 
  c->BuildLegend();
  
  T[n]=fit_func1->GetParameter(1);
  T_err[n]=fit_func1->GetParError(1);

  chi[n]=fit_func1->GetChisquare()/fit_func1->GetNDF();
  chi_err[n]=1.0/sqrt(fit_func1->GetNDF());

  A[n]=fit_func1->GetParameter(2);
  A_err[n]=fit_func1->GetParError(2);

  phi[n]=fit_func1->GetParameter(3);
  phi_err[n]=fit_func1->GetParError(3);

  R_w[n]=fit_func1->GetParameter(4);
  R_err[n]=fit_func1->GetParError(4);
  

  //  R_errs[n]=sqrt(fit_func1->GetParError(4); //leave random walk
  h_wiggletot->SetXTitle("t (ns)");  
  printf("GAin %f prob %f\n",gain,prob);
  gStyle->SetOptFit(1);  

}



void startScan(int n=50, double start=0.){
  double x[500]={0};
  for(int i=0; i<n; i++){
    x[i]=start + i*149; //as bin size is 149
    wiggle(i,x[i]);
    
  } 
  TCanvas* c = new TCanvas("c", "", 600,640);
  c->Divide(1,3);

  TCanvas* c1 = new TCanvas("c1", "", 600,640);
  c1->Divide(1,2);

  TGraphErrors *range_A = new TGraphErrors(n,x,A,errx,A_err);
  //range_A->GetXaxis()->SetTitle("Start Time (ns)");
  range_A->GetXaxis()->CenterTitle();
  range_A->GetXaxis()->SetLabelSize(0.09);  
  range_A->GetXaxis()->SetTickLength(0.07);
  range_A->GetYaxis()->SetLabelSize(0.07);
  range_A->GetYaxis()->SetTitleSize(0.09);
  range_A->GetYaxis()->SetTitle("Asym");
  range_A->GetYaxis()->CenterTitle();
  range_A->SetMarkerStyle(22);
  range_A->SetMarkerSize(0.3);
  range_A->SetMarkerColor(kBlue);
  range_A->SetTitle("");
  //range_A->SetTitle("Stable - Case (NR)");
  //range_A->SetTitle("Exponential perturbation");
  c->cd(1);  
  gPad->SetGridx();
  range_A->Draw("alp");
  
  TGraphErrors *range_T = new TGraphErrors(n,x,T,errx,T_err);
  //range_T->GetXaxis()->SetTitle("Start Time (ns)");
  range_T->GetXaxis()->CenterTitle();
  range_T->GetXaxis()->SetLabelSize(0.09);
  range_T->GetXaxis()->SetTickLength(0.07);
  range_T->GetYaxis()->SetTitleSize(0.09);
  range_T->GetYaxis()->SetLabelSize(0.07);
  range_T->GetYaxis()->SetTitle("#tau");
  range_T->GetYaxis()->CenterTitle();
  range_T->SetMarkerStyle(22);
  range_T->SetMarkerSize(0.3);
  range_T->SetMarkerColor(kBlue);
  range_T->SetTitle("");
  c->cd(2);
  range_T->Draw("alp");


  TGraphErrors *range_phi = new TGraphErrors(n,x,phi,errx,phi_err);
  range_phi->GetXaxis()->SetTitle("Start Time (ns)");
  range_phi->GetXaxis()->CenterTitle();
  range_phi->GetXaxis()->SetLabelSize(0.09);
  range_phi->GetXaxis()->SetTickLength(0.07);
  range_phi->GetYaxis()->SetTitleSize(0.09);
  range_phi->GetYaxis()->SetLabelSize(0.07);
  range_phi->GetYaxis()->SetTitle("#phi");
  range_phi->GetYaxis()->CenterTitle();
  range_phi->SetMarkerStyle(22);
  range_phi->SetMarkerSize(0.3);
  range_phi->SetMarkerColor(kBlue);
  range_phi->SetTitle("");
  c->cd(3);
  range_phi->Draw("alp");

  TGraphErrors *range_R = new TGraphErrors(n,x,R_w,errx,R_err);
  range_R->GetXaxis()->SetTitle("Start Time (ns)");
  range_R->GetXaxis()->SetTitleSize(0.07);
  range_R->GetXaxis()->SetLabelSize(0.09);
  range_R->GetXaxis()->SetTickLength(0.07);
  range_R->GetXaxis()->SetLabelSize(0.09);
  range_R->GetXaxis()->CenterTitle();
  range_R->GetYaxis()->SetLabelSize(0.07);
  range_R->GetYaxis()->SetTitleSize(0.06);
  range_R->GetYaxis()->SetTitle("R");
  range_R->GetYaxis()->CenterTitle();
  range_R->SetMarkerStyle(22);
  range_R->SetMarkerSize(0.3);
  range_R->SetMarkerColor(kBlue);
  range_R->SetTitle("");
  c1->cd(1);
  range_R->Draw("alp");

  TGraphErrors *range_chi = new TGraphErrors(n,x,chi,errx,chi_err);
  range_chi->GetXaxis()->SetTitle("Start Time (ns)");
  range_chi->GetXaxis()->SetTitleSize(0.07);
  range_chi->GetXaxis()->SetLabelSize(0.06);
  range_chi->GetXaxis()->SetTickLength(0.07);
  range_chi->GetXaxis()->SetLabelSize(0.09);
  range_chi->GetXaxis()->CenterTitle();
  range_chi->GetYaxis()->SetLabelSize(0.07);
  range_chi->GetYaxis()->SetTitleSize(0.09);
  range_chi->GetYaxis()->SetTitle("#chi^{2}");
  range_chi->GetYaxis()->CenterTitle();
  range_chi->SetMarkerStyle(22);
  range_chi->SetMarkerSize(0.3);
  range_chi->SetMarkerColor(kBlue);
  range_chi->SetTitle("");
  c1->cd(2);
  range_chi->Draw("alp");


  
}


/* Notes:
Both codes do not show any error in asym and phi which is a bit strange. Daves presentation slide 8 shows errors but his histogram e.g. slide 7 does not??.
Dave's code bugs: var or resolution (sigma) = 0.05/sqrt(E) - should be multiplied I think.
The weight which is obtained from the the erf is not correct I think. So I used the erfc of Antonios code and Dave's pseudo code which he explains in his presentation in the logic of my own code. My results match with Antonios at least in order of magnitude especially for the exponential case.....   
 */
