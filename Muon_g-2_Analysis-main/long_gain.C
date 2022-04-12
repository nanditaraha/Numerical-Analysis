/* Few changes
Fill time = 700 us..
*/
#include <TH1D.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include "TFitResult.h"
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


void wiggle() {
//Variables
  Double_t pi = 3.14159265358979323846;
  Double_t Emax = 3.1;
  Int_t bin = 100;
  //Int_t bin = 1000;
  Double_t eps = 0.001;//for one fill;
  Double_t prob;
  Double_t y;
  Double_t pert_y;
  Double_t mean;
  Double_t E;
  Double_t var;
  //Double_t Nt;
  Double_t fNt;
  Double_t gain = 0.;
  Double_t gain_array[20][5000]={0};
  Double_t turn = 149; //considering 1 turn
  //   Double_t turn = 1; //considering every ns 
  Double_t Norm = 2.e11; //Number in Total data set 
  Double_t tau = 64400.;
  Double_t fAsy = 0.4;
  Double_t R = 0.0;
  //Double_t fphase = pi/2;
  Double_t fphase = 4.5e-6;
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
    //Simulated gain function from long_gain.root
    //TFile *f = TFile::Open("long_gain_exp5_1_2%.root");//exp with 2% fluctuation in gain
    //TFile *f = TFile::Open("long_gain_exp.root");//exp with 0.3% fluctuation in gain - worst case scenario
    //TFile *f = TFile::Open("long_gain_linear.root");//linear
  //TFile *f = TFile::Open("long_gain_cos.root");//cos
  //TFile *f = TFile::Open("long_gain_mixed.root");//mixed
  TFile *f = TFile::Open("junk.root");//trails
    TH1D* h_gain_final = (TH1D*)f->Get("h_gain_final");

  //N function
  TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set
  //TF1* Nfunc = new TF1("Nfunc","[0]*[1]*1.15e4*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);//single fill*20*2000
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
  TF1* Asy = new TF1("Asy","(-8.*x^2+x+1.)/(4.*x^2-5.*x-5.)",0.,1.);
  TH1D* ha = new TH1D("ha","A(y)",bin,0.,1.);
  Asy->SetNpx(10000000);
  Asy->SetLineColor(3);
  ha->Eval(Asy);
  
//Fit Function
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
  //wiggle->SetNpx(10000000);

//probability matrix
  TH2D* p_matrix = new TH2D("p_matrix","Probability matrix",bin,0.,bin,bin,0.,bin);

//wiggle function
//   TH1D* h_wiggle = new TH1D("h_wiggle","wiggle normal",4700, 0., time*turn);
//   TH1D* h_wiggley = new TH1D("h_wiggley","wiggle(y)",4700, 0., time*turn);
  TH1D* h_wiggletot = new TH1D("h_wiggletot","wiggle(t)",time, 0., time*turn);
//   TH1D* h_wiggletot_pert = new TH1D("h_wiggletot_pert","perturbed wiggle(t)",time, 0., time*turn);  
//   TH1D* h_wiggletot_rate = new TH1D("h_wiggletot_rate","Rate vs time",time, 0., time*turn);
//   TGraph* wiggletot_freq = new TGraph();
    
//gain function
  TH2D* h_gain = new TH2D("h_gain", "Gain(t)", time, 0., time*turn, 20, 0.,20.);
  TF1* ATF_gain_func = new TF1("ATF_gain_func","1-[0]*exp(-x/[1])",0., time*turn);
  ATF_gain_func->SetParNames("C","tau");
  ATF_gain_func->SetParameters(eps, tau);
    
  TH1D* h_prob = new TH1D("h_prob", "", bin, 0., 1);

// non perturbed histogram
  for(int it=0; it<time;it++){
    Double_t Ntot=0;
    Double_t Nt=0;
    Double_t del = rd.Gaus(0.,1.);
    
    gain=h_gain_final->GetBinContent(it);//REading from the file for long term gain...
    //h_gain->SetBinContent(it,no_fills,gain);
      //gain = 1;
    for(int i=0.6*bin;i<bin;i++){//threshold is at y=0.6 
      Double_t Nytot=0;
      for(int j=0; j<bin; j++){//this loop is for reconstruction considering probability distributions...
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
	h_prob->SetBinContent(j,prob);
	Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
	Nytot = Nytot+Ny;
	//
      } //end for j
      Ntot = Ntot+Nytot;
    }//end for i
    Nt = Ntot*TMath::Exp(-(it*turn)/tau);
    Nt = (Nt + del*sqrt(Nt));//statistical fluctuations 
    //Nt = (Nt + sqrt(Nt));//statistical fluctuations for Poisson distribution...      
    h_wiggletot->SetBinContent(it,Nt);      
    //h_gain->SetBinContent(it, gain);
  }//end for it
 
 

  TCanvas* c2 = new TCanvas("c2");
  //h_gain_final->Draw();
  h_wiggletot->Draw();
  //h_gain->Draw("colz");
  //fit_func1->Draw("same");
  //fit_func1->Draw();
    TFitResultPtr r=h_wiggletot->Fit("fit_func1","RS","",0,time*turn);//gives correlations b/w parameters..
    printf("%f+%f\n",-0.027-fit_func1->GetParameter(4), 0.09+fit_func1->GetParError(4));
    //r->Print("V");
  //h_wiggletot->Fit("fit_func1","R");
  //h_wiggletot->SetXTitle("t[ns]");
  //h_wiggletot->SetYTitle("Entries");
  gStyle->SetOptFit(1);
  
}
