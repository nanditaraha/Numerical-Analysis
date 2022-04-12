/*Simulating exp gain function for laser triggers at 5 us for SiPM
Trying with getrandom and with fillrandom
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

void simulate(double wt) {

  //Variables
  Double_t fNt;
  Double_t gain = 0.;
  Double_t turn = 149; //considering 1 turn 
  Double_t time = 4700; //700 mus if turn 149 ns
  TRandom3 *random = new TRandom3;
  int n = 100;//simulating with 1000000 pts
  int no_pts = 140; //in our case 140 for 5 us
  double x_low=0.;
  double x_hi=time*turn;
  //function for gain in micro sec units
  TF1* f_gain = new TF1("f_gain","1. + [0]*TMath::Exp(-x/[1])",x_low,x_hi);
  f_gain->SetParNames("eps","tau");
  f_gain->SetParameters(100,64400.0);
  
  double data;
  //TH1D* h_gain[140];
  double sum_f=0;
  TH1D* h_gain_tot = new TH1D("h_gain_tot", "", wt*time, x_low,x_hi); 
  char name[50]="h_gain";
  for(int i=0;i<140;i++){
    //f_gain->SetRange(5000*(i+1)-149,5000*(i+1)+149);
    //TF1* f_gain1 = new TF1("f_gain1","1. + [0]*TMath::Exp(-x/[1])",5000*(i+1)-149,5000*(i+1)+149);
    int start = h_gain_tot->FindBin(5000*(i+1)-149);
    int stop = h_gain_tot->FindBin(5000*(i+1)+149);
    for (Int_t j = 0; j < n; j++) {
      //double pts = random->Uniform(5000*(i+1)-149,5000*(i+1)+149);      
      double data = f_gain->GetRandom(5000*(i+1)-149,5000*(i+1)+149);
      h_gain_tot->Fill(data,gRandom->Uniform(5000*(i+1)-149,5000*(i+1)+149)/(turn*n));
      //h_gain_tot->Fill(data,f_gain->Integral(x_low,x_hi)/(n*140*time/3.));
      //h_gain_tot->Fill(data,wt*f_gain->Integral(5000*(i+1)-149,5000*(i+1)+149)/(turn*n));//works best till now
      //h_gain_tot->FillRandom("f_gain1",n);
      //h_gain_tot->Scale(f_gain->Eval(5000*(i+1))/h_gain_tot->Integral(start,stop));
    }
  }
  f_gain->SetRange(x_low,x_hi);
  //h_gain_tot->Scale(f_gain->Integral(x_low,x_hi)/(n*140*time/3.));
  //h_gain_tot->Scale(f_gain->Eval(5000)/h_gain_tot->GetBinContent(h_gain_tot->FindBin(5000)));
  h_gain_tot->Draw();
  //h_gain_tot->Fit("f_gain","R");
  f_gain->Draw("same");
    
}
