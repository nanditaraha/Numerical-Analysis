/*
 Author - Nandita Raha
 Usage: root -l
 .x simulate.C(parameter_list)
 */
#include <TH1D.h>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TMath.h>
#include "TApplication.h"
#include <TRandom3.h>
#include <TSpline.h>
/*
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
*/
void simulate(Int_t events=2000,Int_t rep=1){
    //Variables
    Double_t eps=0.001; Double_t bins=140; //Int_t events=2000;
    Double_t fNt;
    Double_t gain = 0.;
    TRandom3 *random = new TRandom3;
    //TRandom3 *random = new TRandom3;
    Int_t no_pts = 140; //in our case 140 for 5 us
    Double_t turn = 149; //considering 1 turn
    Double_t time = 4700; //700 mus if turn 149 ns
    //double x_low=0.;
    //double x_hi=time*turn;
    //function for gain in micro sec units
    TF1* f_gain = new TF1("f_gain","1. + [0]*TMath::Exp(-x/[1])",0,time*turn);//exp gain
    //TF1* f_gain = new TF1("f_gain","1.+[0]*x",0,time*turn);//no gain
    //TF1* f_gain = new TF1("f_gain","1. + [0]*(149*4700-x)/(149*4700)",0,time*turn);//linear gain
    //TF1* f_gain = new TF1("f_gain","1. + [0]*TMath::Cos(0.001439*x + [1])",0,time*turn);//cos gain
    //TF1* f_gain = new TF1("f_gain","1. + [0]*TMath::Exp(-x/[1])*TMath::Cos(0.001439*x + [2])",0,time*turn);//mixed gain

    f_gain->SetParNames("eps","tau");//exp
    //f_gain->SetParNames("eps","tau","phase");//mixed
    //f_gain->SetParNames("eps");//linear
    //f_gain->SetParNames("eps","phase");//    , cos
    
    f_gain->SetParameters(eps,64400);//exp
    //f_gain->SetParameters(eps,64400,4.5e-6);//mixed
    //f_gain->SetParameter(0,eps);//linear,
    //f_gain->SetParameters(eps,4.5e-6);//cos
    //f_gain->SetParameter(0,0);//no SiPM gain
    
    //TF1* f_gauss = new TF1("f_gauss","TMath::Gaus(x,[0],[1],0)",0,2);
    TF1* f_gauss = new TF1("f_gauss","gaus",0,2);
    f_gauss->SetLineColor(kGreen);
    
    
    // double data;
    TH1D* h_gain = new TH1D("h_gain_mod", "", bins, 0,7e5);//simulated gain for 2000 events
    TH1D* h_gain_final = new TH1D("h_gain_final", "", time, 0,time*turn);//final simulated gain ratio for 2000 events interpolated

    TH1D* h_gain_ratio = new TH1D("h_gain_ratio", "", bins, 0,7e5);//final gain for all 2000 cycles
    TH1D* h_gain_norm = new TH1D("h_gain_norm", "", bins, 0,7e5);//final gain for all 2000 cycles
    //  double gain=0;
    double pmt=0,pmt1=0;
    double gain_ratio=0, gain_ratio1=0;
    //Expt:-Double # of bins
    TH1D* h_test;// = new TH1D("h_test", "", 2000, 0.6,1.2);
    TH1D* h_test1;
    char name[50];
    //TF1* f_gauss[140];
    double mean,sigma,data,err,data1;
    // simulating the SiPM's
    
    // simulating the SiPM's
    for(int i=1;i<=140;i++){//time bins (gaps in time of 5 us)
        h_test1 = new TH1D("SiPM Gain Simulation", "", 2000, 0.6,3);
        for (Int_t j = 0; j < events; j++) {//Loop of 2000 simulated points to produce ideal gain
            mean = f_gain->Eval(i*5000);//f_gain is our SiPM exponential function at 5, 10, 15... us i.e 5 us or 5000 ns apart
            for (Int_t k = 0; k < rep; k++){
                //data1 = random->Gaus(mean,0.02*mean);//original plot of gain
                data1 = mean*(1+0.02*mean*random->Gaus(0,1));//original plot of gain
                h_test1->Fill(data1);//original
            }
        }//end of cycles of sim pts
        //cout << i<<" " <<f_gain->Eval(i*5000)<<" "<<0.02*f_gain->Eval(i*5000)<<" NO. o f entries "<<h_test1->GetEntries()<<endl;
        f_gauss->SetParameters(mean,0.02*mean);
        //if(i==141) i=140;
        h_test1->Fit("gaus","Q,0");//"0" not to draw and Q quiet mode - not show fit results
        TF1 *fit2 = (TF1*)h_test1->GetFunction("gaus");
        fit2->SetParError(1,fit2->GetParameter(2)/sqrt(events*rep));
        h_gain_norm->SetBinContent(i,fit2->GetParameter(1));
        h_gain_norm->SetBinError(i,fit2->GetParError(1));
        
        
    }
    //End of SiPM gain simu..
    
    //Simulating PMT gain folded with SiPM
    for(int i=1;i<=140;i++){//time bins (gaps in time of 5 us)
        h_test = new TH1D("SiPM Gain Simulation - PMT's folded", "", 2000, 0.6,3);
        for (Int_t j = 0; j < events; j++) {//Loop of 2000 simulated points to produce ideal gain
            pmt=1+0.02*(j-1)/events;//PMT gain is 2% for 2000 cyles. We use gain for 1 cycle
            //printf("PMT %f\n",pmt);
            //pmt=1+0.003*(j-1)/events;//PMT gain is 0.3% for 2000 cyles use 1 cycles gain from worst case scenario of our optical connectors
            //pmt=1+0.001*(j-1)/events;//0.1%
            //pmt=1+0.0005*(j-1)/events;//0.05%
            
            //pmt=1;//NO PMT gain
            mean = f_gain->Eval(i*5000);//f_gain is our SiPM exponential function at 5, 10, 15... us i.e 5 us or 5000 ns apart
            gain_ratio=mean/pmt;
            sigma = 0.02*gain_ratio;//assume 2% spread in sigma..
            for (Int_t k = 0; k < rep; k++){
                //data = random->Gaus(gain_ratio,sigma);
                data = gain_ratio*(1+sigma*random->Gaus(0,1));
                h_test->Fill(data);
            }
        }
        
        f_gauss->SetParameters(gain_ratio,sigma);
        h_test->Fit("gaus","Q");//"0" not to draw and Q quiet mode - not show fit results
        TF1 *fit1 = (TF1*)h_test->GetFunction("gaus");
        fit1->SetParError(1,fit1->GetParameter(2)/sqrt(events*rep));
        
        h_gain->SetBinContent(i,fit1->GetParameter(1));
        h_gain->SetBinError(i,fit1->GetParError(1));
        
    }
    
    
    
    //h_gain_norm->Draw();
    //h_gain_norm->Fit("f_gain");
    //f_gain->Draw("same");
    
    
    //Second loop for finding/simulationg SiPM and PMT's folded together..
    /*
    for(int i=1;i<=140;i++){//time bins (gaps in time of 5 us)
        h_test = new TH1D("SiPM Gain Simulation - PMT's folded", "", 2000, 0.6,3);
        for (Int_t j = 1; j <= events; j++) {//Loop of 2000 simulated points to produce ideal gain
            //pmt=1+0.02*(j-1)/events;//PMT gain is 2% for 2000 cyles. We use gain for 1 cycle
            //printf("PMT %f\n",pmt);
            //pmt=1+0.003*(j-1)/events;//PMT gain is 0.3% for 2000 cyles use 1 cycles gain from worst case scenario of our optical connectors
            //pmt=1+0.001*(j-1)/events;//0.1%
            //pmt=1+0.0005*(j-1)/events;//0.05%
            pmt=1;//NO PMT gain
            mean = f_gain->Eval(i*5000);//f_gain is our exponential function at 5, 10, 15... us i.e 5 us or 5000 ns apart
            //mean = h_gain_norm->GetBinContent(i);//
            gain_ratio=mean/pmt;
            sigma = 0.02*gain_ratio;//assume 2% spread in sigma..
            data = random->Gaus(gain_ratio,sigma);
            h_test->Fill(data);
        }//end of cycles of sim pts
        
        f_gauss->SetParameters(gain_ratio,sigma);
        h_test->Fit("gaus","Q,0");//"0" not to draw and Q quiet mode - not show fit results
        TF1 *fit1 = (TF1*)h_test->GetFunction("gaus");
        fit1->SetParError(1,fit1->GetParameter(2)/sqrt(events));
        h_gain->SetBinContent(i,fit1->GetParameter(1));
        h_gain->SetBinError(i,fit1->GetParError(1));
        //h_test->Delete();
    }
     */
     
    //h_gain_norm->Draw();
    //f_gain->Draw("same");
    //end of time bins loop to simulate gain function
   
    TCanvas* c3 = new TCanvas("c3");
    c3->Divide(1,3);
    c3->cd(1);
    h_gain_norm->GetYaxis()->SetLabelSize(0.07);
    h_gain_norm->GetXaxis()->SetLabelSize(0.07);
    h_gain_norm->SetStats(0);
    //f_gain->Draw();
    
    //h_gain_norm->SetStats(0);
    h_gain_norm->Draw(); //norm => Normal gain of SiPM's only
    //h_test1->Draw();
    //h_gain->SetLineColor(kRed);
    h_gain->SetStats(0);
    //h_gain->Draw("same");
    
    c3->cd(2);
    //h_gain_norm->SetLineColor(kRed);
    h_gain->GetYaxis()->SetLabelSize(0.07);
    h_gain->GetXaxis()->SetLabelSize(0.07);
    h_gain->GetYaxis()->SetTitleSize(0.09);
    h_gain->GetYaxis()->SetTitle("Gain(t)");
    h_gain->GetYaxis()->CenterTitle();
     
    h_gain->SetStats(0);
    h_gain->Draw();
  
    //h_gain_norm->Fit("f_gain","R","",15000,650000);
    //f_gain->Draw("same");
    c3->cd(3);
    TH1F *h_ratio=(TH1F*)h_gain->Clone("h_ratio");
    h_ratio->Divide(h_gain_norm);
    for(int i=0;i<140;i++){
        double x=h_gain_norm->GetBinError(i+1);
        //double x=h_gain->GetBinError(i+1);
        h_ratio->SetBinError(i+1,x);
    }

    h_ratio->GetYaxis()->SetLabelSize(0.07);
    h_ratio->GetYaxis()->SetTitle(" ");
    h_ratio->GetXaxis()->SetLabelSize(0.07);
    h_ratio->GetXaxis()->SetTitleSize(0.09);
    h_ratio->GetXaxis()->SetTitle("Start time (ns)");
    h_ratio->GetXaxis()->CenterTitle();
    h_ratio->SetStats(0);
    h_ratio->Draw();
    double Point_x[140]={0};
    double Point_y[140]={0};
    for(int i=0;i<140;i++){
        Point_x[i]=h_ratio->GetBinCenter(h_ratio->FindBin(i*5000));
        Point_y[i]=h_ratio->GetBinContent(i+1);
        }
    TSpline3 *sp = new TSpline3("Cubic Spline", Point_x, Point_y, 140, "b2e2", 0, 0);
    // "b2e2" together with the last two "0" means that the second derivatives
    // of the begin and end points equal to zero
    sp->SetLineColor(kRed);
    //sp->Draw("lsame");
    for(int i=0;i<time;i++){
        h_gain_final->SetBinContent(i,sp->Eval(i*turn));
        }
    h_gain_final->SetLineColor(kMagenta);
    h_gain_final->Draw("same");
    //TFile f("long_gain_exp.root","recreate");
    TFile f("junk.root","recreate");//only for trails - otherwise give a sensible name..
    //TFile f("long_gain_linear.root","recreate");//linear
    //TFile f("long_gain_cos.root","recreate");//Cos
    //TFile f("long_gain_mixed.root","recreate");//mixed
    h_gain_final->Write();
    
}
