/*
 Author - Nandita Raha
 Usage: root -l
 .L simulate_SiPM.C
 simulate()
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


// Variables used for all functions - thus defined here..


void simulate(Int_t events=2000, Int_t rep=1){
    //Variables
    Double_t eps=0.002; Double_t bins=140; //Int_t events=2000;
    Double_t fNt;
    Double_t gain = 0.;
    TRandom3 *random = new TRandom3;
    gRandom->SetSeed(0);
    //TRandom3 *random = new TRandom3;
    Int_t no_pts = 140; //in our case 140 for 5 us
    Double_t turn = 149; //considering 1 turn
    Double_t time = 4700; //700 mus if turn 149 ns
    //double x_low=0.;
    //double x_hi=time*turn;
    //function for SiPM gain in micro sec units
   
    TF1* f_gain1 = new TF1("f_gain1","1/[0] + [1]*TMath::Exp(-x/[2])",0,time*turn);//exp gain used to fit after adding PMT factor
    TF1* f_gain = new TF1("f_gain","1 + [0]*TMath::Exp(-x/[1])",0,time*turn);//exp gain
    //TF1* f_gain = new TF1("f_gain","1.+[0]*x",0,time*turn);//no gain
    //TF1* f_gain = new TF1("f_gain","1. + [0]*(149*4700-x)/(149*4700)",0,time*turn);//linear gain
    //TF1* f_gain = new TF1("f_gain","1. + [0]*TMath::Cos(0.001439*x + [1])",0,time*turn);//cos gain
    //TF1* f_gain = new TF1("f_gain","1. + [0]*TMath::Exp(-x/[1])*TMath::Cos(0.001439*x + [2])",0,time*turn);//mixed gain

    f_gain->SetParNames("eps","tau");//exp
    f_gain1->SetParNames("PMT gain","eps","tau");//exp
    //f_gain->SetParNames("eps","tau","phase");//mixed
    //f_gain->SetParNames("eps");//linear
    //f_gain->SetParNames("eps","phase");//    , cos
    
    //f_gain1->SetParameters(1.02,eps,64400);//exp
    f_gain1->SetParameters(1,eps,64400);//exp
    f_gain->SetParameters(eps,64400);//exp
    //f_gain->SetParameters(eps,64400,4.5e-6);//mixed
    //f_gain->SetParameter(0,eps);//linear,
    //f_gain->SetParameters(eps,4.5e-6);//cos
    //f_gain->SetParameter(0,0);//no SiPM gain
    
    //TF1* f_gauss = new TF1("f_gauss","TMath::Gaus(x,[0],[1],0)",0,2);
    TF1* f_gauss = new TF1("f_gauss","gaus",0,2);
    f_gauss->SetLineColor(kGreen);
    
    
    // double data;
    TH1D* h_gain = new TH1D("SiPM/PMT gain", "", bins, 0,7e5);//simulated gain for 2000 events
    TH1D* h_gain_final = new TH1D("h_gain_final", "", time, 0,time*turn);//final simulated gain ratio for 2000 events interpolated
    TH1D* h_pmt = new TH1D("PMT gain", "", events, 0,events);//simulated gain for 2000 events
    
    
    TH1D* h_gain_ratio = new TH1D("h_gain_ratio", "", bins, 0,7e5);//final gain for all 2000 cycles
    TH1D* h_gain_norm = new TH1D("SiPM gain", "", bins, 0,7e5);//final gain for all 2000 cycles
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
    for(int i=1;i<=140;i++){//time bins (gaps in time of 5 us)
        h_test1 = new TH1D("SiPM Gain Simulation", "", 2000, 0.6,3);
    	for (Int_t j = 0; j < events; j++) {//Loop of 2000 simulated points to produce ideal gain
            mean = f_gain->Eval(i*5000);//f_gain is our SiPM exponential function at 5, 10, 15... us i.e 5 us or 5000 ns apart
            //mean=1;//Careful: Comment this always - only to check the effect of no SiPM changes....
            for (Int_t k = 0; k < rep; k++){
                //data1 = random->Gaus(mean,0.02*mean);//original plot of gain
                data1 = mean*(1+0.02*mean*random->Gaus(0,1));//original plot of gain
                h_test1->Fill(data1);//original
            }
        }//end of cycles of sim pts (events loop)
        //cout << i<<" " <<f_gain->Eval(i*5000)<<" "<<0.02*f_gain->Eval(i*5000)<<" NO. o f entries "<<h_test1->GetEntries()<<endl;
        f_gauss->SetParameters(mean,0.02*mean);
        //if(i==141) i=140;
        h_test1->Fit("gaus","Q,0");//"0" not to draw and Q quiet mode - not show fit results
        //h_test1->Draw();
        TF1 *fit2 = (TF1*)h_test1->GetFunction("gaus");
        fit2->SetParError(1,fit2->GetParameter(2)/sqrt(events*rep));
        h_gain_norm->SetBinContent(i,fit2->GetParameter(1));
        h_gain_norm->SetBinError(i,fit2->GetParError(1));
    }
    //End of SiPM gain simu. of time bins
    
    //Simulating PMT gain folded with SiPM
    /*
    for(int i=1;i<=140;i++){//time bins (gaps in time of 5 us)
        h_test = new TH1D("SiPM Gain Simulation - PMT's folded", "", 2000, 0.6,3);
        for (Int_t j = 0; j < events+1; j++) {//Loop of 2000 simulated points to produce ideal gain
            //pmt=1.+ 0.02*j/events;//PMT gain is 2% for 2000 cycles. We use gain for 1 cycle = 1e-5
            //printf("PMT %f\n",pmt);
            //pmt=1+0.003*j/events;//PMT gain is 0.3% for 2000 cycles use 1 cycles gain from worst case scenario of our optical connectors
            //pmt=1.02+0.001*j/2000.;//0.1% //skipping 20 events
            pmt=1+0.001*j/events;//0.1%
            //pmt=1.+0.0005*j/events;//0.05%
            
            //pmt=1;//NO PMT gain
            mean = f_gain->Eval(i*5000);//f_gain is our SiPM exponential function at 5, 10, 15... us i.e 5 us or 5000 ns apart
            //mean=1;//Careful: Comment this always - only to check the effect of no SiPM changes....
            gain_ratio=mean/pmt;
            //h_pmt->SetBinContent(j,pmt);
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
        //h_test->Draw();
        //printf("Error: %f\n",fit1->GetParError(1));
        fit1->SetParError(1,fit1->GetParameter(2)/sqrt(events*rep));
        h_gain->SetBinContent(i,fit1->GetParameter(1));
        h_gain->SetBinError(i,fit1->GetParError(1));
        
    }
    
    */
    
    //h_gain_norm->Draw();
    h_gain_norm->Fit("f_gain");//with SiPMs only
    //h_gain->Fit("f_gain1");//with SiPM+PMT
    //h_gain->SetStats(0);
    //h_gain->Draw();
    //h_pmt->Draw();
    //gStyle->SetOptFit(1);
    
}


