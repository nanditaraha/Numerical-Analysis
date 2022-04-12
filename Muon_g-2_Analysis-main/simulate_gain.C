/*
 Author - Nandita Raha
 Usage: root -l
 .L simulate_gain.C
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


void wiggle(double eps, double tau);
void wiggle_correct(double eps, double tau);
void simulate(Int_t events=2000, Int_t rep=1){
    //Variables
    Double_t eps=0.001; Double_t bins=140; //Int_t events=2000;
    Double_t fNt;
    Double_t gain = 0.;
    TRandom3 *random = new TRandom3;
    random->SetSeed(0);

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
                //data1 = mean*(1+0.02*random->Gaus(0,1));//original plot of gain
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
    //printf("Para 1 %f",f_gain->GetParameter(0));
    //wiggle(f_gain->GetParameter(0),f_gain->GetParameter(1));//for fitting h_gain_norm (SiPM only)
    //wiggle(f_gain1->GetParameter(1),f_gain1->GetParameter(2));//for fitting h_gain (SiPM+PMT)
    //printf("Eps %f Tau %f\n", f_gain1->GetParameter(1),f_gain1->GetParameter(2));
    //wiggle(0.001,64400);
    //TCanvas* c2 = new TCanvas("c2");;
    //wiggle_correct(f_gain1->GetParameter(1),f_gain1->GetParameter(2)+f_gain1->GetParError(2));// add subtract error
    //wiggle_correct(f_gain1->GetParameter(1),f_gain1->GetParameter(2));
    //wiggle_correct(f_gain->GetParameter(0),f_gain->GetParameter(1));//SiPM Only
    //wiggle_correct(f_gain->GetParameter(0),f_gain->GetParameter(1)+f_gain->GetParError(1));
    //wiggle_correct(0.001,64400);
}

//The wiggle plot created and fitted....


void wiggle(double eps, double tau) {
    
    //printf("Wiggle: Eps %f Tau %f\n",eps, tau);
    //Variables
    Double_t pi = 3.14159265358979323846;
    Double_t Emax = 3.1;
    Int_t bin = 100;
    //Int_t bin = 1000;
    //Double_t eps1 = 0.001;//for one fill;
    Double_t prob;
    Double_t prob1;
    Double_t y;
    Double_t mean;
    Double_t E;
    Double_t var;
    
    Double_t fNt;
    Double_t gain = 1.;
    Double_t gain1 = 1.;
    //Double_t gain_correct=1;
    Double_t turn = 149; //considering 1 turn
    //   Double_t turn = 1; //considering every ns
    Double_t Norm = 2.e11; //Number in Total data set
    //Double_t Norm = 4e8; //Sometimes the fit fails then use this as it closer to N_0 of fit results especially for 2000 event simulations
    Double_t tau1 = 64400.;
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
    
    //N function
    TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set - eq 3.19
    //TF1* Nfunc = new TF1("Nfunc","7.e10*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set - eq 3.19
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
    TH1D* h_wiggletot_pert = new TH1D("h_wiggletot_pert","perturbed wiggle(t)",time, 0., time*turn);
    //   TH1D* h_wiggletot_rate = new TH1D("h_wiggletot_rate","Rate vs time",time, 0., time*turn);
    //   TGraph* wiggletot_freq = new TGraph();
    TH1D* h_gain = new TH1D("h_gain", "Gain(t)", time, 0., time*turn);
    
    
    TH1D* h_prob = new TH1D("h_prob", "", bin, 0., 1);
    
    // non perturbed histogram
    for(int it=0; it<time;it++){
        Double_t Ntot=0;
        Double_t Nt=0;
        Double_t Ntot1=0;
        Double_t Nt1=0;
        
        Double_t del = rd.Gaus(0.,1.);
        
        //Gain to be corrected......
        //gain=h_gain->GetBinContent(it);//Reading from the file for long term gain...
        //h_gain->SetBinContent(it,no_fills,gain);
        //gain1 = 1;
        gain = 1. + (eps*TMath::Exp((-it*turn)/(tau))); //exponential gain fluctutation
        
        //Theoretical gain perturbation
        for(int i=0.6*bin;i<bin;i++){//threshold is at y=0.6
            Double_t Nytot=0, Nytot1=0;
            for(int j=0; j<bin; j++){//this loop is for reconstruction considering probability distributions for all energy bins...
                y = float(j)/bin;
                y = y/(gain);
                mean = y+(1./(2*bin));
                E = mean*Emax;
                var = 0.05*sqrt(E);
                
                //var = 1e-6;
                Double_t probt=0;
                Double_t probt1=0;
                Double_t Na = hn->GetBinContent(j);
                Double_t A = ha->GetBinContent(j);
                Double_t phi = hp->GetBinContent(j);
                prob = (-TMath::Erfc((((float(i+1))/bin)-mean)/var)+TMath::Erfc(((float(i)/bin)-mean)/var))/(2.);
                probt= probt+prob;
                
                h_prob->SetBinContent(j,prob);
                Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
                Double_t Ny1 = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob1;
                Nytot = Nytot+Ny;
                //
            } //end for j
            Ntot = Ntot+Nytot;
        }//end for i
        Nt = Ntot*TMath::Exp(-(it*turn)/tau1);//Note-the value of tau in data should be fixed to 64400 and not from fitting simulation gain data?
        Nt = (Nt + del*sqrt(Nt));//statistical fluctuations
        
        h_wiggletot->SetBinContent(it,Nt);
        h_gain->SetBinContent(it, gain);
        //if(it<10) printf("INput Bin:%d Gain %f, h_pert %f, h_wiggle %f\n",it,h_gain->GetBinContent(it),h_wiggletot_pert->GetBinContent(it),h_wiggletot->GetBinContent(it));
    }//end for it
    
    fit_func1->SetLineWidth(1);
    //fit_func1->Draw();
    //TFitResultPtr r=h_wiggletot_pert->Fit("fit_func1","RS","",0,time*turn);//gives correlations b/w parameters..
    TFitResultPtr r1=h_wiggletot->Fit("fit_func1","RS","",0,time*turn);//gives correlations b/w parameters..
    //r->Print("V");
    h_wiggletot->Draw();
    gStyle->SetOptFit(1);
}



void wiggle_correct(double eps, double tau)
    {
    
        //Variables
        Double_t pi = 3.14159265358979323846;
        Double_t Emax = 3.1;
        Int_t bin = 100;
        //Int_t bin = 1000;
        Double_t eps1 = 0.001;//for one fill;
        Double_t prob;
        Double_t prob1;
        Double_t y;
        Double_t mean;
        Double_t E;
        Double_t var;
        
        Double_t fNt;
        Double_t gain = 1.;
        Double_t gain1 = 1.;
        //Double_t gain_correct=1;
        Double_t turn = 149; //considering 1 turn
        //   Double_t turn = 1; //considering every ns
        //Double_t Norm = 2.e11; //Number in Total data set
        Double_t Norm = 4.e8; // Gets crazy for 2000 events
        Double_t tau1 = 64400.;
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
        
        //N function
        TF1* Nfunc = new TF1("Nfunc","5.e9*(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);// Total data set - eq 3.19
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
        TH1D* h_wiggletot = new TH1D("h_wiggletot - Corrceted","wiggle(t)",time, 0., time*turn);
        TH1D* h_wiggletot_pert = new TH1D("h_wiggletot_pert","perturbed wiggle(t)",time, 0., time*turn);
        //   TH1D* h_wiggletot_rate = new TH1D("h_wiggletot_rate","Rate vs time",time, 0., time*turn);
        //   TGraph* wiggletot_freq = new TGraph();
        TH1D* h_gain = new TH1D("h_gain", "Gain(t)", time, 0., time*turn);
        
        
        TH1D* h_prob = new TH1D("h_prob", "", bin, 0., 1);
        
        TH1D* h_gain1 = new TH1D("h_gain1", "Gain(t)", time, 0., time*turn);
        for(int it =0; it<time;it++)
        {
            double gain = 1. + (eps*TMath::Exp((-it*turn)/(tau))); //simulated gain ()from simu fit para
            h_gain1->SetBinContent(it, gain);
            
        }
        // non perturbed histogram
        for(int it=0; it<time;it++){
            Double_t Ntot=0;
            Double_t Nt=0;
            Double_t Ntot1=0;
            Double_t Nt1=0;
            
            Double_t del = rd.Gaus(0.,1.);
            
            //Gain to be corrected......
            //gain=h_gain->GetBinContent(it);//Reading from the file for long term gain...
            //h_gain->SetBinContent(it,no_fills,gain);
            gain1 = h_gain1->GetBinContent(it);
            gain = 1. + (eps1*TMath::Exp((-it*turn)/(tau1))); //exponential gain fluctutation - theoretical
            
            //Theoretical gain perturbation
            for(int i=0.6*bin;i<bin;i++){//threshold is at y=0.6
                Double_t Nytot=0, Nytot1=0;
                for(int j=0; j<bin; j++){//this loop is for reconstruction considering probability distributions...
                    y = float(j)/bin;
                    y = y*gain1/gain;
                    mean = y+(1./(2*bin));
                    E = mean*Emax;
                    var = 0.05*sqrt(E);
                    
                    //var = 1e-6;
                    Double_t probt=0;
                    Double_t probt1=0;
                    Double_t Na = hn->GetBinContent(j);
                    Double_t A = ha->GetBinContent(j);
                    Double_t phi = hp->GetBinContent(j);
                    prob = (-TMath::Erfc((((float(i+1))/bin)-mean)/var)+TMath::Erfc(((float(i)/bin)-mean)/var))/(2.);
                    probt= probt+prob;
                    
                    h_prob->SetBinContent(j,prob);
                    Double_t Ny = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob;
                    Double_t Ny1 = Na*(1+A*TMath::Cos(0.001439*it*turn + phi*1.e-3))*prob1;
                    Nytot = Nytot+Ny;
                    //
                } //end for j
                Ntot = Ntot+Nytot;
            }//end for i
            Nt = Ntot*TMath::Exp(-(it*turn)/tau1);//Note-the value of tau in data should be fixed to 64400 and not from fitting simulation gain data
            Nt = (Nt + del*sqrt(Nt));//statistical fluctuations
            
            h_wiggletot->SetBinContent(it,Nt);
            h_gain->SetBinContent(it, gain);
            //if(it<10) printf("INput Bin:%d Gain %f, h_pert %f, h_wiggle %f\n",it,h_gain->GetBinContent(it),h_wiggletot_pert->GetBinContent(it),h_wiggletot->GetBinContent(it));
        }//end for it
        
        fit_func1->SetLineWidth(1);
        //fit_func1->Draw();
        //TFitResultPtr r=h_wiggletot_pert->Fit("fit_func1","RS","",0,time*turn);//gives correlations b/w parameters..
        TFitResultPtr r1=h_wiggletot->Fit("fit_func1","RS","",0,time*turn);//gives correlations b/w parameters..
        //r->Print("V");
        h_wiggletot->Draw();
        gStyle->SetOptFit(1);
    }
