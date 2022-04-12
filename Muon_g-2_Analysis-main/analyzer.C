#include "TSpline.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TStyle.h"

Int_t channelposition(int channel){

  Int_t posch=-1;

  if(channel==0) posch=1;
  if(channel==1) posch=10;
  if(channel==2) posch=19;
  if(channel==3) posch=28;
  if(channel==4) posch=37;

  if(channel==5) posch=2;
  if(channel==6) posch=11;
  if(channel==7) posch=20;
  if(channel==8) posch=29;
  if(channel==9) posch=38;

  if(channel==10) posch=3;
  if(channel==11) posch=12;
  if(channel==12) posch=21;
  if(channel==13) posch=30;
  if(channel==14) posch=39;


  if(channel==15) posch=4;
  if(channel==16) posch=13;
  if(channel==17) posch=22;
  if(channel==18) posch=31;
  if(channel==19) posch=40;

  if(channel==20) posch=5;
  if(channel==21) posch=14;
  if(channel==22) posch=23;
  if(channel==23) posch=32;
  if(channel==24) posch=41;

  if(channel==25) posch=6;
  if(channel==26) posch=15;
  if(channel==27) posch=24;
  if(channel==28) posch=33;
  if(channel==29) posch=42;



  if(channel==30) posch=7;
  if(channel==31) posch=16;
  if(channel==32) posch=25;
  if(channel==33) posch=34;
  if(channel==34) posch=43;

  if(channel==35) posch=8;
  if(channel==36) posch=17;
  if(channel==37) posch=26;
  if(channel==38) posch=35;
  if(channel==39) posch=44;

  if(channel==40) posch=9;
  if(channel==41) posch=18;
  if(channel==42) posch=27;
  if(channel==43) posch=36;
  if(channel==44) posch=45;

  return posch;
  
  
}


void analyzer() {
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  // open data file
  
  //TFile *file = new TFile("run5269/run05269_08.root");

  TFile *file = new TFile("gm2_run01270_35.root");

  TTree *ch[60];
  TString channel;
  TString histo;
  TH2F *h[60];
  vector<short> *trace = 0;

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->Divide(9,5);
  Int_t thech=-1;
  //TCanvas *c2 = new TCanvas("c2","c2",600,400);
  //c2->Divide(6,5);
  
  gStyle->SetPalette(1);


  for (Int_t  nch =0; nch < 45; nch++){

    
    
    channel.Form("IslandDumper/calo25/xtal%d", nch);
   
    ch[nch] = (TTree*)file->Get(channel.Data());

    trace = 0;
    ch[nch]->SetBranchAddress("trace", &trace);

    ch[nch]->GetEntry(0);
  
    histo.Form("h%d", nch);


    int histobins= trace->size();
    if(trace->size()<0) histobins=1000;
    
    h[nch]= new TH2F(histo.Data(),"",histobins,0,histobins,5000,-2500,2500);

    cout << "channel "<< nch << " " << trace->size()<< " " << ch[nch]->GetEntries() << endl;
    
   
    // for (Int_t  k =0; k < ch[nch]->GetEntries() ; k++){

    for (Int_t  k =0; k < 100 ; k++){


      ch[nch]->GetEntry(ch[nch]->GetEntries()-(k+1));

	if(trace->size()<0) continue;
	
      for(size_t i=0; i<trace->size(); i++){
	h[nch]->Fill(i,trace->at(i));
      }

      //      cout << channelposition(nch) <<endl;
      
      // c1->cd(channelposition(nch));

      thech=channelposition(nch);
      c1->cd(thech);
      h[nch]->DrawCopy("COLZ");	
	
      }
  }

  
  file->Close();
  
}

