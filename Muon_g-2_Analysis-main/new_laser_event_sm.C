#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

void new_laser_event_sm()
{
  //TString fileName("july_1_4.root");//all files
  //TString fileName("run01750.root"); 
  
  TString fileName("runs_1791_ro_1869_notrace.root");//New runs with timestamp..
  gStyle->SetPalette(1);

  TH2D *hPedestal = new TH2D("hPedestal", " ",202,0,101,200, 1780, 1790);//in fill
  TH2D *hArea = new TH2D("hArea", " ",202,0,101,2000, 20000, 40000);//in fill  
  TH2D *hArea_ratio = new TH2D("hArea_ratio", " ",202,0,101,100, 0.42, 0.82);//in fill  
  TH2D *hTime_diff = new TH2D("hTime_diff", "",202,0,101,10, 230, 240);//in of fill 
  //TH2D *hIsland = new TH2D("hIsland", " ",202,0,101,50, 0, 50);//in fill

  TFile * file = TFile::Open(fileName,"READ");
  if (!file) { return; }
  TTree *tree;
  file->GetObject("italianoTree/monitors",tree);
  //file->GetObject("italianoTreeOoF/monitors",tree);

  Int_t nentries = (Int_t)tree->GetEntries();
  Double_t pedestal;
  UInt_t channelNum;
  Double_t area;
  Double_t amplitude;
  Double_t area_second_pulse;
  Double_t time_second_pulse;
  Double_t time;
  int eventNum;
  int runNum;
  int islandNum;
  int runStartUnixTimeSeconds;
  
  int time_init;
  int time_final;
  
  double run_hour;

  int first_run=1791;
  int last_run=1869;

  TBranch *bEvt = tree->GetBranch("runStartUnixTimeSeconds");
  bEvt->SetAddress(&runStartUnixTimeSeconds); 

  TBranch *bRun = tree->GetBranch("runNum");
  bRun->SetAddress(&runNum); 
  TBranch *bPed = tree->GetBranch("pedestal");
  bPed->SetAddress(&pedestal);  
  TBranch *bChannel = tree->GetBranch("channelNum");
  bChannel->SetAddress(&channelNum);
  TBranch *bAmp = tree->GetBranch("amplitude");
  bAmp->SetAddress(&amplitude);
  TBranch *bArea = tree->GetBranch("area");
  bArea->SetAddress(&area);
  TBranch *bArea2 = tree->GetBranch("area_second_pulse");
  bArea2->SetAddress(&area_second_pulse);
  TBranch *bTime2 = tree->GetBranch("time_second_pulse");
  bTime2->SetAddress(&time_second_pulse);
  
  TBranch *bTime = tree->GetBranch("time");
  bTime->SetAddress(&time);

  int max=0;

  tree->GetEntry(0);
  time_init=runStartUnixTimeSeconds;
  //cout<<"Tree entry:"<<nentries<<"\n";
  for (Int_t i = 0; i < nentries; i++) {
    //for (Int_t i = 300000; i < nentries; i++) {
    tree->GetEntry(i);
    //if(channelNum==0 && amplitude>500) {
    if(amplitude>500) {
      double amp1,amp2;
      if(channelNum==33) amp1=amplitude;
      if(channelNum==34) amp2=amplitude;
      hPedestal->Fill((runStartUnixTimeSeconds - time_init -3228)/3600.,pedestal);
      hArea->Fill((runStartUnixTimeSeconds - time_init -3228)/3600.,area_second_pulse);
      //hIsland->Fill(,islandNum);
      //hArea_ratio->Fill((runStartUnixTimeSeconds - time_init -3228)/3600.,area_second_pulse/area);
      hArea_ratio->Fill((runStartUnixTimeSeconds - time_init -3228)/3600.,amp2/amp1);
      //if(i<100)
      //cout<<"time"<<(time_second_pulse)*1.25<<"\n";
      hTime_diff->Fill((runStartUnixTimeSeconds - time_init -3228)/3600.,(time_second_pulse-time)*1.25);
      
    }
     //if(i<100) cout<<i+1<<": Run Number:"<<runNum<<" Event:"<<Event[i]<<" Event #:"<<eventNum<<"\n";
   }

  //for(int i=0;i<10;i++) cout<<i+1<<": Run Number:"<<runNum<<" Event:"<<Event[i]<<"Event size:"<<Event.size()<<"\n";
   //hEvt->Draw("colz");

   hPedestal->GetXaxis()->SetTitle("Time (hr)");
   hPedestal->GetYaxis()->SetTitle("Pedestal (ADC)");
   hPedestal->SetStats(0);

   hArea->GetXaxis()->SetTitle("Time (hr)");
   hArea->GetYaxis()->SetTitle("Area (ADC)");
   hArea->SetStats(0);

   hArea_ratio->GetXaxis()->SetTitle("Time (hr)");
   hArea_ratio->GetYaxis()->SetTitle("Ratio of Areas");
   hArea_ratio->SetStats(0);

   hTime_diff->GetXaxis()->SetTitle("Time (hr)");
   hTime_diff->GetYaxis()->SetTitle("Time diff (t_{2}-t_{1}) ns");
   hTime_diff->SetStats(0);

   TProfile *p = hPedestal->ProfileX();
   p->GetXaxis()->SetTitle("Time (hr)");
   p->GetYaxis()->SetTitle("Pedestal (ADC)");
   p->GetYaxis()->SetRangeUser(1785,1787);//for ped
   p->SetStats(0);

   TProfile *a = hArea->ProfileX();
   a->GetXaxis()->SetTitle("Time (hr)");
   a->GetYaxis()->SetTitle("Area (ADC)");
   a->GetYaxis()->SetRangeUser(2260,2350);//for ped
   a->SetStats(0);

   TProfile *ar = hArea_ratio->ProfileX();
   ar->GetXaxis()->SetTitle("Time (hr)");
   ar->GetYaxis()->SetTitle("Ratio of Area");
   ar->GetYaxis()->SetRangeUser(0.95,1.05);//for ped
   ar->SetStats(0);
   //hAmp_ratio->Scale(1/ar->GetMean(2));   

   TProfile *ar_y = hArea_ratio->ProfileY();
   //ar_y->Draw();

   TProfile *t = hTime_diff->ProfileX();
   t->GetXaxis()->SetTitle("Time (hr)");
   t->GetYaxis()->SetTitle("Time diff (t_{2}-t_{1}) ns");
   t->GetYaxis()->SetRangeUser(235.8,236.4);//for ped
   t->SetStats(0);
  
   //Histograms: 
   //gPad->SetGridx();
   //hPedestal->Draw("colz");
   //hArea->Draw("colz");
   //hArea_ratio->Draw("colz");
   //hTime_diff->Draw("colz");
   //hIsland->Draw("colz");

   //Profiles:
   //gPad->SetGridx();
   //p->Draw();
   //a->Draw();

   //Scaling ar
   ar->Scale(1/ar->GetMean(2));
   ar->Draw();
   //t->Draw();
   
   TFile *f1 = new TFile("infill_sm.root","RECREATE");
   //TFile *f1 = new TFile("outfill.root","RECREATE");
   ar->Write("prof_ar");
   t->Write("prof_t");
   p->Write("prof_p");
   hArea->Write("area");
   hArea_ratio->Write("hist_area_ratio");
   hPedestal->Write("hist_ped");
   hTime_diff->Write("hist_time");
   
}
