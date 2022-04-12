//Experimenting different gains triangle
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPolyLine.h"
#include "TCanvas.h"//to stop resize and assign sizes
#include "TVirtualPad.h" //for log scale in of pads gPad->SetGridy()

TCanvas *c;
TGraph *gr; 
//TMarker *r1, *r2,*r3,*r4,*r5,*r6, *r1_1, *r2_1,*r3_1,*r4_1,*r5_1,*r6_1,*r1_2, *r2_2,*r3_2,*r4_2,*r5_2,*r6_2, *r7,*r7_1;

void plot_Results()
{
  //TCanvas *c0 = new TCanvas("c1","multigraph L3",200,10,700,500);
    
    double eps[10]={0,10,20,40,60,100,200,300,400,500};
    double exp1[10]={-0.028,-0.021908,-0.014880,-0.009486,-0.004494,-0.000001,0.001006,0.000368,0.000107,0.000032};
  TGraph *gr1 = new TGraph(10,eps,exp1); gr1->SetLineColor(kBlue);
  
  
  gr1->GetXaxis()->SetTitle("#epsilon");
  gr1->Draw("ac*");
  gPad->SetGridx();
  gPad->SetGridy();
  //c0->BuildLegend();

}

