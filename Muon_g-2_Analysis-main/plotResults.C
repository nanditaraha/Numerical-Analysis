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
  TCanvas *c0 = new TCanvas("c1","multigraph L3",200,10,700,500);
  //c0->SetFrameFillColor(30);
  TMultiGraph *mg = new TMultiGraph();
  /*
//These are the previous results - gain for a fill only
  double eps[7]={-0.003,-0.002,-0.001,0.,0.001,0.002,0.003};
  double linear[7]={35.679, 23.636, 11.744, 0.000000, -11.597,-23.041,-34.338};
  double exp[7]={1202.465,799.587,398.761,0., -396.686,-791.281,-1183.776};
  double cos[7]= {-106.742, -69.965, -34.367, 0.0, 33.188, 65.226,96.137};
  double exp_cos[]={1110.291,735.165,365.092,0.0,-360.223,-715.525,-1066.087};
*/
    //These are results for the long gain
    
    double eps[7]={-0.003,-0.002,-0.001,0.,0.001,0.002,0.003};
    double linear[7]={-29.63,-235.21,1166.29,73.83,521.24,7.89,-258.86};
    double exp[7]={144.02,940.75,777.33,73.84,92.70,-366.76,-653.57};
    double cos[7]= {381.24,1321.74,139.49,73.84,1298.79,-145.9,245.26};
    double exp_cos[7]={1119.29,11.13,576.13,73.84,493.91,-230.86,-662.53};
    
  TGraph *gr1 = new TGraph(7,eps,linear); gr1->SetLineColor(kBlue);
  TGraph *gr2 = new TGraph(7,eps,exp); gr2->SetLineColor(kRed);
  TGraph *gr3 = new TGraph(7,eps,cos); gr3->SetLineColor(kGreen);
  TGraph *gr4 = new TGraph(7,eps,exp_cos); gr4->SetLineColor(kOrange);
  
  mg->Add(gr4); gr1->SetTitle("Linear"); gr4->SetLineWidth(3);
  mg->Add(gr3); gr2->SetTitle("Exp")  ; gr3->SetLineWidth(3);
  mg->Add(gr2); gr3->SetTitle("Cos(x)")    ; gr2->SetLineWidth(3);
  mg->Add(gr1); gr4->SetTitle("Exp and Cos(x)")  ; gr1->SetLineWidth(3);
  
  gr4->GetXaxis()->SetTitle("#epsilon");
  mg->Draw("a");
  gPad->SetGridx();
  gPad->SetGridy();
  c0->BuildLegend();

}

