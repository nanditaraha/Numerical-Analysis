#include <TH1D.h>
#include <TFile.h>

void norm_gains(int kink1=398, int kink2=1900)
{
  Double_t x1, x2, ratio;
  TFile *nwig = TFile::Open("nowiggle_gain.root");
  TFile *wig = TFile::Open("gain.root");

  TH1D *h1 = (TH1D*)wig->Get("Gain");
  TH1D *h2 = (TH1D*)nwig->Get("Gain");

  x1=h1->GetMinimum();
  x2=h2->GetMinimum();
  ratio = x1/x2;

  for (int i=1;i<kink1;i++)
    {
      //x1 = h1->GetBinContent(i);
      //x2 = h2->GetBinContent(i);
      h2->SetBinContent(i,h2->GetBinContent(i)*ratio);
      
    }
  
 
 for (int i=kink1;i<kink2;i++)
   {
     //x1 = h1->GetBinContent(i);
     //x2 = h2->GetBinContent(i);
      h2->SetBinContent(i,h1->GetBinContent(i));
      
   }
 
  h1->Draw();
  h2->Draw("same");

}
