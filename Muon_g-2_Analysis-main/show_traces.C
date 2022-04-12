{

  TCanvas *c = new TCanvas();
  c->Divide(4,6);

  TString condition;
  //TString entry;
  TString g_name;
  //TGraph *gr[24];
  TFile *f = TFile::Open("runs_1791_ro_1869.root");//all runs
  //TFile *f = TFile::Open("calib_all.root");//calib runs from 2053 to 2059
  italianoTree->cd();
  
  //Int_t n1 = monitors->Draw("trace:Iteration$","Entry$==90 &&channelNum==0 && amplitude>500","goff");
  //int n1=monitors->GetSelectedRows();
  //Double_t *v1 = monitors->GetVal(1);
  //Double_t *v2 = monitors->GetVal(0);

  //TGraph *gr = new TGraph(n1, v1, v2);
  //gr->Draw("ap")
  
  for(int i=0;i<24;i++){
    
    g_name.Form("Channel %d",i+1);
    condition.Form("Entry$==%d &&channelNum==%d && amplitude>500",i+90,i);
    Int_t n1 = monitors->Draw("trace:Iteration$",condition.Data(),"goff");
    //int n1=monitors->GetSelectedRows();
    
    Double_t *v1 = monitors->GetVal(1);
    Double_t *v2 = monitors->GetVal(0);

    TGraph *gr = new TGraph(n1, v1, v2);
    gr->SetLineColor(kBlue);
    gr->SetTitle(g_name.Data());
    gr->GetXaxis()->SetRangeUser(0,600);    
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.07);
    c->cd(i+1);
    gr->Draw("alp");

  }

  /*
c->cd(1);
  gr->Draw("ap");
  */
}
