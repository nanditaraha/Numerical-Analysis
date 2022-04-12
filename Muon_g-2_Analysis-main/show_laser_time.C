//Never name functions if they are commond line interpretors - no declarations required
//void show_laser()
{
  //TFile *_file0 = TFile::Open("tree_1270.root");
TFile *_file0 = TFile::Open("runs_1791_ro_1869_notrace.root");
  TGraphErrors *g_amp, *g_chi2;
  TString channel;
  
  italianoTree->cd();
  gStyle->SetPalette(1);
  double mean1[24];
  double chi2_1[24];
  double m_err1[24];
  double c_err1[24];
  double ch[24];
  double ch_e[24]={0};

  for(int i=0;i<24;i++){
    channel.Form("channelNum==%d",i);
    //monitors->Draw("pedestal>>h1(1000,1700,1800)",channel.Data(),"goff");//to not draw 
    monitors->Draw("(time_second_pulse-time)*1.25>>h1",channel.Data(),"goff");//to not draw 
    //monitors->Draw("amplitude>>h1(1000,0,3500)",channel.Data(),"goff");//to not draw 
    TH1F *h1 = (TH1F*)gDirectory->Get("h1");
    //Above step is very important to retrieve the histogram or  the code would break in the next element of the loop;
    mean1[i] =h1->GetMean();
    m_err1[i]=h1->GetRMS();
    ch[i]=i+1;
    //cout<<"Mean "<<mean1[i]<<"\n";
    /*
    if(i==20) {
    h1->Draw();
    break;
    }
    */
  }
  //Amplitudes
  
  g_amp = new TGraphErrors(24,ch,mean1,ch_e,m_err1);
  g_amp->SetTitle("Time difference b/w first and second pulse all LM");
  g_amp->GetXaxis()->SetTitle("Calo #");
  //g_amp->GetYaxis()->SetTitle("Amplitude ratio (A2:A1)");
  g_amp->GetYaxis()->SetTitle("Time difference (T2-T1) ");
  g_amp->SetMarkerColor(kBlue);
  g_amp->SetMarkerStyle(kFullTriangleUp);
  g_amp->Draw("AP");

}
