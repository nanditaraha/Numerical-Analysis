/*
 My undestanding of this: Vector decayTimes (private member - decay drop for a decay. vector of n pulses) and doubles dropPerDecay and recoveryTime is for each pulse.
 An instance of gainDropInFill (drop here) calls operator() which is like a root function.
 
 Laser: 
 drop
Function gainSimulation() draws the function for 100 muons of a fill (just one fill)
  */
const int cycle=2000;//need 2000 trying with less first...
const int fills=16;//need 16 trying with less first...
const int laser=8;
const int nMuon=100;

double dropE = 6.4e-4;
double tau_r=18;
int gap = 160; //This is the gap in time in us arising from the frequency of the laser
int offset = 5; //This is the offset we move the laser from fill to fill in us
double delta_t = 0.5; //+/- Time range about a laser shot in us 
class gainDropInFill {
 public:
  gainDropInFill(vector<double> decayTimesIn, vector<double> dropPerDecayIn,
                double recoveryTimeIn)
      : decayTimes(decayTimesIn),
        dropPerDecay(dropPerDecayIn),
        recoveryTime(recoveryTimeIn) {
    sort(decayTimes.begin(), decayTimes.end());
      }

  double operator()(double* x, double* p) {
    double t = x[0];
    double gain = 1;
    if (t > decayTimes[0]) {
      // drop gain for the first pulse
      gain = 1 - dropPerDecay[0];
      unsigned int i = 1;
      for (; i < decayTimes.size() && decayTimes[i] < t; ++i) {
        // recover gain from pulse i - 1 to pulse i
        gain = 1 -
               (1 - gain) *
                   exp(-(decayTimes[i] - decayTimes[i - 1]) / recoveryTime);
        // drop gain for pulse i
        gain = gain * (1 - dropPerDecay[i]);
      }
      // recover gain from the last pulse before time t to time t
      gain = 1 - (1 - gain) * exp(-(t - decayTimes[i - 1]) / recoveryTime);
    }
    return gain;
  }

 private:
  vector<double> decayTimes;
  // decay drop for one decay (when gain is fully recovered)
  vector<double> dropPerDecay;
  // bias reconvery time constant for one pulse
  double recoveryTime;
};

void energy_dist(int cycle=2000, int fill=16){
  //TFile *f = new TFile::Open("prob.root");
  TFile *f = TFile::Open("prob.root");
  TH1 *Prob = (TH1*)f->Get("Prob");
  TH1D *h = new TH1D("h", "Random energy dist", 3100,0,3100.);
  for(int c=0;c<cycle;c++){
    for(int n=0;n<fill;n++){
      h->Fill(Prob->GetRandom());
    }
  }
  h->Draw();
  Prob->SetLineColor(kRed);
  Prob->SetLineWidth(2);
  //Prob->Scale(h->Integral()/Prob->Integral());
  Prob->Scale(h->GetMaximum()/(2*Prob->GetMaximum()));
  Prob->Draw("same");
}

void gainSimulationTest(int cycle,int fill, int laser) {
  //gRandom->SetSeed(0);
  double mean,sigma,data,err,data1;
  TGraph* g = new TGraph();
  TCanvas* c1 = new TCanvas();
  TF1* f;
  TF1* f1;
  // TF1* f_gauss = new TF1("f_gauss","gaus",0,2);
  TH1D* hf; TH1D* hf1; TH1D* hf2; TH1D *hl;
  TH2D* xhist1 = new TH2D("xhist1","SiPM gain Vs. time - Muons only", 7e2,0,700,50000,0.95,1.0);//final gain for all 2000 cycles
  TH2D* xhist2 = new TH2D("xhist2","SiPM gain Vs. time - Laser and Muons", 7e2,0,700,50000,0.95,1.0);//final gain for all 2000 cycles
  TH2D* xhist3 = new TH2D("xhist3","SiPM gain Vs. time - Laser only", 7e2,0,700,50000,0.95,1.0);//final gain for all 2000 cycles
  TH1D* hist1 = new TH1D("hist1","SiPM gain Vs. time - Laser only", 7e2,0,700);//final gain for all 2000 cycles
  TH1D* hist2 = new TH1D("hist2","SiPM gain Vs. time - All pulses", 7e2,0,700);//final gain for all 2000 cycles
  TRandom3 *random = new TRandom3;
  int nPulses = nMuon;
  TProfile *py;
  bool laser_on = false;
  //vector<double> times;
  //vector<double> dropEnergy;
  //vector<double> times_noLaser;
  //vector<double> dropEnergy_noLaser;
  //vector <vector <vector<double> > >  times_Laser;  //xx(cycles, vector <vector<double> >(fill,extraction_mu()));
  //vector <vector <vector<double> > >  dropEnergy_Laser;   
  vector<double> times_Laser;
  vector<double> dropEnergy_Laser;

  int size=0, sizeL=0;
  for(int c=0;c<cycle;c++){
    for(int n=0;n<fill;n++){
      vector<double> decayTimes;
      vector<double> dropPerDecay;
      vector<double> decayTimes_noLaser;
      vector<double> dropPerDecay_noLaser;
      vector<double> decayTimes_Laser;
      vector<double> dropPerDecay_Laser;
      gRandom->SetSeed(0);
      for (int i = 0; i < nPulses; ++i){ 
	decayTimes.push_back(gRandom->Exp(64));
	dropPerDecay.push_back(dropE);
	decayTimes_noLaser.push_back(gRandom->Exp(64));
	dropPerDecay_noLaser.push_back(dropE);
      }
      for (int i = 0; i < laser; ++i){
	decayTimes.push_back(gap*i+offset*n);
	//dropPerDecay.push_back(dropE*10);
	dropPerDecay.push_back(1e-2);
	//decayTimes_Laser.push_back(gap*i+offset*n);
	//dropPerDecay_Laser.push_back(dropE*10);
	laser_on=true;
      }      	
    
      gainDropInFill drop(decayTimes, dropPerDecay, tau_r);//this is the function call to plot SiPM behaviour - Muons + laser....
      //gainDropInFill drop(decayTimes_Laser, dropPerDecay_Laser, tau_r);//this is the function call to plot SiPM behaviour - Muons + laser....
      gainDropInFill drop_noLaser(decayTimes_noLaser, dropPerDecay_noLaser, tau_r);//this is the function call to plot SiPM behaviour..
      //gainDropInFill drop_Laser(decayTimes_Laser, dropPerDecay_Laser, tau_r);//this is the function call to plot SiPM behaviour..

      f = new TF1(Form("%i pulses", nPulses+laser), drop, -10, 700, 0,
		    "gainDropInFill");
      f->SetTitle(Form("%i pulses; time [#mus]; gain", nPulses+laser));
      f1 = new TF1(Form("%i pulses", nPulses), drop_noLaser, -10, 700, 0,
		  "gainDropInFill");
      f1->SetTitle(Form("%i pulses; time [#mus]; gain", nPulses));
      
      //PLot graph only for 1 cycle 1 fill to test
      /*
      TGraph* g = new TGraph();
      g->SetName(Form("g%i", nPulses));
      for (unsigned int i = 0; i < decayTimes.size(); ++i) {
	g->SetPoint(i, decayTimes[i], f->Eval(decayTimes[i]));
	cout<<"Samples:"<<i<<"="<<decayTimes[i]<<" Gain="<<f->Eval(decayTimes[i])<<" E drop="<<dropPerDecay[i]<<"\n";
      }
      g->SetMarkerStyle(kFullStar);
      g->SetMarkerSize(2);
      g->SetMarkerColor(kBlue);
      f->SetNpx(1000);
      f->Draw();
      g->Draw("p");
      */
      for (unsigned int i = 0; i < decayTimes.size(); ++i){//size should be nPulses+laser = Muon + Laser
	//cout<<"Samples:"<<i<<"="<<decayTimes[i]<<" Gain="<<f->Eval(decayTimes[i])<<" E drop="<<dropPerDecay[i]<<"\n";
	if(i>=nPulses){
	  //
	  /*
	  double time = decayTimes[i];
	  for (unsigned int j = 0; j < decayTimes.size(); ++j)
	    if(time>decayTimes[i]-delta_t && time<decayTimes[i]+delta_t){
	      cout<<"Samples:"<<j<<"="<<decayTimes[j]<<" Gain="<<f->Eval(decayTimes[j])<<" E drop="<<dropPerDecay[j]<<"\n";
	    }
	  */
	  xhist3->Fill(decayTimes[i],f1->Eval(decayTimes[i])); 
	  hist1->SetBinContent(hist1->FindBin(decayTimes[i])-1,f1->Eval(decayTimes[i])); 
	}
	xhist2->Fill(decayTimes[i], f->Eval(decayTimes[i]));
	hist2->SetBinContent(hist2->FindBin(decayTimes[i])-1,f1->Eval(decayTimes[i])); 
      }

      for (unsigned int i = 0; i < decayTimes_noLaser.size(); ++i){ //Moun only 
	xhist1->Fill(decayTimes_noLaser[i],f1->Eval(decayTimes_noLaser[i]));
      }
    
      f->SetNpx(1000);
      f1->SetNpx(1000);
      // f2->SetNpx(1000);
      //f->Draw();
            
      //printf("=============End of fill %d===============\n",n+1);
      //c->Print(Form("gainDrop%iPulses.pdf", nPulses));
    }
    //printf("=============End of cycle %d===============\n",c+1);    
    
  }
  hf = new TH1D("Muons only","Gain versus time - Laser Shots",7e2,0,700);
  hf1 = new TH1D("Laser","Gain versus time - Muons Only",7e2,0,700);
  hf2 = new TH1D("All","Gain versus time - All pulses ",7e2,0,700);
    
  for ( Int_t i=0; i<700; i++){   //Muons only
    xhist1->GetXaxis()->SetRange(i,(i+1));
    double time = xhist1->GetMean(1);
    double gain = xhist1->GetMean(2);
    double gain_err = xhist1->GetRMS(2)/sqrt(xhist1->Integral());
    hf1->SetBinContent(hf1->FindBin(time)-1, gain); // Bin starts from 1 and arrays from 0
    hf1->SetBinError(hf1->FindBin(time)-1,gain_err);
  }
  
  for ( Int_t i=0; i<700; i++){   //Muons + Laser  
    xhist2->GetXaxis()->SetRange(i,(i+1));
    double time = xhist2->GetMean(1);
    double gain = xhist2->GetMean(2);
    double gain_err = xhist2->GetRMS(2)/sqrt(xhist2->Integral());
    hf2->SetBinContent(hf2->FindBin(time)-1, gain); // Bin starts from 1 and arrays from 0  
    hf2->SetBinError(hf2->FindBin(time)-1,gain_err);
  }
  
  /*
  //Trying with a 1-d histogram - it is wrong!!!! 1-d plots just the last pulse
  for ( Int_t i=0; i<700; i++){   //Muons + Laser                                               
    hist2->GetXaxis()->SetRange(i,(i+1));                                                   
    double time = hist2->GetMean(1);                                       
    double gain = hist2->GetMean(2);                                      
    double gain_err = hist2->GetRMS(2)/sqrt(hist2->Integral());                                            
    hf2->SetBinContent(hf2->FindBin(time)-1, gain); // Bin starts from 1 and arrays from 0                                       
    hf2->SetBinError(hf2->FindBin(time)-1,gain_err);                                                       
  }  

  for ( Int_t i=0; i<(int)(700/offset); i++){   //Laser only 
    hist1->GetXaxis()->SetRange(i*offset+3,(i+1)*offset+2);                                                                      
    double time = offset*i;                                                                                           
    double gain = hist1->GetBinContent(hist1->FindBin(time));            
    //double gain_err = hist1->GetRMS(2)/sqrt(hist1->Integral());                                                               
    hf->SetBinContent(hf->FindBin(time)-1, gain);                                                                         
    //hf->SetBinError(hf->FindBin(time)-1,gain_err);                                                                                              
  }  
  */
  for ( Int_t i=0; i<(int)(700/offset); i++){   //Laser only                                                                       
    xhist3->GetXaxis()->SetRange(i*offset+3,(i+1)*offset+2);
    double time = offset*i;
    double gain = xhist3->GetMean(2);
    double gain_err = xhist3->GetRMS(2)/sqrt(xhist3->Integral());
    hf->SetBinContent(hf->FindBin(time)-1, gain);
    hf->SetBinError(hf->FindBin(time)-1,gain_err);
  }
  
  /*
  for ( Int_t i=0; i<(int)(700/offset); i++){   //Laser + Muons  Avg. for laser points      
    xhist2->GetXaxis()->SetRangeUser(i*offset-0.5,i*offset+0.5);
    double time = offset*i; 
    double gain = xhist2->GetMean(2);                                                          
    double gain_err = xhist2->GetRMS(2)/sqrt(xhist2->Integral());  
    hf->SetBinContent(hf->FindBin(time)-1, gain); 
    hf->SetBinError(hf->FindBin(time)-1,gain_err);                                                                                              
  } 
  */
  hf1->GetYaxis()->SetRangeUser(0.988,1.002);
  hf1->GetXaxis()->SetRangeUser(0,300);
  hf1->GetYaxis()->SetTitle("#frac{G}{G_{0}}   ");
  hf1->GetXaxis()->SetTitle("Time (#mu s)");
  hf1->SetStats(0);
  gPad->SetGridx();  
  gPad->SetGridy(); 
  hf1->GetYaxis()->SetRangeUser(0.987,1.002);
  hf1->GetXaxis()->SetRangeUser(0,300);
  hf1->Draw("e");
  
  hf->SetMarkerStyle(kFullStar);
  hf->SetMarkerSize(1);
  hf->SetMarkerColor(kBlue);
  hf->Draw("same");

  hf2->GetYaxis()->SetRangeUser(0.988,1.002);
  hf2->GetXaxis()->SetRangeUser(0,300);
  hf2->GetYaxis()->SetTitle("#frac{G}{G_{0}}   ");
  hf2->GetXaxis()->SetTitle("Time (#mu s)");
  hf2->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hf2->GetYaxis()->SetRangeUser(0.987,1.002);
  hf2->GetXaxis()->SetRangeUser(0,300);
  //hf2->Draw("e");

  TLatex latex;
  latex.SetTextSize(0.037);
  latex.DrawLatex(121.5,.995,Form("Laser %d #mus, %d #mus step",gap,offset));
  //latex.DrawLatex(121.5,.995,"No Laser - 100 mouns each fill");
  latex.DrawLatex(121.5,.994,Form("%d cycles, %d fills each, %d laser shot",cycle,fill,laser));
  latex.DrawLatex(121.5,.992,Form("p_{0} = P_{0} #tau_{#mu} = %0.1f #times 10^{-4}",dropE*1e4));
  
  latex.SetTextColor(kRed);
  latex.DrawLatex(121.5,.993,"#frac{G}{G_{0}} = 1 - #frac{P_{0} n_{0} #tau_{#mu}#tau_{r}}{#tau_{#mu}-#tau_{r}}#left(e^{-t/#tau_{#mu}} - e^{-t/#tau_{r}} #right)");
  TF1 *f3 = new TF1("f3","1-[0]*[1]*64.4*[2]*(exp(-(x/64.4))- exp(-(x/[2])))/(64.4-[2])",0.,160.);
  f3->SetParameters(1e-5,nMuon,tau_r);
  f3->Draw("same");
  
  TFile *file1 = TFile::Open(Form("All_%d_%d.root",gap,offset), "RECREATE");
  //TFile *file1 = TFile::Open(Form("All_%d_%d.root",gap,offset), "RECREATE");
  hf->Write("Laser");
  hf1->Write("Muons");
  hf2->Write("All");
  //hist1->Write("Laser1");
  //hist2->Write("All1");
}

