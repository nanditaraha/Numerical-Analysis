class gainDropInFill {
 public:
  gainDropInFill(vector<double> decayTimesIn, double dropPerDecayIn,
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
      gain = 1 - dropPerDecay;
      unsigned int i = 1;
      for (; i < decayTimes.size() && decayTimes[i] < t; ++i) {
        // recover gain from pulse i - 1 to pulse i
        gain = 1 -
               (1 - gain) *
                   exp(-(decayTimes[i] - decayTimes[i - 1]) / recoveryTime);
        // drop gain for pulse i
        gain = gain * (1 - dropPerDecay);
      }
      // recover gain from the last pulse before time t to time t
      gain = 1 - (1 - gain) * exp(-(t - decayTimes[i - 1]) / recoveryTime);
    }
    return gain;
  }

 private:
  vector<double> decayTimes;
  // decay drop for one decay (when gain is fully recovered)
  double dropPerDecay;
  // bias reconvery time constant for one pulse
  double recoveryTime;
};

void gainSimulation() {
  gRandom->SetSeed(0);
  //int nMuonsArray[] = {10, 20, 50, 100, 300, 1000};
  int nMuonsArray[] = {100};
  int nMuonsArraySize = distance(begin(nMuonsArray), end(nMuonsArray));
  for (int j = 0; j < nMuonsArraySize; ++j) {
    TCanvas* c = new TCanvas();
    int nPulses = nMuonsArray[j];
    vector<double> decayTimes;
    for (int i = 0; i < nPulses; ++i) {
      decayTimes.push_back(gRandom->Exp(64));
    }
    gainDropInFill drop(decayTimes, 0.001, 3.0);
    TF1* f = new TF1(Form("%i pulses", nPulses), drop, -10, 140, 0,
                     "gainDropInFill");
    f->SetTitle(Form("%i pulses; time [#mus]; gain", nPulses));
    
    TGraph* g = new TGraph();
    g->SetName(Form("g%i", nPulses));
    for (unsigned int i = 0; i < decayTimes.size(); ++i) {
      g->SetPoint(i, decayTimes[i], f->Eval(decayTimes[i]));
    }
    g->SetMarkerStyle(kFullStar);
    g->SetMarkerSize(2);
    g->SetMarkerColor(kBlue);
    
    f->SetNpx(1000);
    f->Draw();
    g->Draw("p");
    //c->Print(Form("gainDrop%iPulses.pdf", nPulses));
  }//end of nMuonArray size
}
