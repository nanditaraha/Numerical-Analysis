{
  TFile *f1 = new TFile("analyzed_gm2slac_run01896.root");
  TTree *t1 = (TTree*)f1->Get("t1");

  TFile *f2 = new TFile("Run23_20160604072759.002.root");
  TTree *t2 = (TTree*)f2->Get("t1");
  Float_t px, py, pz;
  Double_t random;
  Int_t ev;
  t1->SetBranchAddress("px",&px);
  t1->SetBranchAddress("py",&py);
  t1->SetBranchAddress("pz",&pz);
  t1->SetBranchAddress("random",&random);
  t1->SetBranchAddress("ev",&ev);
  
  //create two histograms
  TH1F *hpx   = new TH1F("hpx","px distribution",100,-3,3);
  TH2F *hpxpy = new TH2F("hpxpy","py vs px",30,-3,3,30,-3,3);

  //read all entries and fill the histograms
  Long64_t nentries = t1->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    t1->GetEntry(i);
    hpx->Fill(px);
    hpxpy->Fill(px,py);
  }
}
