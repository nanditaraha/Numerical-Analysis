{

 TF1* Nfunc = new TF1("Nfunc","(1./3.)*(x-1.)*(4.*x*x-5*x-5.)*(3./5.)*0.008536",0,1);
 Nfunc->SetNpx(10000000);
 /*TH1D* h1_n = new TH1D("h1_n","N(y)",100,0.,1.);
 h1_n->FillRandom("Nfunc",2e9);
 
 TCanvas* c1 = new TCanvas("c1");
 c1->Divide(1,2);
 c1->cd(1);
 h1_n->Draw();
 c1->cd(2);
 */
 Nfunc->Draw();
}

