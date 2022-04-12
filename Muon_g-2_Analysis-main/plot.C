#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>
#include "TGraphErrors.h"
#include "TH1.h"

const Int_t muon_no = 100;

double biasv(double *xv, double *p){
  double x=*xv; 
  double tau1=p[1];
  //  double tau2=p[2];
  double G0=p[2];
  double x0=p[3];
  //  cout<<"BIAS "<<x0<<" \n";
  double gain=G0;
  
  if(x>p[3]){
    gain=G0*(1.-p[0]);
    gain=1.-(1-gain)*exp(-(x-p[3])/tau1);
    
  }
    return gain;		 
}

 void plot(){
   //   TF1 *xplot = new TF1("xplot",biasv,-5,10,4);
     TF1 *xplot = new TF1("xplot",biasv,-2,1.,4);
     xplot->SetParameters(0.001,3.,1.,0.);
     xplot->Draw();
     xplot->Print();
 }

//void extraction_mu(){
vector<double> extraction_mu(){
  TH1D *hpar2=new TH1D("hpar2","exp",muon_no,0.,700.);
   
  vector<double> xt(muon_no,0);
  //double Ntot = muon_no;
  double tau=64.4;
  for(int i=0;i<muon_no;i++){
    //    x0=gRandom->Gaus()*0.01;
    xt[i]=gRandom->Exp(tau);
    double x0=xt[i];
    //cout<<i<<":"<<xt[i]<<" \n";
    //hpar2->Fill(x0);
    //hpar2->Draw();
  }
  sort(xt.begin(), xt.end());
  //  cout<<"SORT \n";
  //for (int i = 0; i != xt.size(); ++i)
  //cout << xt[i] << " ";
  return xt;
}



void pulse_laser_old(int fill=16, int laser_no=8, int cycles=2000){
    
  //printf("Started \n");
  //static int repeat=200;//muon_no*fill*cycles;
  
  //Define 2d vectors for storing cycle numbers and fill # - the suffix 2dv => a 2-d vector for each cycle and fill number
  vector <vector <vector<double> > >  xx(cycles, vector <vector<double> >(fill,extraction_mu()));
  //Necessary to initialize this way or the code breaks at runtime
  vector <vector <int> > laser(fill, vector <int>(9,0));
  //vector<double> xx=extraction_mu();
  //List of final 1-d vectors for each fill of muon_no muons
  vector <double>  x(muon_no,0);
  vector <double>  energy(muon_no,0.001);
  //double energy;
    TGraph *gr;
    vector <double> G0(muon_no,1.); vector <double> x0(muon_no,0.);
    vector <double> ex(muon_no,0); vector <double> ey(muon_no,0.);
    vector <double> G(muon_no,1.);
  
    TCanvas* c2 = new TCanvas("c2");
    c2->DrawFrame(-2., 0.95, 150., 1, "Global Title;X Axis Title;Y Axis Title");

    
  for(int j=0;j<fill;j++)
    for (int l=0; l<=laser_no; l++){
      laser[j][l]=80*l+5*j;
      //cout<<laser[j][l]<<"\n";
    }
  //cout<<laser[3][2]<<"\n";

  for(int k=0;k<cycles;k++){

    for(int j=0;j<fill;j++){

      xx[k][j]=extraction_mu();
      
      for(int i=0;i<muon_no;i++){
          x[i]=xx[k][j][i];
          energy[i]=0.001;
	  ex[i]=0;
          //cout<<energy[i];<<"\n";
          //energy=0.001;
          for (int l=0; l<=laser_no; l++){
              x[laser[j][l]]=laser[j][l];
              //energy=0.005;
	  //energy[laser[j][l]]=0.005;
	}
	
	//cout<<"Samples:"<<i<<"="<<x[i]<<"\n";
      }//Muon loop ends here as to initialize..
      //Next we sort in time with energy corresponding to 0.005 for lasers.
      //cout<<x.size()<<"\n";
      //Make a pair of time and energy first and then sort them
      sort(x.begin(),x.end());
      
        
      //2nd loop over muons to assign energies with lasers
      for(int i=0;i<muon_no;i++){
          energy[i]=0.001;
          for (int l=0; l<=laser_no; l++){
              if(x[i]==laser[j][l]) energy[i]=0.005;
            
          }
          cout<<"Samples:"<<i<<"="<<x[i]<<" E="<<energy[i]<<"\n";
      }// end of 2nd muon loop
        
      for(int i=0;i<muon_no;i++){
          if(i==0) {
              TF1 *xplot = new TF1("xplot",biasv,0.,x[1],4);
	      //TH1D *xhist = new TH1D("xhist","",700,0,700);
              xplot->SetParameters(energy[i],3.,1.,x[0]);
              // xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
              xplot->Draw("same");
	      x0[i+1]=x[0];
	      //     G[i]=1;
	      G0[i+1]=1.;
	      ey[i+1]=sqrt(G0[i+1]);
	  }
	  else {
	    if(i<muon_no-1){
	      double xmin=x0[i]*0.99;
	      TF1 *xplot = new TF1("xplot",biasv,xmin,x[i],4);
	      xplot->SetParameters(energy[i],3.,G0[i],x0[i]);
        
	      G[i]=xplot->Eval(x[i]);
	      x0[i+1]=x[i];
	      G0[i+1]=G[i];
	      ey[i+1]=1./sqrt(G0[i+1]);
          cout<<i<<" G0= "<<G0[i]<<" x0= "<<x0[i]<<" G= "<<G[i]<<" x= "<<x0[i+1]<<" \n";
	      //if(j%2==0)xplot->SetLineColor(kBlue);
          //if(j%3==0)xplot->SetLineColor(kGreen);
          //if(j%5==0)xplot->SetLineColor(kBlack);
	      xplot->Draw("R+,same");
	    } 
	  }
	  
	}
    //int n=muon_no;
    //gr = new TGraph(n,&x0[0],&G0[0]);// addresses cos we use vectors
    //gr->GetXaxis()->SetRangeUser(-2.,200.);
    //gr->GetYaxis()->SetRangeUser(0.95,1.);
    //gr->Draw("AC");
	cout<<"***************End of fill "<<j+1<<" ****************"<<"\n";
        //c2->Update();
      }//End of fill
    cout<<"***************End of cycle "<<k+1<<" ****************"<<"\n";
  }//end of cycle
  int n=muon_no;
  gr = new TGraph(n,&x0[0],&G0[0]);// addresses cos we use vectors
  gr->GetXaxis()->SetRangeUser(-2.,200.);
  gr->GetYaxis()->SetRangeUser(0.95,1.);
  gr->Draw("*, same");
  
  c2->Update();
  TFile *f = TFile::Open("pulse_laser_1000.root", "RECREATE");
  gr->Write("Laser");
  }
  
  
 
void pulse_laser_only(int fill=16, int laser_no=8, int cycles=20){
  //This function is only for the laser pulses fired every 80 us - as the name suggests
  TGraph *gr;
  TProfile *pf;
  TH1D *hf;;
  vector <double> fill_time(muon_no, 0);
  vector <double> fill_time_avg(muon_no, 0);
  vector <double> cycle_time(muon_no, 0);
  vector <double> cycle_time_avg(muon_no, 0);
  
  
  vector <Double_t> fill_gain(muon_no, 0);
  vector <Double_t> fill_gain_avg(muon_no, 0);
  vector <Double_t> cycle_gain(muon_no, 0);
  vector <Double_t> cycle_gain_avg(muon_no, 0);
  
  
  vector <Double_t> fill_gain_sigma(muon_no, 0);
  vector <Double_t> cycle_gain_sigma(muon_no, 0);
  vector <Double_t> time(muon_no, 0);
  vector <Double_t> gain(muon_no, 0);  
  
  vector <vector <vector<double> > >  xx(cycles, vector <vector<double> >(fill,extraction_mu()));
  //Necessary to initialize this way or the code breaks at runtime
  vector <vector <int> > laser(fill, vector <int>(9,0));
  //vector<double> xx=extraction_mu();
  //List of final 1-d vectors for each fill of muon_no muons
  vector <double>  x_all(muon_no,0);//All muons
  vector <double>  x(muon_no,0);//Muons only for laser shots at 5, 10, 15 etc. us
  vector <double>  energy(muon_no,0.001);
  int muon_index=0;
  vector <double> G0(muon_no,1.); vector <double> x0(muon_no,0.);
  vector <double> ex(muon_no,0); vector <double> ey(muon_no,0.);//for errors
  vector <double> G(muon_no,1.);
  
  TCanvas* c2 = new TCanvas("c2");
  c2->DrawFrame(-2., 0.95, 150., 1, "Global Title;X Axis Title;Y Axis Title");
  
  TH1D *xhist1 = new TH1D("xhist1","xxx",100,0,700);
  TH1D *xhist2 = new TH1D("xhist2","xxx",100,0,700);

  for(int k=0;k<cycles;k++){
    for(int i=0;i<muon_no;i++){
      fill_time[i]=0;
      fill_gain[i]=0;
    }
    for(int j=0;j<fill;j++){
      for (int l=0; l<=laser_no; l++){
	laser[j][l]=80*l+5*j;
	xx[k][j]=extraction_mu();
      
	for(int i=0;i<muon_no;i++){
	  //xx[k][i]=extraction_mu();
	  //if(i==80*l+5*j){
	    x_all[i]=xx[k][j][i];
	    energy[i]=0.005;
	    ex[i]=0;
	    if(int(x_all[i])==80*l+5*j)
	       x[i]=xx[k][j][i];

	    //Next we sort in time with energy corresponding to 0.005 for lasers.
	    //cout<<x.size()<<"\n";
	    //Make a pair of time and energy first and then sort them
	    sort(x.begin(),x.end());
	    
	    //}
	    //cout<<energy[i];<<"\n";
          //energy=0.001;
	/* for (int l=0; l<=laser_no; l++){
              x[laser[j][l]]=laser[j][l];
              //energy=0.005;
	      energy[laser[j][l]]=0.005;
	  }
	*/
	cout<<"Samples:"<<i<<"="<<x[i]<<"\n";
      }//Muon loop ends here as to initialize..
      
     
      for(int i=0;i<muon_no;i++){
          if(i==0) {
              TF1 *xplot = new TF1("xplot",biasv,0.,x[1],4);
	      //TH1D *xhist = new TH1D("xhist","",700,0,700);
              xplot->SetParameters(energy[i],3.,1.,x[0]);
              // xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
              //xplot->Draw("same");
	      x0[i+1]=x[0];
	      //     G[i]=1;
	      G0[i+1]=1.;
	      ey[i+1]=sqrt(G0[i+1]);
	  }
	  else {
	    if(i<muon_no-1){
	      double xmin=x0[i]*0.99;
	      TF1 *xplot = new TF1("xplot",biasv,xmin,x[i],4);
	      xplot->SetParameters(energy[i],3.,G0[i],x0[i]);
        
	      G[i]=xplot->Eval(x[i]);
	      x0[i+1]=x[i];
	      G0[i+1]=G[i];
	      ey[i+1]=1./sqrt(G0[i+1]);
	      //cout<<i<<" G0= "<<G0[i]<<" x0= "<<x0[i]<<" G= "<<G[i]<<" x= "<<x0[i+1]<<" \n";
	      //if(j%2==0)xplot->SetLineColor(kBlue);
          //if(j%3==0)xplot->SetLineColor(kGreen);
          //if(j%5==0)xplot->SetLineColor(kBlack);
	  //xplot->Draw("R+,same");
	    } 
	  }
	 fill_time[i]+=x0[i];
	 fill_gain[i]+=G0[i];
	 //fill_time.push_back(x0[i]);
	 
	 fill_time_avg[i]=fill_time[i]/(j+1);
	 fill_gain_avg[i]=fill_gain[i]/(j+1);  
      }
    //int n=muon_no;
    //gr = new TGraph(n,&x0[0],&G0[0]);// addresses cos we use vectors
    //gr->GetXaxis()->SetRangeUser(-2.,200.);
    //gr->GetYaxis()->SetRangeUser(0.95,1.);
    //gr->Draw("AC");
      //cout<<"***************End of fill "<<j+1<<" ****************"<<"\n";
      }
      // c2->Update();
    }//End of fill
    
    cout<<"***************End of cycle "<<k+1<<" ****************"<<"\n";
    for(int i=0;i<muon_no;i++){
	  cycle_time[i]+=fill_time_avg[i];
	  cycle_time_avg[i]=cycle_time[i]/(k+1);
	  
	  cycle_gain[i]+=fill_gain_avg[i];
	  cycle_gain_avg[i]=cycle_gain[i]/(k+1);

	  fill_gain_sigma[i]=sqrt((fill_gain_avg[i]-cycle_gain_avg[i])*(fill_gain_avg[i]-cycle_gain_avg[i])/(k+1));	     
	  //cout<<"Muon no:"<<i<<" Fill time " << fill_time[i] <<" Size of fill time "<<fill_time.size()<<" Avg Fill time " << fill_time_avg[i]<<" Avg cycle time " << cycle_time_avg[i]<<"//////////"<<"\n";
	}
	

    
  }//end of cycle
  /*
  int n=muon_no;
  gr = new TGraph(n,&x0[0],&G0[0]);// addresses cos we use vectors
  gr->GetXaxis()->SetRangeUser(-2.,200.);
  gr->GetYaxis()->SetRangeUser(0.95,1.);
  gr->Draw("*, same");
  */
 hf = new TH1D("hf","Hist of time versus gain",700,0,700);
    for ( Int_t i=0; i<cycle_time_avg.size(); i++){
      //cout<<"Muon no:"<<i<<" Err Before: "<<fill_gain_sigma[i]<<"//////////"<<"\n";
      
      hf->SetBinContent(hf->GetBin(cycle_time_avg[i]), cycle_gain_avg[i]);
      //pf->SetBinError( pf->FindBin(cycle_time_avg[i]),fill_gain_sigma[i]);
      //pf->Fill(cycle_time_avg[i], cycle_gain_avg[i],fill_gain_sigma[i]);
      //cout<<"Muon no:"<<i<<" Err Now: "<<fill_gain_sigma[i]<<"//////////"<<"\n";
      hf->SetBinError(hf->GetBin(cycle_time_avg[i]),fill_gain_sigma[i]);
	
      //cout<<"Muon no:"<<i<<" Avg cycle time " << cycle_time_avg[i]<<" Avg cycle gain " << cycle_gain_avg[i]<<" Err: "<<hf->GetBinError(i)<<"//////////"<<"\n";
      
      //pf->Fill(fill_time_avg[i], fill_gain_avg[i],1);
    }
    hf->Draw("e");
  
  c2->Update();
  //TFile *f = TFile::Open("pulse_laser_1000.root", "RECREATE");
  //gr->Write("Laser");
  }
  
  


void sample(int fill=16, int cycles=10){
    
    TGraph *gr;
    TProfile *pf;
    TH1D *hf;;
    vector <double> fill_time(muon_no, 0);
    vector <double> fill_time_avg(muon_no, 0);
    vector <double> cycle_time(muon_no, 0);
    vector <double> cycle_time_avg(muon_no, 0);
    
    
    vector <Double_t> fill_gain(muon_no, 0);
    vector <Double_t> fill_gain_avg(muon_no, 0);
    vector <Double_t> cycle_gain(muon_no, 0);
    vector <Double_t> cycle_gain_avg(muon_no, 0);
   
    
    vector <Double_t> fill_gain_sigma(muon_no, 0);
    vector <Double_t> cycle_gain_sigma(muon_no, 0);
    vector <Double_t> time(muon_no, 0);
    vector <Double_t> gain(muon_no, 0);
    
    vector <vector <vector<double> > >  xx(cycles, vector <vector<double> >(fill,extraction_mu()));
    //Necessary to initialize this way or the code breaks at runtime
    vector <vector <int> > laser(fill, vector <int>(9,0));
    //vector<double> xx=extraction_mu();
    //List of final 1-d vectors for each fill of muon_no muons
    vector <double>  x(muon_no,0);
    vector <double>  energy(muon_no,0.001);
    //double energy;

    int muon_index=0;

    vector <double> G0(muon_no,1.); vector <double> x0(muon_no,0.);
    vector <double> ex(muon_no,0); vector <double> ey(muon_no,0.);
    vector <double> G(muon_no,1.);
    
    TCanvas* c2 = new TCanvas("c2");
    c2->DrawFrame(-2., 0.98, 300., 1, "Global Title;X Axis Title;Y Axis Title");
    
    TH1D *xhist1 = new TH1D("xhist1","xxx",100,0,700);
    TH1D *xhist2 = new TH1D("xhist2","xxx",100,0,700);
    //TH1D *xhist = new TH1D("xhist","xxx",100,0,700);
    TF1 *AFT = new TF1("AFT","1-0.045*(TMath::Exp(-x/64.4) - TMath::Exp(-x/3.))",0,700);
    for(int k=0;k<cycles;k++){
     	for(int i=0;i<muon_no;i++){
	      fill_time[i]=0;
	      fill_gain[i]=0;
	    }
        for(int j=0;j<fill;j++){         
            xx[k][j]=extraction_mu();

            for(int i=0;i<muon_no;i++){
	      x[i]=xx[k][j][i];
	      energy[i]=0.001;
	      ex[i]=0;
	      cout<<"Samples:"<<i<<"="<<x[i]<<" E="<<energy[i]<<"\n";
	      if(i==0) {
		TF1 *xplot = new TF1("xplot",biasv,0.,x[1],4);
		//TH1D *xhist = new TH1D("xhist","",100,0,x[1]);
		xplot->SetParameters(energy[i],3.,1.,x[0]);
		//xhist->Eval(xplot);
		//       xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
		//xplot->Draw("same");
		x0[i+1]=x[0];
		//     G[i]=1;
		G0[i+1]=1.;
		//ey[i+1]=sqrt(G0[i+1]);
		//fill_time.push_back(x0[i]);
	      }
	      else {
		if(i<muon_no-1){
		  double xmin=x0[i]*0.99;
		  TF1 *xplot = new TF1("xplot",biasv,xmin,x[i],4);
		  //TH1D *xhist = new TH1D("xhist","",100,xmin,x[i]);
		  xplot->SetParameters(energy[i],3.,G0[i],x0[i]);
                  //xhist->Eval(xplot);
		  G[i]=xplot->Eval(x[i]);
		  x0[i+1]=x[i];
		  G0[i+1]=G[i];
		  ey[i+1]=(1-G0[i+1])/sqrt(k+1);
		  //cout<<i<<" G0= "<<G0[i]<<" x0= "<<x0[i]<<" G= "<<G[i]<<" x= "<<x0[i+1]<<" \n";
		  //if(j%2==0)xplot->SetLineColor(kBlue);
		  //if(j%3==0)xplot->SetLineColor(kGreen);
		  //if(j%5==0)xplot->SetLineColor(kBlack);
		  xplot->SetLineColor(kBlue);
		  //xplot->Draw("R+,same");
		  //xhist
		  //xhist->Draw("same");
		} 
	      }
	      fill_time[i]+=x0[i];
	      fill_gain[i]+=G0[i];
	      //fill_time.push_back(x0[i]);
	      
	      fill_time_avg[i]=fill_time[i]/(j+1);
	      fill_gain_avg[i]=fill_gain[i]/(j+1);

	      int n=muon_no;
	      gr = new TGraph(n,&x[0],&G[0]);
	      gr->GetXaxis()->SetRangeUser(-2.,200.);
	      gr->GetYaxis()->SetRangeUser(0.95,1.);
	      gr->SetMarkerStyle(10);
	      // gr->Draw("AP");

	      //cout<<"Muon no:"<<i<<" Fill time " << fill_time[i] <<" Size of fill time "<<fill_time.size()<<" Avg Fill time " << fill_time_avg[i]<<" Avg cycle time " << cycle_time_avg[i]<<"//////////"<<"\n";
	      }//Muon loop for function plotting
            //cout<<"***************End of fill "<<j+1<<"****************"<<"\n";
	   
        }//End of fill	      
	for(int i=0;i<muon_no;i++){
	  cycle_time[i]+=fill_time_avg[i];
	  cycle_time_avg[i]=cycle_time[i]/(k+1);
	  
	  cycle_gain[i]+=fill_gain_avg[i];
	  cycle_gain_avg[i]=cycle_gain[i]/(k+1);

	  fill_gain_sigma[i]=sqrt((fill_gain_avg[i]-cycle_gain_avg[i])*(fill_gain_avg[i]-cycle_gain_avg[i])/(k+1));	     
	  //cout<<"Muon no:"<<i<<" Fill time " << fill_time[i] <<" Size of fill time "<<fill_time.size()<<" Avg Fill time " << fill_time_avg[i]<<" Avg cycle time " << cycle_time_avg[i]<<"//////////"<<"\n";
	}
	
        //cout<<"***************End of cycle "<<k+1<<" ****************"<<"\n";
	
    }//end of cycle
    //Only in case you want a graph.....
    /*
    int n=muon_no;
    gr = new TGraphErrors(n,&x[0],&G[0],&ex[0],&ey[0]);
    gr->GetXaxis()->SetRangeUser(-2.,200.);
    gr->GetYaxis()->SetRangeUser(0.95,1.);
    gr->SetMarkerStyle(10);
    gr->Draw("*");
    gr->SetLineColor(kBlue);
    //gr->Draw("AC*");
    */
    //AFT->Draw("same");
    //pf = new TProfile("pf","Profile of time versus gain",700,0,700,"S");
    hf = new TH1D("hf","Hist of time versus gain",700,0,700);
    for ( Int_t i=0; i<cycle_time_avg.size(); i++){
      //cout<<"Muon no:"<<i<<" Err Before: "<<fill_gain_sigma[i]<<"//////////"<<"\n";
      
      hf->SetBinContent(hf->GetBin(cycle_time_avg[i]), cycle_gain_avg[i]);
      //pf->SetBinError( pf->FindBin(cycle_time_avg[i]),fill_gain_sigma[i]);
      //pf->Fill(cycle_time_avg[i], cycle_gain_avg[i],fill_gain_sigma[i]);
      //cout<<"Muon no:"<<i<<" Err Now: "<<fill_gain_sigma[i]<<"//////////"<<"\n";
      hf->SetBinError(hf->GetBin(cycle_time_avg[i]),fill_gain_sigma[i]);
	
      //cout<<"Muon no:"<<i<<" Avg cycle time " << cycle_time_avg[i]<<" Avg cycle gain " << cycle_gain_avg[i]<<" Err: "<<hf->GetBinError(i)<<"//////////"<<"\n";
      
      //pf->Fill(fill_time_avg[i], fill_gain_avg[i],1);
    }
    hf->Draw("e");
    c2->Update();
    TFile *f = TFile::Open("sample_only_100.root", "RECREATE");
    //gr->Write("MySample");
    
}



void sample_old(){
    
    double par[4]={0.001,3.0,1.,0.};
    vector<double> xx=extraction_mu();
    
    double x[muon_no];
    for(int i=0;i<muon_no;i++){
        x[i]=xx[i];
        //x[i]=5*i;
    }
    for (int i = 0; i<muon_no; ++i) cout<<"Sample: "<< x[i] << "\n ";
    
    //  vector<double> G0(10,1),x0(10,0.),ex(10,0.),ey(10,0.);
    double G0[muon_no]={1.},x0[muon_no]={0.},ex[muon_no]={0},ey[muon_no]={0.};
    double G[muon_no]={1.};
    //  cout<<"SAMPLE "<<x[]<<" \n";
    //  double x[]=extraction_mu();
    TCanvas* c2 = new TCanvas("c2");
    c2->DrawFrame(-2., 0.65, 150., 1.0002, "Global Title;X Axis Title;Y Axis Title");
    for(int i=0;i<muon_no;i++){
        if(i==0) {
            TF1 *xplot = new TF1("xplot",biasv,0.,x[1],4);
            xplot->SetParameters(0.001,3.,1.,x[0]);
            //       xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
            xplot->Draw("same");
            x0[i+1]=x[0];
            //     G[i]=1;
            G0[i+1]=1.;}
        else {
            if(i<muon_no-1){
                double xmin=x0[i]*0.99;
                TF1 *xplot = new TF1("xplot",biasv,xmin,x[i],4);
                xplot->SetParameters(0.001,3.,G0[i],x0[i]);
                
                G[i]=xplot->Eval(x[i]);
                x0[i+1]=x[i];
                G0[i+1]=G[i];
                cout<<i<<" G0= "<<G0[i]<<" x0= "<<x0[i]<<" G= "<<G[i]<<" x= "<<x0[i+1]<<" \n";
                // TF1 *xplot = new TF1("xplot",biasv,x[i],x[i+1],4);
                // xplot->SetParameters(0.001,3.,G[i],x[i]);
                xplot->SetLineColor(kBlue);
                xplot->Draw("R+,same");
            } 
        }
        
        
    }
    
    int n=muon_no;
    TGraph* gr = new TGraph(n, x0, G0); // because we use vectors here
    gr->GetXaxis()->SetRangeUser(-2.,200.);
    gr->GetYaxis()->SetRangeUser(0.95,1.);
    gr->Draw("AC");
}


