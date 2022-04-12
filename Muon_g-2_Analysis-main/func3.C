#include <iostream>
#include <algorithm>
#include <vector>
#include "TGraphErrors.h"
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
    TH1D *hpar2=new TH1D("hpar2","exp",100,0.,700.);
    vector<double> xt(100,0);
    double Ntot = 100;
    double tau=64.4;
    for(int i=0;i<100;i++){
        //    x0=gRandom->Gaus()*0.01;
        xt[i]=gRandom->Exp(tau);
        double x0=xt[i];
        cout<<xt[i]<<" \n";
        hpar2->Fill(x0);
        hpar2->Draw();
    }
    sort(xt.begin(), xt.end());
    //  cout<<"SORT \n";
    //for (int i = 0; i != xt.size(); ++i)
        //cout << xt[i] << " ";
    return xt;
}

//This needs to be completely written........

void pulse_laser(int fill=16, int laser_no=8, int cycles=2000){
    
    vector<double> xx=extraction_mu();
    vector <double> x(100,0);
    vector <int> laser;
    vector <double> energy(100,0);
    
    vector <double> x_sort;
    vector <double> energy_sort;
    //double it[100];
    pair <double,double> apair;
    vector< pair<double,double>  > time_energy; //This forms a pair of energy and time so that sorting with time does not change the corresponding energy
    
    for(int i=0;i<cycles;i++)
        for(int j=0;j<100;j++)
            x[j]=xx[j];
            
    for(int i=0;i<100;i++) {
        //x[i]=xx[i];
        energy[i]=0.001;
    }
    
     for (int j=0; j<fill; j++) {
        for (int k=0; k<=laser_no; k++)
            if(80*k+5*j<100)
            laser.push_back(80*k+5*j);//=80*k+5*j;
     }
    
    for(int i=0;i<laser.size();i++)
    {
        x[laser[i]]=laser[i];
        energy[laser[i]]=0.005;
        //cout<<"IN Loop Sample:"<<i << "="<<x[laser[i]]<<" Energy "<<energy[laser[i]]<<" laser index"<<laser[i]<<"\n ";
    }
    
    for(int j=0;j<x.size();j++)
    {
        apair=make_pair(x[j],energy[j]);
        time_energy.push_back(apair);
    }
    
    sort(time_energy.begin(),time_energy.end());
    
    
    for (vector<pair<double,double> >::iterator it2 = time_energy.begin(); it2 != time_energy.end(); ++it2) {
        apair = *it2;
        //cout << "(" << apair.first << "," << apair.second << ") ; "<<"\n";
        x_sort.push_back(apair.first);
        energy_sort.push_back(apair.second);
        //i++;
    }
    
    
    double G0[100]={1.},x0[100]={0.},ex[100]={0},ey[100]={0.};
    double G[100]={1.};
    
    for (int i = 0; i<100; ++i) cout<<"Sample:"<<i << "="<<x_sort[i]<<" Energy "<<energy_sort[i]<<"\n ";
    
    //cout<<"SAMPLE "<<x[]<<" \n";
    //  double x[]=extraction_mu();
    TCanvas* c2 = new TCanvas("c2");
    c2->DrawFrame(-2., 0.98, 150., 1, "Global Title;X Axis Title;Y Axis Title");
    
    for(int i=0;i<100;i++){
        if(i==0) {
            TF1 *xplot = new TF1("xplot",biasv,0.,x_sort[1],4);
            xplot->SetParameters(energy_sort[i],3.,1.,x_sort[0]);
            //       xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
            xplot->Draw("same");
            x0[i+1]=x_sort[0];
            //     G[i]=1;
            G0[i+1]=1.;}
        else {
            if(i<99){
                double xmin=x0[i]*0.99;
                TF1 *xplot = new TF1("xplot",biasv,xmin,x_sort[i],4);
                xplot->SetParameters(energy_sort[i],3.,G0[i],x0[i]);
                
                G[i]=xplot->Eval(x_sort[i]);
                x0[i+1]=x_sort[i];
                G0[i+1]=G[i];
                cout<<i<<" G0= "<<G0[i]<<" x0= "<<x0[i]<<" G= "<<G[i]<<" x= "<<x0[i+1]<<" \n";
                // TF1 *xplot = new TF1("xplot",biasv,x[i],x[i+1],4);
                // xplot->SetParameters(0.001,3.,G[i],x[i]);
                if(i==laser[i])
                    xplot->SetLineColor(kBlue);
                xplot->Draw("R+,same");
            } 
        }
        
        
    }


}

void laser_only(int fill=16, int laser_no=8){
    
    double par[4]={0.001,3.0,1.,0.};

    double x[100];
    
    for (int j=0; j<fill; j++) {
        for (int k=0; k<=laser_no; k++)
            if(80*k+5*j<100)
                x[80*k+5*j]=80*k+5*j;
    }

    double G0[100]={1.},x0[100]={0.},ex[100]={0},ey[100]={0.};
    double G[100]={1.};
    //  cout<<"SAMPLE "<<x[]<<" \n";
    //  double x[]=extraction_mu();
    TCanvas* c2 = new TCanvas("c2");
    c2->DrawFrame(-2., 0.65, 150., 1.0002, "Global Title;X Axis Title;Y Axis Title");
    for(int i=0;i<100;i++){
        if(i==0) {
            TF1 *xplot = new TF1("xplot",biasv,0.,x[1],4);
            xplot->SetParameters(0.001,3.,1.,x[0]);
            //       xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
            xplot->Draw("same");
            x0[i+1]=x[0];
            //     G[i]=1;
            G0[i+1]=1.;}
        else {
            if(i<99){
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
}

void sample(){

  double par[4]={0.001,3.0,1.,0.};

  vector<double> xx=extraction_mu();

  double x[100];
   for(int i=0;i<100;i++){
     x[i]=xx[i];
     //x[i]=5*i;
       }
    for (int i = 0; i<100; ++i) cout<<"Sample: "<< x[i] << "\n ";
    
    //  vector<double> G0(10,1),x0(10,0.),ex(10,0.),ey(10,0.);
  double G0[100]={1.},x0[100]={0.},ex[100]={0},ey[100]={0.};
  double G[100]={1.};
  //  cout<<"SAMPLE "<<x[]<<" \n";
  //  double x[]=extraction_mu();
  TCanvas* c2 = new TCanvas("c2");
   c2->DrawFrame(-2., 0.65, 150., 1.0002, "Global Title;X Axis Title;Y Axis Title");
  for(int i=0;i<100;i++){
    if(i==0) {
       TF1 *xplot = new TF1("xplot",biasv,0.,x[1],4);
       xplot->SetParameters(0.001,3.,1.,x[0]); 
       //       xplot->GetHistogram()->GetXaxis()->SetRangeUser(-2.,10.);
       xplot->Draw("same");
     x0[i+1]=x[0];
     //     G[i]=1;
    G0[i+1]=1.;}
     else {
       if(i<99){
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
    
    int n=100;
    TGraph* gr = new TGraph(n,x,G0);
    gr->GetXaxis()->SetRangeUser(-2.,200.);
    gr->GetYaxis()->SetRangeUser(0.993,1.);
}


