#include <iostream>
#include <algorithm>
#include <vector>
#include "TGraphErrors.h"

//This the laser distribution function...
double laser_func(double *xv, double *p){
    double x=*xv;
    double y;
    if(x<p[0])
        y = 0;
    //if(x==p[0])
    else if (x>(p[0]+p[1]))
        y = 0;
    else
       y= p[2];
    return y;
}



//This needs to be completely written........
void pulse_laser(){
    //TF1 *xplot = new TF1("xplot",laser_func,xmin,x[i],1);
    //xplot->SetParameter(0,0.005);
    TCanvas* c2 = new TCanvas("c2");
    c2->DrawFrame(-0, -0, 700000, 0.006, "Global Title;X Axis Title;Y Axis Title");
    for(int i=0;i<140;i++){//time bins (gaps in time of 5 us)
        //for (Int_t j = 0; j < 2000; j++) {
            //TF1 *xplot = new TF1("xplot",laser_func,0,10000,2);
            //xplot->SetParameters(5000,500);
    //if(i==0){
            TF1 *f1 = new TF1("f1",laser_func,5000*i,5000*(i+2),3);
            f1->SetNpx(2000);
            f1->SetParameters(5000*(i+1)-5,5,0.005);
            f1->Draw("same");
    //}
    /*else {
            TF1 *f1 = new TF1("f1",laser_func,9990,20000,3);
            f1->SetNpx(1000);
            f1->SetParameters(9995,5,0.005);
            //f1->SetLineColor(kBlue);
            f1->Draw("R+,same");

            //xplot->Draw();
            //xplot->Draw("R+,same");
            //double data = xplot->Eval(i*5000);
            //h_pmt->SetBinContent(j,pmt);
            //sigma = 0.02*gain_ratio;//assume 2% spread in sigma..
            //h_test->Fill(data);
    }*/
    //}

    }
}
