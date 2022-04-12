
#include <TLatex.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <fstream>

#include<math.h>

void larmor();
void thomas();
TCanvas *c1 = new TCanvas("c1","c1",0,0,50,100);
TLatex Tl; //Tl.SetTextFont(43); Tl.SetTextSize(20);
Double_t dy = 0.3;

void talk()
{
    Tl.SetTextFont(43); Tl.SetTextSize(20);
    //larmor();
    thomas();
    
}

void larmor(){
    Tl.DrawLatex(0.1,3*dy,"#mu = - #frac{g e}{2 m} L");
    Tl.DrawLatex(0.1,2*dy,"L sin#theta #frac{#Delta#phi}{#Delta t} = - #frac{g e}{2 m} L B sin#theta");
    Tl.DrawLatex(0.1,dy+0.1,"#omega_{L} = - #frac{g e}{2 m} B");
}

void thomas(){
    Tl.DrawLatex(0.1,dy,"#Rightarrow #omega_{T} = #frac{v^{2}}{c^{2}} #left(#frac{#gamma^{2}}{#gamma + 1}#right) #frac{e}{#gamma m} B = (#gamma -1)#frac{e}{#gamma m} B");
    Tl.DrawLatex(0.1,2*dy,"|#vec{a} x #vec{v}| = #frac{e}{#gamma m} v^{2} B");
    Tl.DrawLatex(0.1,3*dy,"a = #frac{e}{#gamma m} v B");
}