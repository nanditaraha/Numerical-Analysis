/*
Simple code that finds the time (in mins) required for laser calibration with a certain frequency of laser pulses as input
5 micro sec sampling pt.
 */

#include <TRandom.h>
#include <TMath.h>
#include "TApplication.h"
#include <TRandom3.h>
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
//#include<TAttParticle.h>

//Input laser frequency in kHz and number of events for required precision in uncertainty
void findTime(double freq=10, int events=2000) {
  //Variables
  double laser_pulse_time; //Convert laser frequncy from kHz to us
  int no_pts;// # of points in one fill
  int no_fills=0;
  double time;// in sec
  //laser_pulse_time = 1.*1e3/freq;
  no_pts = TMath::Nint(700*freq/1000.) ; // Used the usual round-off
  int fills = TMath::Nint(140./no_pts) ;
  //no_pts = 700./laser_pulse_time ;
  //int fills = 140./no_pts ;
  double fill_gap_t = 11e-3; //ms to s
  double bunch_gap_t = 66.7e-3;//ms to s - bunch means 4 events + gap time
  double cycle_gap_t = 1.06; //s
 
for(int no_fills=0;no_fills<fills;no_fills++){
    double tot_time;
    int fill_gap_no = 0;
    int bunch_gap_no=0;
    int cycle_gap_no = 1;

    fill_gap_no = no_fills%4;
    //if(no_fills>15)
    cycle_gap_no=(no_fills+1)/16 + 1;
    bunch_gap_no=(no_fills-16*(cycle_gap_no-1))/4;
    if((no_fills+1)%16==0) time = (cycle_gap_no-1)*cycle_gap_t;
    else
        time = ((fill_gap_no+1)*fill_gap_t + bunch_gap_no*bunch_gap_t + (cycle_gap_no-1)*cycle_gap_t);
        //printf("%d: Fill# %d, Bunch # %d, cycle# %d\n",no_fills+1,fill_gap_no+1,bunch_gap_no,cycle_gap_no-1);
        //printf("\t Fill_t %f, Bunch_t %f, cycle_t %f time %f\n",(fill_gap_no+1)*fill_gap_t, bunch_gap_no*bunch_gap_t, (cycle_gap_no-1)*cycle_gap_t,time);
  }
 	time = time*events/60.;//converting time to mins depending on no. of events
	printf("Pulse time %f and no of pts %d and no of fills %d and time %f\n",laser_pulse_time,no_pts,fills, time); 
}
