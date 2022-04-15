## Numerical analysis and Monte Carlo simulations of the energy gain fluctuation in SiPMs
This study derives how the gain fluctuates in the SiPM used for reading the electron calorimeters in the Muon g-2 experiment. The numerical solutions 
for this gain is then simulated and finally compared with the experimental data.

### Gain variation in a fill: analytical approach
The distribution of the electron/positron as a function of electron energy E (i.e. the wiggle plot), decay asymmetry A and a constant phase term φ that depends on the drift time of the electrons is given by,</br>
<img width="290" alt="Screen Shot 2022-04-11 at 2 35 09 PM" src="https://user-images.githubusercontent.com/27436642/162806505-d9efa68e-3a65-4b84-bddd-010c89462e8b.png"></br>
where N, A and φ depend on the energy. Expanding this</br>
<img width="327" alt="Screen Shot 2022-04-11 at 2 37 31 PM" src="https://user-images.githubusercontent.com/27436642/162806914-067718df-82a5-4147-a2aa-a787a313af15.png"></br>
where n0 is the number of pluses in a fill and θ(t − t<sub>i</sub>) is the Heaviside step function to include pulses sorted in increasing order of time.
We obtain the average gain function by averaging over all times t<sub>i</sub> (which are assumed to be the same) and the entire energy range. This is given by,</br>
<img width="468" alt="Screen Shot 2022-04-11 at 2 39 15 PM" src="https://user-images.githubusercontent.com/27436642/162807192-d163e4aa-7baf-4cec-adba-5341cc359095.png"></br>
Gain functions for τ<sub>k</sub> = 1 μs, τ<sub>k</sub> = 3 μs, τ<sub>k</sub> = 10 μs and τ<sub>k</sub> = 20 μs corresponding to the green, blue, red and black colors respectively. </br>
<p align = "center">
<img width="406" alt="Screen Shot 2022-04-11 at 1 54 55 PM" src="https://user-images.githubusercontent.com/27436642/162800199-af362a3b-5794-42d7-bf4a-0e6cc1cfcf80.png"></br>
 Fig. 1
</p>

### Simulation of gain function
We produced a ***Monte Carlo simulation*** to study the bias voltage (BV) sagging of the SiPMs (Silicon Photo Multipliers) . It is the convolution of a single energy drop with the time distribution of the positrons. The cumulative gain (for n<sub>o</sub> pulses) is,</br>
<img width="341" alt="Screen Shot 2022-04-11 at 7 51 16 PM" src="https://user-images.githubusercontent.com/27436642/162851283-85d8c498-a536-4e42-8d15-9e4bef4554dc.png">     </br>
An eaxmple of a single muon pulse (left) and 100 pulses (right) in a fill is shown below:</br>
<img width="892" alt="Screen Shot 2022-04-11 at 7 56 15 PM" src="https://user-images.githubusercontent.com/27436642/162851751-b079ffeb-2392-4a07-abe8-20deb881530f.png">
<p align = center> Fig.2</p>
The final results matched well with small desktop experimental results as shown:</br> 
<p align = "center">
<img width="376" alt="Screen Shot 2022-04-11 at 7 42 42 PM" src="https://user-images.githubusercontent.com/27436642/162850620-556680a7-bd75-43c9-adf4-2eb090b10b27.png">
<img width="417" alt="Screen Shot 2022-04-10 at 6 45 55 AM" src="https://user-images.githubusercontent.com/27436642/162614447-6c587309-45aa-411b-b2c1-192a5fa95a6e.png"></br>
Fig. 3
</p>

------------------------------------------------------------------------------------------------------------------------

### Instructions for the code:
This just uses C++ and ROOT. Make sure to install these. I have briefly described the three important modules that produce the above plots.

***The gainSimulation.C code:***</br>
This a simple C file using root as:</br>
$ root -l gainSimulation.C</br>
This produces pdfs with the following number of muons in a fill:</br>
1, 20, 50, 100, 300, 1000. </br>
Plots of fig 2 are produced running this code with 1 and 100 muons (saved in the pdfs). </br

 ***The laser_gainSimulation.C code:***</br> 
The default gainsimulationTest() function takes 50 cycles, 16 fills and 5 laser shots in each cycle.</br> 
This gives reasonble results. You can run reduce the number of cyles for a quick check or increase it to reproduce fig. above.</br>
$ root -l laser_gainSimulation.C</br>
$ root [1] gainSimulationTest()</br>
Creates a root file *noWiggle_fac10_g320_m100.root* with a drop-factor (p<sub>o</sub>) of 10, gap $tau;<sub>r</sub>= 320 us, and n<sub>o</sub>=100 muons in a fill.</br>
Left panel of fig.3 is the same for 2000 cycles.

***The fit_sipm_wiggle.C code:***
This takes the root file *noWiggle_fac10_g320_m100.root* as fits it using the final formula for simulation gain G(t)/G<sub>o</sub>. Left panel of fig.3 shows this fit.</br>
$ root -l fit_sipm_wiggle.C

------------------------------------------------------------------------------------------------------------------------
