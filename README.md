# Numerical analysis of gain fluctuation in a SiPM
This study derives how the gain fluctuates in the SiPM used for reading the electron calorimeters in the Muon g-2 experiment. The numerical solutions 
for this gain is then simulated and finally compared with 

## Gain variation in a fill: semi-analytical approach
The distribution of the electron/positron as a function of electron energy E (i.e. the wiggle plot), decay asymmetry A and a constant phase term φ that depends on the drift time of the electrons is given by,
<p align = "center">
<img width="546" alt="Screen Shot 2022-04-11 at 1 54 55 PM" src="https://user-images.githubusercontent.com/27436642/162800199-af362a3b-5794-42d7-bf4a-0e6cc1cfcf80.png"></br>
 Fig. 1
</p>
## Simulation of gain function
Silicon Photo-Multipliers (SiPMs) were used to read data from the calorimeters and the gain fluctuation (or function) was modelled based on 
bias-voltage sagging (due to decay of current in capacitors and other circuit elements) - this matched well with small desktop experimental 
results as shown:</br> 
<img width="348" alt="Screen Shot 2022-04-10 at 6 41 23 AM" src="https://user-images.githubusercontent.com/27436642/162614319-5cb05518-582e-4437-9447-fe8fc9adaefb.png">
<img width="437" alt="Screen Shot 2022-04-10 at 6 45 55 AM" src="https://user-images.githubusercontent.com/27436642/162614447-6c587309-45aa-411b-b2c1-192a5fa95a6e.png"></br>

**The laser_gainSimulation.C code:**
The default gainsimulationTest() function takes 2000 cycles, 16 fills and 2 laser shots in each cycle.</br> 
This gives reasonble results. You can run reduce the number of cyles for a quick check. </br>
$ root -l laser_gainSimulation.C</br>
$ root [1] gainSimulationTest()</br>

**The gainSimulation.C code:**
This produces pdfs with the following number of muons in a fill:</br>
1, 20, 50, 100, 300, 1000. </br>
Fig 1 and 2 are produced running this code with 1 and 100 muons (saved in the pdfs). </br>

**The fit_simp_wiggle.C code:**
varyDrop10_m100_c500_f16_80_5.root
