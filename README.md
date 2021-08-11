# The project
 
The GW_Rdown project aims to produce accurate ringdown (RD) models that describe the final state of binary black hole (BBH) mergers. 
The BH RD emerges as the late trail of radiation usually represented as a sum of damped sinusoids and their modelling is essential to produce tests of general relativity GR in its most extreme regime.  

# Structure

In this repository one can find: 
* codes to generate the models 
* data describing the QNM spectrum 
* data resulting from fits informed by numerical relativiy (NR) waveforms. 

Structure:
* ./codes :  -- Mathematica. It contains two classes of codes, both written in Mathematica. Rdown.nb and Rdown.m: these two codes contain the main functions                                      used to produce the RD models and ansätze (OvertoneModel), to compute the QNM frequencies and damping times (\[Omega]lmnPy) from                                    the python code [1] or to estimate the final mass and spin from the GW strain (FitRingdownGrid). 
                                   Moreover, it also includes a code (KerrQNMOvertones) to compute the QNM spectrum from scratch based on the original code from 
                                   [2-4] and available at [5]. 
* ./data :    -- /QNMdata. It contains the QNM data (tables) for the co-rotating lmn = 228 and lmn = 229 and 22(8-9) counter-rotating modes.
                   -- /NRFits. NR fit results obtained by reference [6].                                  
       
                                   

             
