# The project
 
The GW_Rdown project aims to produce accurate ringdown (RD) models that describe the final state of binary black hole (BBH) mergers. 
The BH RD emerges as the late trail of radiation usually represented as a sum of damped sinusoids and their modelling is essential to test general relativity (GR) in its most extreme regime.  

# Prerequisites

To fully and succesufully use all codes stored here, one needs to have the following codes and versions installed:
* Python >= 3.0
* Mathematica >= 11.0
* The qnm package [7].

# Structure

In this repository one can find: 
* codes to generate the models.
* data describing the QNM spectrum.
* data resulting from fits informed by numerical relativiy (NR) waveforms. 

Structure:
* ./codes :  
  * Mathematica. It contains two classes of codes, both written in Mathematica. Rdown.nb and Rdown.m: these two codes contain the main functions                                      used to produce the RD models and ansätze (OvertoneModel), to compute the QNM frequencies and damping times (\[Omega]lmnPy) from                                    the python code [7] or to estimate the final mass and spin from the GW strain (FitRingdownGrid). 
                                   Moreover, it also includes a code (KerrQNMOvertones) to compute the QNM spectrum from scratch based on the original code from 
                                   [2-4] and available at [5]. 
  * Python. Same version of the KerrQNMOvertones code in python language.
                                   
* ./data :    
  * ./QNMdata. It contains the QNM data (tables) for the co-rotating lmn = 228 and lmn = 229 and 22(8-9) counter-rotating modes.
  * ./NRFits. NR fit results obtained by reference [6].                                  

# Description of KerrQNMOvertones.

This code is based on the original code from [2 - 4] available at [5]; please do refer to these original work and note that the code available at [5] also provides a computation of the angular and radia wavefunctions of the QNMs (not included here).
The specificities of this version are described in [6] and include:
* the use of Leaver' s inversions [1] in the calculation of the continued fractions for a more stable recovery of overtones
* much larger numbers of steps in the approximations to the continued fractions allowed by the direct implementation of the secant method instead of Mathematica' s
       memory - consuming built - in root - finding algorithm
* and a progressive increase of this number of  steps over successive iterations along with a convergence criterion, that can used for modes for which Leaver' s method is less efficient such as near the algebraically special Schwarzschild mode \[Omega] = - 2 I.
   
The units convention al so differs from [1 - 5] where th e mass an d maximal dimensionless spin ar e se t to 1/2. This version requires two close but distinct initial guesses for \[Omega]lmn. Units are set such that the mass of the black hole is 1.
Note that this version has not been extensively tested for s != -2 nor for l != 2
                                   
# References
* [1] E. W. Leaver, "An analytic representation for the quasi-normal modes of Kerr black holes".  Proc. R. Soc. Lond. A 402: 285 (1985).

* [2] E. Berti, V. Cardoso and C. M. Will, "On gravitational-wave spectroscopy of massive black holes with the space interferometer LISA".  Phys. Rev. D 73:064030 (2006).  [arXiv: gr-qc/0512160]

* [3] E. Berti, V. Cardoso and M. Casals, "Eigenvalues and eigenfunctions of spin-weighted spheroidal harmonics in four and higher dimensions".  Phys. Rev. D 73:024013 (2006);  Erratum -- ibid. 73:109902 (2006).
	[arXiv: gr-qc/0511111]

* [4] E. Berti, V. Cardoso and A. O. Starinets, "TOPICAL REVIEW: Quasinormal modes of black holes and black branes".  Class. Quant. Grav. 26:163001 (2009).  [arXiv: 0905.2975]

* [5] E. Berti, https://pages.jh.edu/~eberti2/ringdown/

* [6] F. Jiménez Forteza and P. Mourier, "High-overtone fits to numerical relativity ringdowns: beyond the dismissed n=8 special tone".  (2021).  [arXiv: 2107.11829]

* [7] https://pypi.org/project/qnm/ 

             
