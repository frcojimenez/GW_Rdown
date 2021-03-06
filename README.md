# The project
 
The GW_Rdown project aims at producing accurate ringdown (RD) models that describe the final state of binary black hole mergers. <br>
The RD emerges as the late trail of radiation usually represented as a sum of damped sinusoids and its modelling is essential to test general relativity (GR) in its most extreme regime.  


# Prerequisites

To fully and succesufully use all codes stored here, one needs to have the following codes and versions installed:
* Python ≥ 3.0 along with the qnm package [8-9] and the numpy package;
* Mathematica ≥ 11.0;
* (optional, only if one wishes to use the Fortran version of the KerrQNMOvertones code) A Fortran compiler compatible at least with the Fortran 2003 standard. Compatibility with the Fortran 90 or 95 standards is sufficient provided a few lines of the code are disabled (see the Readme file in the ./codes/Fortran folder).

The KerrQNMOvertones code is available in three versions that respectively require only the use of Fortran, only the use of Python (without the need for the qnm package), or only the use of Mathematica.


# Structure

In this repository one can find: 
* codes to generate the models;
* data describing the quasi-normal mode (QNM) spectrum;
* data resulting from fits informed by numerical relativity (NR) waveforms;

all from reference [7].<br>

**Please refer to the above paper if using these codes and results; please also refer to the original works [2-4] describing the Mathematica code available at [5-6] on which the codes KerrQNMOvertones.nb, .ipynb and .f90 are based, if using the latter codes.**


Structure:
* ./codes :  
  * /Mathematica. It contains two classes of codes, both written in Mathematica. Rdown.nb and Rdown.m: these two (equivalent) codes contain the main functions                                      used to produce the RD models and ansätze (OvertoneModel), to compute the QNM frequencies and damping times (\\[Omega]lmnPy) from                                    the python code qnm [8-9] complemented by our tables for the (l=2,m=2,n=8,9) modes ([7]), or to estimate the final mass and spin
                               from the GW strain (FitRingdownGrid). 
                                   Moreover, it also includes a notebook (KerrQNMOvertones.nb) to compute the QNM spectrum from scratch based on the original code 
                               from [2-4] and available at [5-6]. 
  * /Python. Same version of the KerrQNMOvertones code translated into Python language.
  * /Fortran. Same version of the KerrQNMOvertones code translated into Fortran 90/Fortran 2003 language for a vastly reduced computing time. It is provided along with a Readme file and two example input parameter files (.nml).
                                   
* ./data :    
  * /QNMdata. It contains the QNM data (tables) for the lmn = 228 and lmn = 229 co-rotating and lmn = 22(8-9) counter-rotating modes, as decribed in reference [7].
  * /NRFits. NR fit results obtained by reference [7].                                  

# The data/NRFits folder.

You will find there a list of BBH_SXS_index folders, each one corresponding to one of the 620 SXS NR ringdown fits produced in [7]. In each folder you will find three files:

* mass_spin.dat: It contains the normalised final mass and the final dimensionless spin of the given simulation, read off from the simulation's metadata.
* datafile_nmin0_nmax12.dat: It contains the results of the fits in the following format: mass, spin , mismatch, epsilon, Bayesian Information Criterion (BIC).
* datafile_BFAmplitudes_nmin0_nmax12.mx. It contains the results of the best fit amplitudes for each NR case. You need Mathematica to load the .mx files.

# Description of KerrQNMOvertones.

This code allows for the computation of the (l,m,n) QNM complex frequency ωlmn and angular separation constant Almn for a spin-weigth-s perturbation of a Kerr (or Schwarzschild) black hole using Leaver's method [1].

The code is based on the original Mathematica code from [2-4] available at [5-6]; please do refer to these original works and note that the code available at [5-6] also provides a computation of the angular and radial wavefunctions of the QNMs (not included here).

The specificities of this version are described in [7] and include:
* the use of Leaver's inversions [1] in the calculation of the continued fractions for a more stable recovery of overtones;
* much larger numbers of steps in the approximations to the continued fractions allowed by the direct implementation of the secant method instead of Mathematica's
       memory-consuming built-in root-finding algorithm;
* and a progressive increase of this number of steps over successive iterations along with a convergence criterion, that can be of use for modes for which Leaver's method is less efficient such as near the algebraically special Schwarzschild mode ω = - 2 *i*.
   
The units convention also differs from [1-6] where the mass and maximal dimensionless spin are set to 1/2: in KerrQNMOvertones, units are set such that the mass of the black hole is 1. <br>
This version requires two close but distinct initial guesses for ωlmn.

Note that this version has *not* been extensively tested for s ≠ -2 nor for l ≠ 2.
    
    
    
# References

* [1] E. W. Leaver, "An analytic representation for the quasi-normal modes of Kerr black holes".  Proc. R. Soc. Lond. A **402**:285 (1985).

* [2] E. Berti, V. Cardoso and C. M. Will, "On gravitational-wave spectroscopy of massive black holes with the space interferometer LISA".  Phys. Rev. D **73**:064030 (2006).  [arXiv: gr-qc/0512160]

* [3] E. Berti, V. Cardoso and M. Casals, "Eigenvalues and eigenfunctions of spin-weighted spheroidal harmonics in four and higher dimensions".  Phys. Rev. D **73**:024013 (2006);  Erratum -- ibid. **73**:109902 (2006).  [arXiv: gr-qc/0511111]

* [4] E. Berti, V. Cardoso and A. O. Starinets, "TOPICAL REVIEW: Quasinormal modes of black holes and black branes".  Class. Quant. Grav. **26**:163001 (2009).  [arXiv: 0905.2975]

* [5] V. Cardoso, https://centra.tecnico.ulisboa.pt/network/grit/files/ringdown/

* [6] E. Berti, https://pages.jh.edu/~eberti2/ringdown/

* [7] F. Jiménez Forteza and P. Mourier, "High-overtone fits to numerical relativity ringdowns: beyond the dismissed n=8 special tone".  (2021).  [arXiv: 2107.11829]

* [8] L. C. Stein, "qnm: A Python package for calculating Kerr quasinormal modes, separation constants, and spherical-spheroidal mixing coefficients",  J. Open Source Softw. **4**:1683 (2019). [arXiv:1908.10377]

* [9] https://pypi.org/project/qnm/
