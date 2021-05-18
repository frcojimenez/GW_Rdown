The folders

    ./data/BL-5-CSVs
    ./data/ratio1-CSVs
    ./data/ratio3-CSVs

contain comma-separated-value (CSV) files with columns for
1. simulation time
2. horizon area
3. mode 0
4. mode 1, ...


The quantities are the shear, multipoles, and the "xi scalar" on the horizons:

    AH - the apparent horizon which relaxes asymptotically to the final BH horizon
    inner - inner common horizon
    top - S1, the (smaller) horizon on larger z-values
    bot - S2, the (larger) horizon on smaller z-values


There are files for multiple resolutions:

    1/h = 60  extending up to T = 50/1.3 = 38.46153846153846
    1/h = 120 extending up to T = 50/1.3 = 38.46153846153846
    1/h = 180 extending up to T = 50/1.3 = 38.46153846153846
    1/h = 240 extending up to T = 50/1.3 = 38.46153846153846
    1/h = 480 extending up to T = 20/1.3 = 15.384615384615383
    1/h = 960 extending up to T = 7/1.3  = 5.384615384615384

All values and times are scaled to an `M_ADM=1` simulation. A comparison of
the different numerical results can be found in
`BL-5-compare_resolutions.pdf`, from which we can infer that the `res=240`
results are reliable up to about `l=12` for shear/multipoles and probably
`l=8` for the xi scalar (c.f. rows 26, 56, and 85 in the PDF).


In the scaled times, the common horizon forms at:

    T_bifurcate = 1.05739


The start time for late time fits in arXiv:2006.03940 is:

    T0 = 20/1.3 = 15.384615384615383
# GW_Rdown
