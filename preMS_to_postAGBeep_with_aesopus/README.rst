********************************
preMS_to_postAGBeep_with_aesopus
********************************

The inlists in this test suite case are designed to evolve a 2.4Msun model from the zero age main sequence to the post asymptotic giant branch (AGB). We use the equal evolutuionary point (EEP) definition of Dotter (2016), whereby a model enters the post AGB when its H rich envelope comprises 20% of the current total stellar mass. 

Running the test suite case as is, i.e.:

./clean; ./mk; ./rn

will evolve the model from the 41st thermal pulse to the post-AGB EEP. This takes approximately X minutes on 6 cores. To reproduce the entire run, comment out the 'pre_test_suite_case_termination' condition in the extras_finish_step subroutine of src/run_star_extras.f90. Copy the script /docs/run_prems_eep into this directory and run using ./run_prems_eep (~X days on 6 cores).

This directory contains low temperature opacity tables built using AESOPUS 2.0 (see Marigo, Aringer, Girardi and Bressan 2022, also Marigo and Aringer 2009), built for reference metallicities of Z=0.0134,  Z=0.004. Tables custom to the metallicity regime being modelled are important for AGBs, given the significant increase in opacity they can undergo due to the third dredge up. The settings for these tables and method for implemention into MESA are described at the end of this file, and the sixth instrument paper: Jermyn et al. (2023). 

A successful run is indicated by: 

-The termination code 't_extras_finish_step' printed to your terminal. Final stats are also printed: tp-agb duration, change in h-exhausted core mass and the intial to final isotopic ratios: c12/c13, n14/n15, c/o. These are interesting to look at on the TP-AGB as indicators of third dredge up and hot bottom burning (in intermediate mass stars). A solar scaled model with start ZAMS evolution with c/o ~ 0.5. This will reduce to c/o ~ 0.3 after the first dredge up on the red giant branch. Further increases in the c/o ratio then on the TP-AGB phase indicate the model has undergone some third dredge up, where primary c12 has been mixed into the envelope. Also indicated by an increasing c12/c13 ratio. Hot bottom burning occurs in higher mass AGBs (~4-5Msun, depending on Z), where temperatures at the base of the convective envelope are high enough for CNO burning. In this case, the surface n14/n15 ratio increases, both c12/c13 and c/o will decrease.

-Output file, `agb_stats.dat', printed during AGB evolution. This contains information about individual thermal pulses (core mass, stellar age, surface isotopic ratios).  

-Plotting window comparing the evolution of the current test suite case model with the model in /docs/freedman_run (same inlists, but freedman low temperature opacities). Shown is the difference between the surface opacity and isotopic ratios during the AGB. 

********************************

From TP=41 to EEP: X hours on 6 cores. From preMS to EEP: X hours on 6 cores. 

This test suite case is based on work developed for Cinquegrana, Joyce and Karakas (2022, 2023). Our timestep and solver controls are modified from the work of R. Farmer (see Farmer et al. (2016), Laplace et al. (2021) and the test suite case: 12M_pre_ms_to_core_collapse). We based our thermal pulse counter off the c13_pocket test suite case.  

AESOPUS in MESA: 

1. Build tables at: http://stev.oapd.inaf.it/cgi-bin/aesopus. The settings we used are as follows: X. 

2. Save each table in a text file and copy into $MESA_DIR/kap/preprocessor/AESOPUS. The README file within that directory describes how to convert your tables into a .h5 file. 

3. Copy .h5 file into the current working directory. 

Last-Updated: October-2023 (mesa r23051) by G. Cinquegrana
