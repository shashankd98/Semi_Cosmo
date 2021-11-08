This work studies the cosmological evolution of baryonic gas in an idealized isolated halo. 

Authors: Shashank Dattathri (shashank.dattathri@gmail.com) and Prateek Sharma (prateek@iisc.ac.in). 

Some of the codes used to generate the figures from the paper are uploaded here. The full data will be shared on reasonable request to the authors. The simulations were run using PLUTO v4.4.  

Description of codes: <br />
/Codes/DK14_density_potential.py : creates the density-potential pair for a DK14 density profile.  <br />
/Codes/mass_accretion_histories.py : creates the mass accretion history profiles for the dark matter halos.  <br />
/Codes/plane.py : plots the fundamental plane <br />

PLUTO files: <br />
/PLUTO/Non_rad : Non-radiative run for M_0 = 10^{14} M_\odot <br />
/PLUTO/Cooling_flow : Cooling flow run for M_0 = 10^{14} M_\odot <br />
/PLUTO/AGN_feedback : AGN feedback run for M_0=10^{14} M_\odot, \epsilon = 10^{-4} <br />

Data files: <br />
/Data/1e12_1em5/: data for the M0=1e12 Msun, epsilon=1e-5 run (figure 8) <br />
/Data/1e14_1em4/: data for the M0=1e14 Msun, epsilon=1e-4 run (figure 9) <br />
/Data/Baryon_fraction/: data for baryon fraction for various AGN feedbback runs (figure 15) <br />
/Data/Cycles/: data for cool core cycles for the fiducial feedback run (figure 11) <br />
/Data/Inflow_accretion/: data for mass inflow rate (Mdot_in) and mass accretion rate (Mdot_acc) for different efficiencies and halo masses (figures 12 and 13) <br />
/Data/Self_similar/: data for analyzing the self-similarity of density profiles for non-radiative and AGN feedback runs (figures 5 and 16 respectively) <br />

Data analysis (pyPLUTO) files:  <br />
(the variable wdir should be pointed to the working directory) <br />
/Visualization/coolcores.py : plots the time-averaged density and standard deviation for the AGN feedback run.  <br />
/Visualization/baryon_fraction.py : calculates the baryon fraction <br />
/Visualization/snapshots.py : plots snapshots for the 2D runs <br />
/Visualization/visualize.py : plots relavent quantitites (density, pressure, velocity, temperature, t_cool/t_ff, Mdot) for 1D runs <br />
/Visualization/visualize2d.py : plots relavent quantitites (density, pressure, velocity, temperature, t_cool/t_ff, Mdot) for 2D runs <br />
