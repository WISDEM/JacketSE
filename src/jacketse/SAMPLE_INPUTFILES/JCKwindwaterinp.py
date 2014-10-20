#WATER INFO FOR LOAD CALCULATIONS: SEE MPUTILIZATION
#WIND WAVE LOAD INFORMATION just for load calculation and plot of
#utilization at the end.. Adapted from Matlab's version of windwaterinp.m input

#THIS IS FOR 5MW REFERENCE MONOPILE/or JACKET IN 30 m Water Depth SIte for Offshore 2013
# water   -structure contraining the following
water.wlevel=30.;# water surface z from base of structure which may be below seabed [m]:IT MAY CHANGE
water.wdpth=30.;# water depth [m]
water.T=12; #[s] 50-yr wave period, minimum period for that wave  -Brent gave me the max wave period
water.hw=17.48/2.; #[m] the wave 1/2-height (amplitude not peak-to-peak) 50-yr K-13 deep water site
water.psi=0.; #[deg]  water angle relative to downwind
water.rho=1025;# water density  [kg/m^3]
water.Cd=1.25; # Drag coefficient
water.Cm=2; # Added mass coefficient
# wind     -structure contraining the following
wind.al_shear=0.14;# power law exponent wind from IEC 61400-3 under power production NTM
wind.psi=0.;# [deg] wind angle relative to downwind. positive according to RHR with positive z direction.
wind.rho=0.;#1.225; # density of air  [kg/m^3] 0 if no tower wind load is included (thrust from rotor still on)
wind.U50HH=70.; # [m/s]    50-yr return 3-s gust [m/s] :From IEC Class I
wind.Cd=water.Cd;# Drag Coefficient for tower in air flow
wind.HH=90.; #[m] Hub-height above MSL
#________________________________________________________________________#