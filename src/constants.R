# ----- parameters to infer -----
# variable parameters
r_0 = 1
J = 0.42
c = 0

# ----- hyperparameters -----
### ARRAY DYNAMICS SIMULATION ###
# input:         -
# output:       constants
# description:  list all constants used for simulation
# author:       HR

# number of iterations
TIME = 2e03
dt = 1  # data sampling rate
X = 50  # matrix dimension

cT = TIME/2  # time after which ligand is spiked in

epsilon = 1e-06  # some noise

# conversion factor between FRET ratio and activity
lambda = 0.1

# K_on and K_off - Asp
# # Vladimirov et al., PLOS Comp. Biol. 2008
# Kon = 12  # uM
# Koff = 1.7  # uM

# MeAsp
Kon = 3000  # uM
Koff = 20  # uM

# express concentration in units of Koff
# simplification for Kon >> Koff
# Kon = Inf
# Koff = 1

# methylation
Mtot = 48  # total number of methylation sites per CSU
k0 = 1
k1 = 1/12

# reaction constants of methylation
# Vladimirov et al., PLOS Comp. Biol. 2008
kR = 0.0182
kB = 0.0364

# # Boltzmann constant and temperature
# k = 1.38064852e-23
# Temp = 293.15

