# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:50:12 2021

@author: XPS
"""


import numpy as np
Network_type = 'BA' # The options are: 'WS'/ 'ER'/ 'BA'/ 'Villages'

# Initialising the basic networks
N_init = 1000 # population size
I_0 = 10 # initial infectious number of ppl (they need to show symptoms to be discovered at first, so they are already infectious)
Reff_range = 3 # period of taking average of the secondary infections for calculating R_eff
m_BA = 5 # Number of edges to attach from a new node to existing nodes in BA

# Initialising the villages - pop size different from the above
m_BA = 5 # input for BA group generations
pop = [int(i) for i in np.array([1607000,2869200,999000,1644400,1882800])/200] # Population of the villages (sub-regions)
n_connect_villages = int(np.mean(pop)/20) # Number of connections between each villages (currently signed as a fourth of the population)

# 26th march is the time first lockdown comes into force 
# First case at 11st Feb.
# Data of 91 infection is at 6th March (assumen scaled 100 times down)
end_day = 18+25 # days from the first case
start_day = 18+6 # shift the real data to the left to the number of scaled cases
N_real = 8.982e6 # London population


# Infection rate distribution for edges - beta: infection rate
beta_mu, beta_sigma = 0.2, 0.1 # mean and standard deviation
# Other rates (recovery, mortality, reinfection) are replaced by distributions, but can still be calculated (will be printed in the console)
# This means that we can get individual patterns in control by altering the distributions, and excempted ppl from not moving into the next stage forever


# Contact tracing: threshold to be quarantine if their neighbor are tested positive (weighted_beta > contact_thre)
contact_thre = 1/2*beta_mu
# Remained number of contacts if ppl are self-isolated at home
remain_contact = 3


# Recovery vs. mortality ratio:
ratio_r, ratio_m = 98.655, 1.345 # London
# RECOVERY PERIOD
mu_r, sigma_r = 18, 12.6 # Western Ethiopia
# MORTALITY PERIOD
mu_m, sigma_m = 13, 10.37 # London


# LATENT PERIOD (China)
mu_l, sigma_l = 5.5, 3.5 
# INCUBATION PERIOD (Canada)
mu_ib, sigma_ib = 6.74, 2.42 
# HOSPITALISATION PERIOD
# Hospitalization should change both the recovery & death rate. The recovery period shortened and death period lengthened.
# It should not be available for everyone. Need discussion
mu_h, sigma_h = 19, 12.5
# REINFECTION PERIOD (global)
mu_ri, sigma_ri = 60, 10


# Assigning colors for the network plots
colors = {"S":"lightblue","E":"orange","I":"red","R":"green","D" : "black"}


r_fraction = ratio_r / (ratio_r + ratio_m) # fraction of recovering ppl from pandemic
m_fraction = 1 - r_fraction # fraction of mortality ppl from pandemic
# gamma: recovery rate
# mortality: mortality rate
gamma = 1/mu_r * r_fraction
mortality = 1/mu_m * m_fraction
print('Strong ppl percentage (will recover): ', r_fraction*100, '%')
print('Vulnerable ppl percentage (will die): ', m_fraction*100, '%')
print('Recovery rate (gamma) = ', gamma)
print('Mortality rate (mu) = ', mortality)