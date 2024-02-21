# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:45:01 2021

@author: XPS
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
import networkx as nx
import random
import matplotlib.animation as animation
from IPython.display import HTML,display
from itertools import combinations
from transmission_functions import *
from attributes import *
import timeit
import pandas as pd
import igraph as ig

t_start = timeit.default_timer()
ig_net = ig.load(r"testnetwork.gml")
network = ig_net.to_networkx()
N = nx.Graph.number_of_nodes(network)


# Start transmission iterations
days = 100 # days of iteration

# The current isolation method is that the least weighted connections will be cut off, except the first several ppl (number set below). The connections will be reestablished when they are recovered
if_isolate =  0 # Chosing whether to stay at home or not when ppl get infected - will be very slow
test_ct_quar = 0 # Testing, contact tracing and quarantine of close contacts
daily_net = 0 # Plots the network everyday

pos = nx.spring_layout(network)
S = [N-I_0]; E=[0]; I = [I_0]; R =[0]; D=[0]; R_eff=[] # initialize the number of each group 

if daily_net == True:
    network_plot(network, N, colors, pos)
    plt.figure(0,figsize=(10,10))
    plt.tight_layout()
    plt.axis('off')
    plt.title('Day 0 (initial network)')
    plt.show()

# Refresh removed edges every iteration
for node in network:
    network.nodes[node]['removededges'] = [] # contains beta values

for day in range(days):
    print('day',day+1)
    network,sec_inf = transmission(network, N, if_isolate, test_ct_quar, beta_mu, contact_thre, remain_contact)
    if daily_net == True:
        # Plotting the network for each day
        network_plot(network, N, colors, pos)
        plt.figure(day,figsize=(10,10))
        plt.tight_layout()
        plt.axis('off')
        plt.title('Day %d' %(day+1))
        plt.show()
    
    R_eff.append(sec_inf)
    x,y,z,p,q = count(network)
    S.append(x)
    E.append(y)
    I.append(z)
    R.append(p)
    D.append(q)

# Stop timing here
t_stop = timeit.default_timer()
print("time =",(t_stop-t_start)/60, 'min')
# Case fatality rate (possibility to die if an individual got infected)
CFR = D[-1]/(N-S[-1])
print('CFR =',CFR)
print('Required immuned portion of people for herd inmmunity from theory:', 1-1/R_0)


#%% Plotting the cumulative changes of all the groups
x = list((range(0,days + 1)))
#x = list((range(0,11)))
plt.figure(2)
plt.plot(x,S,label='Susceptible')
plt.plot(x,E,label='Exposed (infected but not infectious)')
plt.plot(x,I,label='Infectious')
plt.plot(x,R,label='Recovered')
plt.plot(x,D,label='Dead')
plt.xlabel('Days')
plt.ylabel('Population')
plt.title("Population=%d" %(N) )
plt.legend()
plt.grid()
plt.show()


#%% Plotting against the London data
n_day = end_day - start_day # Number of days plotted

# Reading the real data
df = pd.read_csv("data/data_2021-Nov-07.csv", usecols=["cumCasesBySpecimenDate"])
data = df[-end_day : -start_day]
x = list((range(0,n_day)))

plt.figure(3)
plt.plot(x, np.flip(np.array(data)) * N/N_real, 'black', alpha=1, lw=2, label='Real data (London before lockdown)')
#plt.plot(x,S,label='Susceptible')
#plt.plot(x,E,label='Exposed (infected but not infectious)')
plt.plot(x,I[:n_day],label='Infectious')
plt.plot(x,R[:n_day],label='Recovered')
plt.plot(x,D[:n_day],label='Dead')
plt.xlabel('Days')
plt.ylabel('Population')
plt.title("Population=%d" %(N) )
plt.yscale("log")
plt.legend()
plt.grid()
plt.show()


#%% Calculating the R_eff values averaged in Reff_range number of days
iter_num = int(days/Reff_range)
Reff = []
Reff_err = []
for i in range(iter_num):
    l_lim = i*Reff_range
    h_lim = (i+1)*Reff_range
    denom = 0
    numer = 0
    for j in range(l_lim,h_lim):
        denom += len(R_eff[j])
        numer += sum(R_eff[j])
    if denom == 0:
        final_r = 0
    else:
        final_r = numer/denom
    Reff_err.append(np.sqrt(final_r))
    Reff.append(final_r)

Reff_err = np.array(Reff_err)
Reff = np.array(Reff)
Reff_x = np.linspace(0,iter_num-1,iter_num)*Reff_range
plt.figure(4)
plt.plot(Reff_x,Reff)
#plt.errorbar(Reff_x,Reff,Reff_err,ls='none')
plt.fill_between(Reff_x, Reff - Reff_err, Reff + Reff_err,
                 color='gray', alpha=0.2)
plt.xlabel('Day')
plt.ylabel('Effective reproductive number R_eff')
plt.title("Total population=%d  Calculated every %d days "%(N,Reff_range))
plt.show()


#%% Growth rate of infectious ppl

growth_rate = np.exp(np.diff(np.log(I))) - 1 # exponential growth rate of cases per day (lambda)
plt.figure(5)
plt.plot(x[:-1],growth_rate)
plt.xlabel('Day')
plt.ylabel('Growth rate r')
plt.title("Total population=%d"%(N))
plt.show()

# doubling_time = [np.log(2)/r for r in growth_rate]
# plt.figure(4)
# plt.plot(x[:-1],doubling_time)
# plt.xlabel('Day')
# plt.ylabel('Doubling time')
# plt.title("R_0=%r  Total population=%d"%(R_0,N))
# plt.show()