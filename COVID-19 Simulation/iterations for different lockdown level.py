# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 09:48:50 2021

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

# Start transmission iterations
days = 50 # days of iteration

# The current isolation method is that the least weighted connections will be cut off, except the first several ppl (number set below). The connections will be reestablished when they are recovered
if_isolate =  0 # Chosing whether to stay at home or not when ppl get infected - will be very slow
test_ct_quar = 0 # Testing, contact tracing and quarantine of close contacts


# Simulating lockdown
for node in network:
    isolate(network, node, remain_contact)
    network.nodes[node]["isolate"] = True

for remain_contact in [3,5,7,9]:
    print('Remained contact:',remain_contact)
    ig_net = ig.load(r"testnetwork.gml")
    network = ig_net.to_networkx()
    N = nx.Graph.number_of_nodes(network)
    network = copy.deepcopy(network_init)
    
    # Refresh removed edges every iteration
    S = [N-I_0]; E=[0]; I = [I_0]; R =[0]; D=[0]; R_eff=[] # initialize the number of each group 
    
    for node in network:
        network.nodes[node]['removededges'] = [] # contains beta values
    
    for day in range(days):
        network,sec_inf = transmission(network, N, if_isolate, test_ct_quar, beta_mu, contact_thre, remain_contact)
        
        R_eff.append(sec_inf)
        x,y,z,p,q = count(network)
        I.append(z)
    x = list((range(0,days + 1)))
    plt.plot(x,I,label='Remained contacts = %d'%(remain_contact))
    
plt.xlabel('Days')
plt.ylabel('Population')
plt.title("Mean infection rate=%r  Mean recovery rate=%.5f mortality rate=%.5f \n\
          Population=%d Latent period=%d"
          %(beta_mu, gamma, mortality, N, mu_l) )
plt.legend()
plt.grid()
plt.show()


# Stop timing here
t_stop = timeit.default_timer()
print("time =",(t_stop-t_start)/60, 'min')

