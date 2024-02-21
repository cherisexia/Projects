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
iter_num = 20 # Number of iterations for plotting the errors

# The current isolation method is that the least weighted connections will be cut off, except the first several ppl (number set below). The connections will be reestablished when they are recovered
if_isolate =  0 # Chosing whether to stay at home or not when ppl get infected - will be very slow
test_ct_quar = 0 # Testing, contact tracing and quarantine of close contacts


S_list = np.zeros([iter_num, days+1])
E_list = np.zeros([iter_num, days+1])
I_list = np.zeros([iter_num, days+1])
R_list = np.zeros([iter_num, days+1])
D_list = np.zeros([iter_num, days+1])
for i in range(iter_num):
    print('Iteration time:',i,'/',iter_num)
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
        S.append(x)
        E.append(y)
        I.append(z)
        R.append(p)
        D.append(q)
    S_list[i,:] = S
    E_list[i,:] = E
    I_list[i,:] = I
    R_list[i,:] = R
    D_list[i,:] = D

ss = np.zeros([3,days+1])
ee = np.zeros([3,days+1])
ii = np.zeros([3,days+1])
rr = np.zeros([3,days+1])
dd = np.zeros([3,days+1])

matrix_list = [ss,ee,ii,rr,dd]
all_list = [S_list,E_list,I_list,R_list,D_list]

for i in range(5):
    matrix = matrix_list[i]
    raw_list = all_list[i]
    for day in range(days+1):
        matrix[0,day] = sum(raw_list[:,day])/iter_num
        matrix[1,day] = min(raw_list[:,day])
        matrix[2,day] = max(raw_list[:,day])
        
    

x_list = list((range(0,days + 1)))
labels = ['Susceptible','Exposed','Infectious','Recovered','Dead']
colors = ["tab:blue","tab:orange","tab:red","tab:green","tab:brown"]
for i in range(5):
    matrix = matrix_list[i]
    plt.plot(x_list, matrix[0,:],label=labels[i],color=colors[i])
    plt.fill_between(x_list, matrix[1,:], matrix[2,:], color=colors[i], alpha=0.2)

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

