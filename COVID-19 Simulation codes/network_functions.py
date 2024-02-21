# -*- coding: utf-8 -*-

import numpy as np
import copy
import networkx as nx
import random
from itertools import combinations


# Note: we need to make sure every node has an edge connected to it! Lonely ppl disturbs the metrics and adds on simulation time
def basic_networks(model_type, init_I, N, m): # returns a network. Input: network type, initialised infected number of ppl
    if model_type == 'WS':
        n = N    #Number of nodes
        k = 10      #Number of neighbors
        p = 0.8   #probability of rewiring each edge
        init = nx.watts_strogatz_graph(n, k, p)
        for i in range(N):
            init.nodes[i]['state'] = 'S'
            init.nodes[i]['secinf'] = 0
            init.nodes[i]['Iday'] = 0
            init.nodes[i]['Rday'] = 0
        L = np.random.randint(0, N, size=init_I) # randomly select a init_I number of infectious ppl
        for j in L:
            init.nodes[j]['state'] = 'I' 
        
    
    if model_type == 'ER':
        n = N
        p = 0.03
        init = nx.erdos_renyi_graph(n,p)
        for i in range(N):
            init.nodes[i]['state'] = 'S'
            init.nodes[i]['secinf'] = 0
            init.nodes[i]['Iday'] = 0
            init.nodes[i]['Rday'] = 0
        L = np.random.randint(0, N, size=init_I) # randomly select a init_I number of infectious ppl
        for j in L:
            init.nodes[j]['state'] = 'I' 
    
    if model_type == 'BA':
        seed = None # Seed for random generator (optional)
        init = nx.barabasi_albert_graph(N,m)
        for i in range(N):
            init.nodes[i]['state'] = 'S'
            init.nodes[i]['secinf'] = 0
            init.nodes[i]['Iday'] = 0
            init.nodes[i]['Rday'] = 0
        L = np.random.randint(0, N, size=init_I) # randomly select a init_I number of infectious ppl
        for j in L:
            init.nodes[j]['state'] = 'I' 
    
    return init



def network_villiges(init_I, pop, base_contact, n_connect_villages):
    Network_list = []
    node_group_list = []
    first_node = 0
    for i in range(len(pop)):       # generate "villages" by looping through each village population
        n_pop = pop[i]
        G = nx.barabasi_albert_graph(n_pop, base_contact)
        Network_list.append(G)
        last_node = first_node + n_pop
        node_group_list.append(list(range(first_node, last_node)))
        first_node = last_node
    
    init = nx.disjoint_union_all(Network_list)
    combs = combinations(node_group_list,2)
    
    for comb in list(combs): # Going through each combination of the village node numbers
        for i in range(n_connect_villages): # Connect each village with n_connect_group different new edges
            u = np.random.choice(comb[0]) # Pick one end of the new edge randomly in the first group
            v_list = comb[1]
            u_neighbors = np.array(list(init.edges(u)))[:,1]
            # print(v_list)
            # print(u_neighbors)
            v_legit = []
            for v in v_list: # Pick the other end of the new edge in the second group that is not previously connected already
                if v not in u_neighbors:
                    v_legit.append(v)
                else:
                    pass
                    # print('Someone is already connected! - ',v)
            v = np.random.choice(v_legit)
            # print(u,v)
            init.add_edge(u,v)
    
    N_c = len(init)
    # print('Population is ' + str(N_c))
    
    for i in range(N_c):
        init.nodes[i]['state'] = 'S'
        init.nodes[i]['secinf'] = 0
        init.nodes[i]['Iday'] = 0
        init.nodes[i]['Rday'] = 0
    L = np.random.randint(0, N_c, size=init_I) # randomly select a init_I number of infectious ppl
    for j in L:
        init.nodes[j]['state'] = 'I' 
    
    return init



# Note: Some values are lower han zero, but that is fine because they make the distribution normal and we can see them as zeros really
# Normal distribution for all parameters, different for every node. Input: network, population size and other parameters
def normal_dist(G, n, r_fraction, beta_mu, beta_sigma, mu_r, sigma_r, mu_m, sigma_m, mu_l, sigma_l, mu_ib, sigma_ib, mu_h, sigma_h, mu_ri, sigma_ri):
    
    LT_dist = []
    IB_dist = []
    
    # LATENT PERIOD
    s_l = np.random.normal(mu_l, sigma_l, n)
    # INCUBATION PERIOD
    s_ib = np.random.normal(mu_l, sigma_l, n) + np.random.normal(mu_ib-mu_l, sigma_ib+sigma_l, n) 
    # HOSPITALISATION
    s_h = np.random.normal(mu_h, sigma_h, n)
    # REINFECTION
    s_ri = np.random.normal(mu_ri, sigma_ri, n)
    
    # RECOVERY PERIOD
    s_r = np.random.normal(mu_r, sigma_r, n)
    # MORTALITY PERIOD
    s_m = np.random.normal(mu_m, sigma_m, n)
    
    
    for node in G: # Adding weighted period lengths to each node
        G.nodes[node]['dayLT'] = s_l[node]
        G.nodes[node]['dayIB'] = s_ib[node]
        IB_dist.append(G.nodes[node]['dayIB'])
        LT_dist.append(G.nodes[node]['dayLT'])
        G.nodes[node]['dayHOS'] = s_h[node]
        G.nodes[node]['dayReinf'] = s_ri[node]
        
        # if if_isolate == True or test_ct_quar == True:
        G.nodes[node]["isolate"] = False # self-quarantine
        G.nodes[node]['iso14day'] = 0 # Our model starts the quarantine (forced by contact tracing) on the day 0
        G.nodes[node]["quarantine"] = False # Forced complete quarantine
        # else:
        #     pass
        
        # distribute recover and mortality periods according to everyone's fate:
        x = random.random()
        if x < r_fraction:
            G.nodes[node]['fate'] = True
            G.nodes[node]['dayRC'] = s_r[node]
            G.nodes[node]['dayD'] = 0
        else:
            G.nodes[node]['fate'] = False
            G.nodes[node]['dayD'] = s_m[node]
            G.nodes[node]['dayRC'] = 0
        
    
    # Adding weighted infection rate attribute to each edge
    n_edges = G.number_of_edges()
    w = np.random.normal(beta_mu, beta_sigma, n_edges)
    k = 0
    for i, j in G.edges():
        G[i][j]['weightedbeta'] = w[k]
        k+=1
    
    return G, LT_dist, IB_dist


