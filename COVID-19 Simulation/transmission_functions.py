# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 01:09:05 2021

@author: XPS
"""
import numpy as np
import copy
import networkx as nx
import random


def network_plot(network, N, colors, pos):
    # Plotting the input network
    states= nx.get_node_attributes(network, 'state') 
    color=[colors[states[i]] for i in range(N)]
    edgewidth = [d['weightedbeta']*3 for (u,v,d) in network.edges(data=True)]
    nx.draw_networkx_edges(network, pos, width=edgewidth, edge_color="tab:blue")
    nx.draw_networkx_nodes(network, pos=pos, node_color = color, node_size=10000/N)
    nx.draw_networkx_labels(network, pos=pos, font_color='white', font_size=500/N, font_weight='bold')
    # nx.draw(network, node_color = color, with_labels=True, font_weight='bold')


# Creates secondary infection according to the weighted beta of the edges
def secondary_inf(G, node):
    for neighbor in G.neighbors(node):
        if G.nodes[neighbor]['new_state'] == False: # exempt the new S ppl (turn from R) to be infected on the same day
            x = random.random()
            if G.nodes[neighbor]['state'] == 'S' and x < G[node][neighbor]['weightedbeta']:
                G.nodes[neighbor]['state'] =  'E'
                # if this node was newly infected on this day, then this node should not loop through 'E' logic since Iday = 0
                G.nodes[neighbor]['new_state'] = True
                G.nodes[node]['secinf'] += 1 # added on secondary infection count for each node


# This function defines how ppl would isolate themselves, used in contact tracing and self-isolate after showing symptoms
def isolate(G, node, remain_contact):
    G.nodes[node]["isolate"] = True
    edges = G.edges(node, data=True)
    num_edges = len(edges)
    if num_edges > 0: # Some nodes have zero edges going into it (cut because of other ppl)
        # num = np.random.randint(0,remain_contact+1) # random number of the most connected ppl to be remained (0 to  remain_contact)
        num = remain_contact
        if num_edges > num:
            sorted_edges = sorted(edges, key=lambda t: t[2].get('weightedbeta', 1))
            removededges = sorted_edges[:num_edges-num]
            G.nodes[node]['removededges'] = copy.deepcopy(removededges) # the connections to be removed
            for i in removededges:
                G.edges[i[:2]]['weightedbeta'] = 0 # prevent secondary infection from these edges


# Forced testing + quarantine + contact tracing when infectious ppl get symptoms
def quarantine(G, node, remain_contact, contact_thre):
    # print('Starting the quarantine:',node)
    edges = G.edges(node, data=True)
    G.nodes[node]["quarantine"] = True # prevent the forced quarantined ppl from going through the algorithm again
    # No need to ensure that the infected ppl quarantine for more than 14 days even if they recover earlier
    # because they should be unable to transmit the pandemic either way
    if len(edges) > 0: # Some nodes have zero edges going into it (cut because of other ppl)
        # Contact tracing and removing the edges of the close contacts
        close_contacts = []
        for neighbor in G.neighbors(node):
            if G[node][neighbor]['weightedbeta'] > contact_thre:
                close_contacts.append(neighbor)
        # print(node,close_contacts)
        for c_contact in close_contacts:
            if G.nodes[c_contact]["isolate"] == False:
                # print('Starting the isolation:',c_contact)
                isolate(G,c_contact, remain_contact) # Isolate themselves for at least 14 day
                G.nodes[c_contact]['iso14day'] = 1
        # Testing and quarantine themselves after showing symptoms
        G.nodes[node]['removededges'] = copy.deepcopy(edges) # the connections to be removed
        for i in edges:
            G.edges[i[:2]]['weightedbeta'] = 0 # prevent secondary infection from these edges


# def vaccination(G,node):
    


def restore_edges(G, node):
    G.nodes[node]["isolate"] = False
    removededges = G.nodes[node]['removededges']
    num_removed = len(removededges)
    if num_removed != 0:
        for removed_edge in removededges:
            G.edges[removed_edge[:2]]['weightedbeta'] = removed_edge[2]['weightedbeta']



def transmission(G, N, if_isolate, test_ct_quar, beta_mu, contact_thre, remain_contact):
    # returns a network after one iteration of an input network
    
    secinf = [] # Empty list to collect secondary infection created for each Infectious node
    
    for i in range(N):
        G.nodes[i]['new_state'] = False
    
    # Start each loop with no ppl of new state - used to exempt ppl from a second change of state
    # We exempt the second change in states because it is biased depending on the node order, also the R number and other rates will be biased too
    # e.g. The infectious nodes can infect both nodes before (A) and after (B) them in the order as long as they are connected, but only B have the chance to be looped again.
    #      This will make B's Iday += 1, which is unfair for A.
    
    # Loop through ppl of all the states and take logic changes
    for node in G:
        
        if test_ct_quar == True:
            if G.nodes[node]['iso14day'] >= 1: # Ppl who are in the quarantine has this value >= 1
                if G.nodes[node]['iso14day'] == 15: # Restoring the edges when they reach the 15th day of quarantine
                    # print('Ending the isolation:',node)
                    restore_edges(G, node)
                    G.nodes[node]['iso14day'] = 0
                else:
                    G.nodes[node]['iso14day'] += 1
        else: pass
        
        # Exposed ppl (infected but not infectious)
        if G.nodes[node]['state'] == "E" and G.nodes[node]['new_state'] == False:
            G.nodes[node]['Iday'] += 1                      

            # Node becomes infectious after the latent period (different for each mode)
            if G.nodes[node]['Iday'] >= G.nodes[node]['dayLT']:
                G.nodes[node]['state'] = "I"
                G.nodes[node]['new_state'] = True
                
                # Secondary infection starts on the same day exposed ppl become infectious
                secondary_inf(G, node)
        
        # Infectious ppl
        elif G.nodes[node]['state'] == "I":
            G.nodes[node]['Iday'] += 1
            
            # Infected ppl either recover or die when they have reached their recovery/death day - allows the secondary infection to be reported
            # Recovering infectious ppl
            if G.nodes[node]['Iday'] >= G.nodes[node]['dayRC'] and G.nodes[node]['fate'] == True: # fate need to be examined, because dayRC of dying ppl is 0
                # print('recovered:',node)
                G.nodes[node]["state"] = 'R'
                secinf.append(G.nodes[node]['secinf']) # secondary infection reported
                if test_ct_quar == True:
                    # print('Ending the quarantine:',node)
                    G.nodes[node]["quarantine"] = False # Tested positive to be forced quarantine
                # Reestablish the weights of the lost edges
                if if_isolate == True:
                    restore_edges(G, node)
            
            # Dying infectious ppl
            elif G.nodes[node]['Iday'] >= G.nodes[node]['dayD'] and G.nodes[node]['fate'] == False: # fate need to be examined, because dayD of recovering ppl is 0
                # print('dead:', node)
                G.nodes[node]["state"] = 'D'
                secinf.append(G.nodes[node]['secinf']) # secondary infection reported
            
            # Other ppl stays infectious
            else:
                # Infectious ppl may take tests and be forced to completely quarantine. Their close contacts will also quarantine themselves for 14 days.
                if test_ct_quar == True and G.nodes[node]["quarantine"] == False and G.nodes[node]['Iday'] >= G.nodes[node]['dayIB']:
                    quarantine(G,node,remain_contact, contact_thre)
                            
                # Infectious ppl may isolate themselves after showing symptoms (after incubation period): 
                elif if_isolate == True and G.nodes[node]['Iday'] >= G.nodes[node]['dayIB']:
                    secondary_inf(G, node) # still capable of infecting the close contacts at home
                    G.nodes[node]['iso14day'] = 0 # avoid being out of self-isolation after 14 days if they have been contact traced before
                    if G.nodes[node]["isolate"] == False:
                        isolate(G,node, remain_contact)
                    
                else:
                    # print('Secondary infection by:', node)
                    # Secondary infection happens if nodes stay infectious and are not completely quarantined
                    for neighbor in G.neighbors(node):
                        if G.nodes[neighbor]['new_state'] == False: # exempt the new S ppl (turn from R) to be infected on the same day
                            x = random.random()
                            if G.nodes[neighbor]['state'] == 'S' and x < G[node][neighbor]['weightedbeta']:
                                G.nodes[neighbor]['state'] =  'E'
                                # if this node was newly infected on this day, then this node should not loop through 'E' logic since Iday = 0
                                G.nodes[neighbor]['new_state'] = True
                                G.nodes[node]['secinf'] += 1 # added on secondary infection count for each node
                                
                                # If latent period = 0 is possible - the neighbors are able to infect others on the first day being infectious
                                # Node becomes infectious after the latent period (different for each mode)
                                if G.nodes[neighbor]['Iday'] >= G.nodes[node]['dayLT']:
                                    G.nodes[neighbor]['state'] = "I"
                                    G.nodes[neighbor]['new_state'] = True
                                    
                                    # Secondary infection starts on the same day exposed ppl become infectious
                                    secondary_inf(G, neighbor)
    
        # Recovered ppl
        elif G.nodes[node]['state'] == "R":
            G.nodes[node]['Rday'] += 1

            # The recovered ppl will be able to be reinfected (turn susceptible) after the reinfection date 
            if G.nodes[node]['Rday'] >= G.nodes[node]['dayReinf']:
                G.nodes[node]['Iday'] = 0 # renewed for reinfection count
                G.nodes[node]['state'] = 'S'
                G.nodes[node]['new_state'] == True
                
    
    return G,secinf



def count(G): # counts the number of each group 
    S = 0; E = 0; I = 0; R = 0; D = 0
    for i in G.nodes:
        if G.nodes[i]["state"] == "S":
            S += 1
        elif G.nodes[i]["state"] == "E":
            E += 1
        elif G.nodes[i]["state"] == "I":
            I += 1
        elif G.nodes[i]["state"] == "R":
            R += 1
        elif G.nodes[i]["state"] == "D":
            D += 1
    return S,E,I,R,D
