# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:26:49 2021

@author: XPS
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
col_list = ["cumCasesBySpecimenDate"]
df = pd.read_csv("data_2021-Nov-07.csv", usecols=col_list)


# 26th march is the time first lockdown comes into force (first data at 11st Feb.)
days = 18+25
N = 8.982e6 # London population
data = df[-(days):]


# SIR
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0, M0 = 1, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - M0
# Average weighted infection rate beta, and mean recovery rate gamma, (in 1/days), morality rate morality
beta, gamma, mortality = 0.24, 4./100 ,1./1000     # Note that this beta was scaled with the contact rate r


# SEIR
# Initial number of infected and recovered individuals, I0 and R0. Exposed (infected but not infectious) and dead individials, E0 and M0 
I02, E02, R02, M02 = 1, 0, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S02 = N - I02 - R02 - M02 - E02
# Average weighted infection rate beta, and mean recovery rate gamma, (in 1/days), morality rate morality, rate to turn into infectious epsilon
beta2, gamma2, mortality2, epsilon2 = 0.3, 4./100 ,1./1000, 1     # Note that this beta was scaled with the contact rate r


R_0 = beta/(gamma+mortality) # average weighted infection rate * average removing rate
print('R_01 = ',R_0)
R_02 = beta2/(gamma2+mortality2) # average weighted infection rate * average removing rate
print('R_02 = ',R_02)



# A grid of time points (in days)
t = np.linspace(0, days, days)

# The SIR model differential equations.
def deriv_SIR(y, t, N, beta, gamma, mortality):
    S, I, R, M = y
    dSdt = - beta * S * I / N
    dIdt = beta * S * I / N - (gamma+mortality) * I
    dRdt = gamma * I
    dMdt = mortality * I
    return dSdt, dIdt, dRdt, dMdt


# The SIR model differential equations.
def deriv_SEIR(y, t, N, beta, gamma, mortality, epsilon):
    S, E, I, R, M = y
    dSdt = - beta * S * I / N
    dEdt = beta * S * I / N - epsilon * E
    dIdt = epsilon * E - (gamma+mortality) * I
    dRdt = gamma * I
    dMdt = mortality * I
    return dSdt, dEdt, dIdt, dRdt, dMdt


# Initial conditions vector
y0 = S0, I0, R0, M0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv_SIR, y0, t, args=(N, beta, gamma, mortality))
S, I, R, M = ret.T

# Initial conditions vector
y02 = S02, E02, I02, R02, M02
# Integrate the SEIR equations over the time grid, t.
ret2 = odeint(deriv_SEIR, y02, t, args=(N, beta2, gamma2, mortality2, epsilon2))
S2, E2, I2, R2, M2 = ret2.T



# Plot the data for both models and the real data
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, np.flip(np.array(data)), 'black', alpha=0.5, lw=2, label='Real data (London before lockdown)')
ax.plot(t, I, 'r', alpha=0.5, lw=2, label='Cumulated cases from SIR')
ax.plot(t, I2, 'b', alpha=0.5, lw=2, label='Cumulated cases from SEIR')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number')
#ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
# plt.title("Infection rate=%r  Recovery rate=%r mortality rate=%r \n Population=%d"%(beta, gamma, mortality,N))
plt.yscale("log")
plt.title("Analytical Models & Real Data")
plt.savefig("plots of SIR, SEIR and real world", dpi = 300)
plt.show()

print('R_0 of SIR =',R_0)
print('R_0 of SEIR =',R_02)
#%% R_eff plot for both models - Only plot with higher days! Otherwise the change is too small

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, R_0/N*S, 'r', alpha=0.5, lw=2, label='Cumulated cases from SIR')
ax.plot(t, R_02/N*S2, 'b', alpha=0.5, lw=2, label='Cumulated cases from SEIR')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number')
#ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
# plt.title("Infection rate=%r  Recovery rate=%r mortality rate=%r \n Population=%d"%(beta, gamma, mortality,N))
plt.title("Effective R Number of the Analytical Models")
plt.legend()
plt.show()
