#!/usr/bin/env python
# coding: utf-8

# http://hartzellprop.com/wp-content/uploads/Hartzell-Propeller-Others-Full-Catalog.pdf
# 
# https://www.airforce-technology.com/projects/aermacchisf260traine/
# 
# https://m-selig.ae.illinois.edu/props/propDB.html
# 
# https://www.airliners.net/aircraft-data/aermacchi-f-260/3
# 

# In[1]:


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from math import pi


# In[2]:


def cl(alpha):
    
    from math import pi
    
    return 2*pi*(alpha)

def cd(cl):
    
    return 0.012 + 0.05 * cl**2


# In[101]:


# Test polar
cl(np.deg2rad(10))
cd(1.0)

# Prandtl tip loss
F = 1

solidity = 0.2

# Foil angle
beta = 20

# Airspeed
V_inf = 300 # km/h
V_inf = V_inf * 1000/3600
V_inf

rho = 1.0

# Rotation
omega = 2700 # rpm
omega = omega * (2 * pi) / 60
omega

# Lenghts
D = 2.032 # meters

r_tip = D / 2
r_hub = r_tip * 0.2

(r_tip, r_hub)


# In[4]:


from scipy.optimize import newton_krylov


# In[5]:


def axial_induction(phi):
    
    # Convert to radians
    _beta, _phi = np.deg2rad([beta, phi])
    
    # Compute angle of attack
    _alpha = _beta - _phi
    
    # Polar
    _cl = cl(_alpha)
    _cd = cd(_cl)
    
    num = 4.0 * F * (np.sin(_phi))**2.0
    den = solidity * (_cl * np.cos(_phi) - _cd * np.sin(_phi))
    
    frac = num / den - 1
    
    return 1.0 / frac


# In[93]:


def tangential_induction(phi):
    
    # Convert to radians
    _beta, _phi = np.deg2rad([beta, phi])
    
    # Compute angle of attack
    _alpha = _beta - _phi
    
    # Polar
    _cl = cl(_alpha)
    _cd = cd(_cl)
    
    num = 4.0 * F * np.sin(_phi) * np.cos(_phi)
    den = solidity * (_cl * np.sin(_phi) - _cd * np.cos(_phi))
    
    frac = num / den + 1
    
    return 1.0 / frac


# In[7]:


axial_induction(2)
tangentail_induction(2)


# In[94]:


def residual(phi, r):
    
    _phi = np.deg2rad(phi)
    
    _v = V_inf / omega / r
    
    ax = axial_induction(phi)
    tng = tangential_induction(phi)
    
    return np.tan(_phi) - _v * ((1 + ax) / (1-tng)) 


# In[14]:


from functools import partial


# In[148]:


solutions = []

beta = 40  
sol = 50

for radius in np.linspace(r_hub, r_tip):

    _residual = partial(residual, r = radius)
    try:
        sol = newton_krylov(_residual, sol)
    except Exception as ex:
        solutions.append([radius, np.nan, beta])

    solutions.append([radius, sol.tolist(), beta])


# In[149]:


results_df = pd.DataFrame(solutions, columns = ['r', 'phi', 'beta']).set_index('r')


# In[150]:


ax = results_df.pivot(columns = 'beta', values = 'phi').plot()

ax.set_title('Incidence angle');
ax.set_ylabel('$\\phi$');
ax.set_xlabel('$r$');
ax.grid(True);


# # Integration

# In[151]:


def compute_inductions(x0, radius):
    
    # Create residual function with current position
    _residual = partial(residual, r = radius)
    
    # Attempt solution
    try:
        sol = newton_krylov(_residual, x0)
    except Exception as ex:
        sol = np.nan
        
    return sol
    


# In[152]:


(r_tip - r_hub) / 0.05


# In[153]:


dr = 0.01
N  = np.floor((r_tip - r_hub) / dr)

T_hat = [0.0]
Q_hat = [0.0]

phi = 50

beta = 40

q_inf = 4 * pi * rho * V_inf**2.0

r_space = np.linspace(r_hub, r_tip, N)

idx = 0
for r in r_space[1:]:
    
    # Compute induction angle
    phi = compute_inductions(phi, r)
    
    # Compute induction coefficients
    axi = axial_induction(phi)
    tng = tangential_induction(phi)
    
    # Compute forcing terms
    F_T = (r+dr)**1.0 * (1 + axi) * axi * F
    F_Q = (r+dr)**3.0 * (1 + axi) * tng * F
    
    
    T_hat.append(T_hat[idx] + dr * F_T)
    Q_hat.append(Q_hat[idx] + dr * F_Q)
    
    idx +=1

T_hat = np.array(T_hat)
Q_hat = np.array(Q_hat)
    
# Give proper dimensions
T = np.array(T_hat)
Q = np.array(Q_hat)

T *= q_inf
Q *= q_inf * omega


# In[154]:


fig, axes = plt.subplots(figsize=(8,5))

axes.plot(r_space, T_hat, label = '$\hat{T}$');
axes.plot(r_space, Q_hat, label = '$\hat{Q}$');

axes.set_title('Dimensionless thrust vs. torque')
axes.legend();
axes.grid(True);


# In[155]:


fig, axes = plt.subplots(ncols = 2, figsize=(15,5))

axes[0].plot(r_space, T);
axes[1].plot(r_space, Q);

axes[0].set_title('Thrust');
axes[1].set_title('Torque');

for ax in axes: ax.set_xlabel('$r$');
for ax in axes: ax.grid(True);
