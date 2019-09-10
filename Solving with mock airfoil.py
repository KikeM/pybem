#!/usr/bin/env python
# coding: utf-8

# In[31]:


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from math import pi


# In[52]:


def cl(alpha):
    
    from math import pi
    
    return 2*pi*(alpha)

def cd(cl):
    
    return 0.012 + 0.05 * cl**2


# In[122]:


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

# Rotation
omega = 2700 # rpm
omega = omega * (2 * pi) / 60
omega

# Lenghts
D = 2.032 # meters

r_tip = D / 2
r_hub = r_tip * 0.2

(r_tip, r_hub)


# In[123]:


from scipy.optimize import newton_krylov


# In[124]:


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


# In[125]:


def tangentail_induction(phi):
    
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


# In[126]:


axial_induction(2)
tangentail_induction(2)


# In[127]:


def residual(phi, r):
    
    _phi = np.deg2rad(phi)
    
    _v = V_inf / omega / r
    
    ax = axial_induction(phi)
    tng = tangentail_induction(phi)
    
    return np.tan(_phi) - _v * ((1 + ax) / (1-tng)) 


# In[128]:


from functools import partial


# In[129]:


r_hub


# In[130]:


r_tip


# In[ ]:





# In[131]:


plt.plot([_residual(element) for element in np.linspace(-20,20)])


# In[132]:


sol = 50

solutions = []

for radius in np.linspace(r_hub, r_tip):

    print(radius, sol)
    _residual = partial(residual, r = radius)
    try:
        sol = newton_krylov(_residual, sol)
    except Exception as ex:
        print(ex)
    
    solutions.append([radius, sol.tolist()])


# In[133]:


pd.DataFrame(solutions, columns = ['r', 'phi']).set_index('r').plot()


# In[ ]:





# In[ ]:




