# pybem


[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.org/KikeM/pybem.svg?branch=master)](https://travis-ci.org/KikeM/pybem)

Blade Element Method implementation for propeller calculations.

## Quickstart

```python
import numpy as np


from pybem import BladeElementMethod

bem = BladeElementMethod()
    

# Quick polar
# -----------------------------
def cl(alpha):

    from math import pi
    
    return 2*pi*(alpha)

def cd(cl):
    
    return 0.012 + 0.05 * cl**2

alpha = np.linspace(-20, 20, 50)
alpha_r = np.deg2rad(alpha)

cl_alpha = cl(alpha_r)
cd_polar = cd(cl_alpha)

# -------------------------------

# Airspeed
V_inf = 300 # km/h

# Rotation
omega = 2700 # rpm

# Lenghts
D = 2.032 # meters
r_tip = D / 2
r_hub = r_tip * 0.2

# Propeller
# ---------------------------------
beta = np.linspace(60, 40)
r    = np.linspace(r_hub, r_tip)
# ---------------------------------

# Load data
bem.load_airfoil(alpha, cl_alpha, cd_polar)
bem.load_flight_conditions(V_inf, omega, altitude = 10000)
bem.load_propeller(r_hub, r_tip, r, beta)

# Test!
bem.compute_inflow_angle(r_tip)
```
