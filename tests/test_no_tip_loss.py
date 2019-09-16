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

alpha = np.linspace(-40, 40, 50)
alpha_r = np.deg2rad(alpha)

cl_alpha = cl(alpha_r)
cd_polar = cd(cl_alpha)

# -------------------------------

# Airspeed
V_inf = 40/1000 # km/h

# Rotation
omega = 50 # rpm

# Lenghts
D = 0.237 # meters
r_tip = D / 2
r_hub = r_tip * 0.15

# Propeller
# ---------------------------------
beta = np.linspace(45, 20)
r    = np.linspace(r_hub, r_tip)
# ---------------------------------

# Load data
bem.load_airfoil(alpha, cl_alpha, cd_polar)
bem.load_flight_conditions(V_inf, omega, altitude = 10000)
bem.load_propeller(r_hub, r_tip, r, beta, 4)
bem.set_tip_loss(False)

# Test!
phi = bem.compute_inflow_angle(r_tip)

loads = bem.compute_loads(dr = 0.0001)

assert(loads['T_hat'][-1]>0)