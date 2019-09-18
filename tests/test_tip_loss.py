from math import pi

import numpy as np
import pandas as pd

from pybem import BladeElementMethod


def get_polar():
    # Read and parse polar
    polar_df = pd.read_csv(r'pybem\data\naca2415.csv')
    polar_df = polar_df.set_index('Alpha')[['Cl', 'Cd']]

    # Extract numpy arrays
    alpha    = polar_df.index.values
    cl_alpha = polar_df['Cl'].values
    cd_polar = polar_df['Cd'].values
    
    return alpha, cl_alpha, cd_polar

bem = BladeElementMethod()
    
# Polar
alpha, cl_alpha, cd_polar = get_polar()

# Lenghts
D = 0.152 * 2 # meters
r_tip = D / 2.0
r_hub = r_tip * 0.15 

# Propeller
# ---------------------------------
r    = np.linspace(r_hub, r_tip)
_r = r / r_tip
beta = np.linspace(50, 20, 50)
chord = (1.1 - _r**4.0  * np.exp(_r - 1)) * r
# ---------------------------------

# Load data
bem.load_airfoil(alpha, cl_alpha, cd_polar)

J = 0.1
bem.load_similarity(J = J)

# bem.load_flight_conditions(V_inf, omega, altitude = 10000)
prop = bem.load_propeller(dist_r=r, dist_beta=beta, dist_chord = chord, n_blades = 4)
bem.set_tip_loss(True)

loads = bem.compute_loads(dr = 0.01)

assert(loads['C_T'][-1]>0)
