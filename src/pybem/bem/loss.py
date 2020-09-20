"""Prandtl loss models to account for hub and tip vortices.
"""

from math import pi

import numpy as np


def func(phi, x, a):
    """Auxiliary function for the generalization of the tip losses.

    Parameters
    ----------
    phi : float
        Incidence angle in radians.
    x : float
        Dimensionless radius r / R
    a : float
        Constant to distinguish between tip and hub loss.

    Returns
    -------
    float
    """
    NUM = a - x
    DEN = 2.0 * np.sin(phi) * x

    return NUM / DEN


def compute_tip_loss(phi, r, B):
    """Compute Prandtl tip loss coefficient.

    Parameters
    ----------
    phi : float
        Incidence angle in radians.
    r : float
        Dimensionless radius, r / R
    B : int
        Number of blades

    Returns
    -------
    F_tip: float
    """
    f_tip = func(phi=phi, x=r, a=1.0)
    f_tip *= B

    F_tip = 2.0 * np.arccos(np.exp(-f_tip)) / pi

    return F_tip


def compute_hub_loss(phi, r, B, r_hubR):
    """Compute Prandtl hub loss coefficient.

    Parameters
    ----------
    phi : float
        Incidence angle in radians.
    r : float
        Dimensionless radius, r / R
    B : int
        Number of blades
    r_hubR : float
        Dimensionless location of the hub radius, R_hub / R

    Returns
    -------
    F_hub: float
    """
    f_hub = func(phi=phi, x=r, a=r_hubR)
    f_hub *= -B

    F_hub = 2.0 * np.arccos(np.exp(-f_hub)) / pi

    return F_hub