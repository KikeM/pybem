from math import pi

import numpy as np
import pandas as pd
from fluids.atmosphere import ATMOSPHERE_1976
from scipy.interpolate import interp1d


class FlightConditions():

    def __init__(self, airspeed, omega, altitude):
        """
        Parameters
        ----------
        V: float
            airspeed in km/h
        omega: float
            Blade angular speed
        altitude: float
            Height from the ground, in meters. 
        """

        self.atmosphere = ATMOSPHERE_1976(Z = altitude)

        self.v     = airspeed * (1000.0 / 3600.0)
        self.omega = omega * (2 * pi) / 60.0

    @property
    def rho(self):
        return self.atmosphere.rho

class Airfoil():
    """
    Airfoil aerodynamic properties.
    """

    def __init__(self, alpha, polar_cl, polar_cd):

        # Store the data
        self.polar_cl = polar_cl
        self.polar_cd = polar_cd

        # Create interpolants
        self.interpolant_cl = interp1d(alpha,    polar_cl)
        self.interpolant_cd = interp1d(polar_cl, polar_cd)

    def cl(self, alpha):
        """
        Lift coefficient.

        Parameters
        ----------
        alpha: float

        Returns
        -------
        cl: float
            Interpolation based on lift polar.
        """
        self.alpha = alpha
        _cl        = self.interpolant_cl(alpha).item()
        self._cl   = _cl

        return _cl

    def cd(self, cl):
        """
        Drag coefficient.

        Parameters
        ----------
        cl: float

        Returns
        -------
        cd: float
            Interpolation based on drag polar.
        """
        self._cl = cl
        _cd      = self.interpolant_cd(cl).item()
        self._cd = _cd

        return _cd


class Propeller():

    def __init__(self, r_hub:float, r_tip:float, r_dist, beta_dist, chord_dist):

        # Get blade bounds
        self.r_hub = r_hub
        self.r_tip = r_tip

        # Get airfoil properties
        self.dist_twist = beta_dist
        self.dist_chord = chord_dist
        self.dist_r     = r_dist

        self.interpolant_twist = interp1d(r_dist, beta_dist)
        self.interpolant_chord = interp1d(r_dist, chord_dist)

        self.r      = None
        self._beta  = None
        self._chord = None

    def beta(self, r):
        """
        Airfoil twist at location r.

        Parameters
        ----------
        r: float
            Location based on r_loc.

        Returns
        -------
        beta: float
        """
        self.r     = r
        
        _beta      = self.interpolant_twist(r).item()
        
        self._beta = _beta

        return _beta

    def chord(self, r):
        """
        Airfoil chord at location r.

        Parameters
        ----------
        r: float
            Location based on r_loc.

        Returns
        -------
        chord: float
        """
        self.r      = r
        
        _chord      = self.interpolant_chord(r).item()
        
        self._chord = _chord

        return _chord
