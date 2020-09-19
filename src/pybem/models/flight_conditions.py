from math import pi

import numpy as np
import pandas as pd
from fluids.atmosphere import ATMOSPHERE_1976


class FlightConditions:
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

        self.atmosphere = ATMOSPHERE_1976(Z=altitude)

        self.v = airspeed * (1000.0 / 3600.0)
        self.omega = omega * (2 * pi) / 60.0

    @property
    def rho(self):
        return self.atmosphere.rho