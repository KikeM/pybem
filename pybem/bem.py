
import numpy as np
from scipy.optimize import newton_krylov
from pybem import Airfoil, Propeller, FlightConditions


class BladeElementMethod():

    def __init__(self):

        # Flight conditions
        self.flight = None

        # Airfoil
        self.airfoil = None

        # Propeller
        self.propeller = None

        self.solidity = 0.5
        self.F = 1

        self.beta = None
        self.r = None

    def load_airfoil(self, alpha, polar_cl, polar_cd):

        self.airfoil = Airfoil(alpha, polar_cl, polar_cd)

    def load_flight_conditions(self, airspeed, omega, altitude):

        self.flight = FlightConditions(airspeed, omega, altitude)

    def load_propeller(self, r_hub:float, r_tip:float, r_loc, beta_dist):

        self.propeller = Propeller(r_hub, r_tip, r_loc, beta_dist)

    def induction_axial(self, phi, beta):

        F        = self.F
        solidity = self.solidity

        # Convert to radians
        _phi, _beta = np.deg2rad([phi, beta])
        
        # Compute angle of attack
        alpha = beta - phi
        
        # Polar
        _cl = self.airfoil.cl(alpha)
        _cd = self.airfoil.cd(_cl)
        
        num = 4.0 * F * (np.sin(_phi))**2.0
        den = solidity * (_cl * np.cos(_phi) - _cd * np.sin(_phi))
        
        frac = num / den - 1
        
        return 1.0 / frac

    def induction_tangential(self, phi, beta):
        
        # Convert to radians
        _beta, _phi = np.deg2rad([beta, phi])
        
        # Compute angle of attack
        alpha = beta - phi

        F        = self.F
        solidity = self.solidity

        # Polar
        _cl = self.airfoil.cl(alpha)
        _cd = self.airfoil.cd(_cl)
        
        num = 4.0 * F * np.sin(_phi) * np.cos(_phi)
        den = solidity * (_cl * np.sin(_phi) - _cd * np.cos(_phi))
        
        frac = num / den + 1
        
        return 1.0 / frac

    def compute_inflow_angle(self, r, x0 = 50):

        self.beta = self.propeller.beta(r)
        self.r    = r

        try:
            sol = newton_krylov(self._residual, x0)
        except Exception as ex:
            print(ex)
            sol = np.array(np.nan)

        return sol.item()

    def _residual(self, phi):

        _phi = np.deg2rad(phi)

        ind_ax  = self.induction_axial(phi, self.beta) 
        ind_tan = self.induction_tangential(phi, self.beta)

        _v = self.flight.v / self.flight.omega / self.r

        return np.tan(_phi) - _v * ((1.0 + ind_ax) / (1.0 - ind_tan))

        
