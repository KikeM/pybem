from math import pi

import numpy as np
from scipy.optimize import newton_krylov

from pybem import Airfoil, FlightConditions, Propeller


class BladeElementMethod():

    def __init__(self):

        # Flight conditions
        self.flight = None

        # Airfoil
        self.airfoil = None

        # Propeller
        self.propeller = None

        self.solidity = 0.5
        self.F        = 1

        self.beta = None
        self.r    = None
        self.dr   = None

        self.T = None
        self.Q = None

        self.T_hat = None
        self.Q_hat = None

        self.q_inf = None

        self.include_tip_loss = None

    def load_airfoil(self, alpha, polar_cl, polar_cd):

        self.airfoil = Airfoil(alpha, polar_cl, polar_cd)

    def load_flight_conditions(self, airspeed, omega, altitude):

        self.flight = FlightConditions(airspeed, omega, altitude)

    def load_propeller(self, r_hub:float, r_tip:float, r_loc, beta_dist, B = 1):

        self.propeller = Propeller(r_hub, r_tip, r_loc, beta_dist)
        
        # Load number of blades
        self.B = B

    def set_tip_loss(self, flag = True):
        self.include_tip_loss = flag

    def induction_axial(self, r, phi, beta):

        F        = self.compute_tip_loss(r, phi)
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

    def induction_tangential(self, r, phi, beta):
        
        # Convert to radians
        _beta, _phi = np.deg2rad([beta, phi])
        
        # Compute angle of attack
        alpha = beta - phi

        F        = self.compute_tip_loss(r, phi)
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

    def compute_loads(self, dr = 0.05):
        """
        Compute blade loads.

        Parameters
        ----------
        dr: float
            Spacing between stations.

        Returns
        -------
        tuple: thrust, torque

        """
        _r_hub = self.propeller.r_hub
        _r_tip = self.propeller.r_tip
        
        N  = np.floor((_r_tip - _r_hub) / dr)
        
        omega = self.flight.omega

        print(f"Using {N} stations")

        T_hat = [0.0]
        Q_hat = [0.0]

        # TODO: create property
        q_inf = 4 * pi * self.flight.atmosphere.rho * (self.flight.v)**2.0

        r_space = np.linspace(_r_hub, _r_tip, N)
        phi_space = []
        F_space = []


        idx = 0
        for r in r_space[0:]:

            # Compute induction angle
            phi = self.compute_inflow_angle(r)
            
            # Compute induction coefficients
            axi = self.induction_axial(r, phi, self.beta)
            tng = self.induction_tangential(r, phi, self.beta)
            
            # Tip loss
            F = self.compute_tip_loss(r, phi)
            
            # Compute forcing terms
            F_T = (r+dr)**1.0 * (1 + axi) * axi * F
            F_Q = (r+dr)**3.0 * (1 + axi) * tng * F

            T_hat.append(T_hat[idx] + dr * F_T)
            Q_hat.append(Q_hat[idx] + dr * F_Q)
            
            # Save state
            F_space.append(F)
            phi_space.append(phi)
            
            idx +=1

        T_hat = np.array(T_hat[1:])
        Q_hat = np.array(Q_hat[1:])
            
        # Give proper dimensions
        T = np.array(T_hat)
        Q = np.array(Q_hat)

        T *= q_inf
        Q *= q_inf * omega

        self.T = T
        self.Q = Q

        self.T_hat = T_hat
        self.Q_hat = Q_hat

        self.dr = dr
        self.N  = N
        self.q_inf = q_inf
        
        # Pack up results
        result = dict()

        result['r']     = r_space
        result['T']     = T
        result['Q']     = Q
        result['T_hat'] = T_hat
        result['Q_hat'] = Q_hat
        result['F']     = F_space
        result['phi']   = phi_space

        return result

    def _residual(self, phi):

        _phi = np.deg2rad(phi)

        ind_ax  = self.induction_axial(self.r, phi, self.beta) 
        ind_tan = self.induction_tangential(self.r, phi, self.beta)

        _v = self.flight.v / self.flight.omega / self.r

        return np.tan(_phi) - _v * ((1.0 + ind_ax) / (1.0 - ind_tan))

    def compute_tip_loss(self, r, phi):
        """
        Prandtl tip loss coefficient.

        Parameters
        ----------
        r: float

        phi: float
            In degrees.

        Returns
        -------
        F: float
        """
        # Check correct use
        if self.include_tip_loss is None:
            raise ValueError('You need to invoke the set_tip_loss method first!')

        if self.include_tip_loss == False:
            return 1.0
        else:
            _phi = np.deg2rad(phi)

            _B = self.B
            _R = self.propeller.r_tip
            
            N = _B * (_R - r)
            D = 2 * r * np.sin(_phi)

            f_tip = N / D

            return 2.0 * np.arccos(np.exp(-f_tip)) / pi

    def compute_hub_loss(self, r, phi):
        """
        Hub loss coefficient.

        Parameters
        ----------
        r: float

        phi: float
            In degrees.

        Returns
        -------
        F: float
        """
        _phi = np.deg2rad(phi)

        _B = self.B
        _R = self.propeller.r_hub
        
        N = _B * (r - _R)
        D = 2 * r * np.sin(_phi)

        f_hub = N / D

        return 2.0 * np.arccos(np.exp(-f_hub)) / pi
