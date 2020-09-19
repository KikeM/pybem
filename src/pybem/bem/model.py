from math import pi

import numpy as np
from scipy.optimize import newton_krylov

from pybem import Airfoil, FlightConditions, Propeller


class BladeElementMethod:
    """Blade Element Method implementation.
    This class solves the BEM equations to compute the performance of 
    a propeller. 

    Parameters
    ----------
    flight
    airfoil
    propeller
    J
    C_T
    C_Q
    """
    def __init__(self):

        # Flight conditions
        self.flight = None

        # Airfoil
        self.airfoil = None

        # Propeller
        self.propeller = None

        self.F = 1

        self.beta = None
        self.r = None
        self.dr = None
        self.D = None

        self.J = None

        self.C_T = None
        self.C_Q = None

        self.include_tip_loss = None

    def load_airfoil(self, alpha, polar_cl, polar_cd):
        """Load airfoil polar.

        Parameters
        ----------
        alpha: np.array
            Angle of attack in degrees
        polar_cl: np.array
        polar_cd: np.array
        """

        self.airfoil = Airfoil(alpha, polar_cl, polar_cd)

    def load_similarity(self, J):
        """Load advance ratio similarity coefficient.

        Parameters
        ----------
        J: float
        """
        self.J = J

    def load_flight_conditions(self, airspeed, omega, altitude):

        self.flight = FlightConditions(airspeed, omega, altitude)

    def load_propeller(self, dist_r, dist_beta, dist_chord, n_blades=2):
        """Create nondimensional propeller.

        Parameters
        ----------
        dist_r: np.array
        dist_beta: np.array
        dist_chord: np.array
        n_blades: int

        Returns
        -------
        Propeller object
        """
        # Get geometry
        r_hub = dist_r[0]
        r_tip = dist_r[-1]

        _D = 2.0 * r_tip

        self.D = _D

        # Non-dimensional radius distribution
        _dist_r = dist_r / _D

        # Create propeller
        self.propeller = Propeller(r_hub, r_tip, _dist_r, dist_beta, dist_chord)

        # Load number of blades
        self.B = n_blades

        return self.propeller

    def set_tip_loss(self, flag=True):
        """Activate Prandtl tip loss model.

        Parameters
        ----------
        flag: bool
        """
        self.include_tip_loss = flag

    def compute_induction_axial(self, r, phi, beta):

        _F = self.compute_tip_loss(r, phi)
        _solidity = self.solidity(r)

        # Convert to radians
        _phi, _beta = np.deg2rad([phi, beta])

        # Compute angle of attack
        alpha = beta - phi

        # Polar
        _cl = self.airfoil.cl(alpha)
        _cd = self.airfoil.cd(_cl)

        NUM = 4.0 * _F * (np.sin(_phi)) ** 2.0
        DEN = _solidity * (_cl * np.cos(_phi) - _cd * np.sin(_phi))

        frac = NUM / DEN - 1

        return 1.0 / frac

    def compute_induction_tangential(self, r, phi, beta):

        # Convert to radians
        _beta, _phi = np.deg2rad([beta, phi])

        # Compute angle of attack
        alpha = beta - phi

        F = self.compute_tip_loss(r, phi)
        solidity = self.solidity(r)

        # Polar
        _cl = self.airfoil.cl(alpha)
        _cd = self.airfoil.cd(_cl)

        NUM = 4.0 * F * np.sin(_phi) * np.cos(_phi)
        DEN = solidity * (_cl * np.sin(_phi) - _cd * np.cos(_phi))

        frac = NUM / DEN + 1

        return 1.0 / frac

    def solidity(self, r):
        """Local element solidity.

        Parameters
        ----------
        r: float

        Returns
        -------
        float
        """
        _D = self.D
        _B = self.B

        _chord = self.propeller.chord(r)

        NUM = _B * _chord
        DEN = 2.0 * pi * (r * _D)

        return NUM / DEN

    def compute_inflow_angle(self, r, x0=50):
        """Solve nonlinear inflow angle equation.

        Parameters
        ----------
        r: float
            r / D parameter

        x0: float
            Initial condition

        Returns
        -------
        float
            If NO convergence, returns np.nan
        """
        self.beta = self.propeller.beta(r)
        self.r = r

        try:
            sol = newton_krylov(self._residual, x0)
        except Exception as ex:
            print(ex)
            sol = np.array(np.nan)

        return sol.item()

    def compute_loads(self, dr=0.01):
        """Compute blade loads.

        Parameters
        ----------
        dr: float
            Spacing between stations.

        Returns
        -------
        tuple: thrust, torque

        """
        _J = self.J
        _D = self.D

        _r_hub = self.propeller.r_hub / _D
        _r_tip = self.propeller.r_tip / _D

        # Number of steps
        N = np.floor((_r_tip - _r_hub) / dr)

        # Create nondimensional radius distribution
        r_space = np.linspace(start=_r_hub, stop=_r_tip, num=N)

        # Initial condition
        C_T = [0.0]
        C_Q = [0.0]

        # Initial condition is 20% larger than the twist angle
        phi0 = np.arctan(_J / _r_hub)
        phi0 = np.rad2deg(phi0)

        phi_space = [phi0]
        F_space = [1.0]

        idx = 0  # Index to control numpy arrays
        for r in r_space[:-1]:

            # Compute induction angle
            phi = self.compute_inflow_angle(r, phi_space[idx])

            # Compute induction coefficients
            axi = self.compute_induction_axial(r, phi, self.beta)
            tng = self.compute_induction_tangential(r, phi, self.beta)

            # Tip loss
            _F = self.compute_tip_loss(r, phi)

            # Compute Euler **implicit** derivative
            F_T = 4.0 * pi * _J ** 2.0 * (r + dr) ** 1.0 * (1.0 + axi) * axi * _F
            F_Q = 4.0 * pi * _J ** 1.0 * (r + dr) ** 3.0 * (1.0 + axi) * tng * _F

            # Integrate
            C_T.append(C_T[idx] + dr * F_T)
            C_Q.append(C_Q[idx] + dr * F_Q)

            # Save state
            F_space.append(_F)
            phi_space.append(phi)

            idx += 1

        C_T = np.array(C_T)
        C_Q = np.array(C_Q)

        self.C_T = C_T
        self.C_Q = C_Q

        self.dr = dr
        self.N = N

        # Pack up results
        result = dict()

        result["r"] = r_space
        result["C_T"] = C_T
        result["C_Q"] = C_Q
        result["F"] = F_space
        result["phi"] = phi_space

        return result

    def _residual(self, phi):
        """Nonlinear inflow angle equation residual.

        Parameters
        ----------
        phi: float
            Inflow angle, in deg

        Returns
        -------
        float
            Equation residual
        """
        _phi = np.deg2rad(phi)

        _J = self.J
        _r = self.r

        # Compute incidence coefficients
        ind_ax = self.compute_induction_axial(_r, phi, self.beta)
        ind_tan = self.compute_induction_tangential(_r, phi, self.beta)

        # Compute residual
        return np.tan(_phi) - (_J / _r) * ((1.0 + ind_ax) / (1.0 - ind_tan))

    def compute_tip_loss(self, r, phi):
        """Prandtl tip loss coefficient.

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
            raise ValueError("You need to invoke the set_tip_loss method first!")

        if self.include_tip_loss == False:

            return 1.0

        else:

            _phi = np.deg2rad(phi)

            # Recover geometry
            _B = self.B
            _D = self.D
            _R = self.propeller.r_tip

            _r = r * _D

            NUM = _B * (_R - _r)
            DEN = 2.0 * _r * np.sin(_phi)

            f_tip = NUM / DEN

            return 2.0 * np.arccos(np.exp(-f_tip)) / pi

    def compute_hub_loss(self, r, phi):
        """Hub loss coefficient.

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
