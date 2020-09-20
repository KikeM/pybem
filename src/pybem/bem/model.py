from functools import partial
from math import pi

import numpy as np
from scipy.optimize import newton_krylov

from .loss import compute_hub_loss, compute_tip_loss


class BladeElementMethod:
    """Blade Element Method implementation.

    This class solves the BEM equations to compute the performance and
    load distribution along a propeller.

    Parameters
    ----------
    J : float
        Advance ratio.
    propeller
    flight
    tip_loss : bool
    hub_loss : bool

    Attributes
    ----------
    J : float
        Advance ratio.
    propeller
    flight
    tip_loss : bool
    hub_loss : bool
    C_T : np.array
    C_Q : np.array
    """

    N_SECTIONS = 100

    def __init__(self, J, propeller, flight, tip_loss=False, hub_loss=False):

        self.J = J

        self.flight = flight
        self.propeller = propeller

        self.tip_loss = tip_loss
        self.hub_loss = hub_loss

        # Dimensionless thrust and torque distribution
        self.C_T = None
        self.C_Q = None

    def solve(self):

        r_min = self.propeller.radii[0]

        # Create dimensionless radius distribution
        r_dist = np.linspace(start=r_min, stop=1.0, num=self.N_SECTIONS)

        phi0 = self.propeller.compute_beta(r_min)
        for r in r_dist:

            # Solve inflow angle.
            phi = self.compute_inflow_angle(r=r, phi0=phi0)

            phi0 = phi

    def _residual(self, phi, r):
        """Nonlinear inflow angle equation residual.

        Parameters
        ----------
        phi : float
            Incidence angle, in degrees.
        r : float
            Dimensionless radius, r / R.
        beta : float
            Twist angle.

        Returns
        -------
        residual : float
        """

        propeller = self.propeller

        # Compute angle of attack
        beta = propeller.compute_beta(r)
        alpha = beta - phi

        # Compute lift and drag coefficients
        cl = 1.0
        cd = 0.0

        # Compute incidence coefficients
        sigma = propeller.compute_solidity(r)

        a = self.compute_axial_coefficient(r=r, phi=phi, cl=cl, cd=cd, sigma=sigma, F=F)
        b = self.compute_angular_coefficient(
            r=r, phi=phi, cl=cl, cd=cd, sigma=sigma, F=F
        )

        # Compute residual
        J = self.J

        SUM_1 = np.sin(phi) / (1.0 + a)
        SUM_2 = np.cos(phi) / (J * r * (1.0 - b))

        res = SUM_1 - SUM_2

        return res

    def compute_loss(self, r, phi):

        propeller = self.propeller

        F = 1.0  # Assume no loss

        if self.tip_loss == True:

            F_tip = compute_tip_loss(B=propeller.B, r=r, phi=phi)

            F *= F_tip

        if self.hub_loss == True:

            r_hubR = propeller.R_hub / propeller.R
            F_hub = compute_hub_loss(B=propeller.B, r=r, phi=phi, r_hubR=r_hubR)

            F *= F_hub

        return F

    @staticmethod
    def compute_ct(cl, cd, phi):
        """Compute section tangential force coefficient.

        Parameters
        ----------
        cl : float
        cd : float
        phi : float
            Incidence angle in radians.

        Returns
        -------
        ct : float
        """

        ct = cl * np.cos(phi) - cd * np.sin(phi)

        return ct

    @staticmethod
    def compute_cn(cl, cd, phi):
        """Compute section angular force coefficient.

        Parameters
        ----------
        cl : float
        cd : float
        phi : float
            Incidence angle in radians.

        Returns
        -------
        cn : float
        """

        cn = cl * np.sin(phi) + cd * np.cos(phi)

        return cn

    def compute_axial_coefficient(self, r, phi, cl, cd, sigma, F):

        ct = self.compute_ct(cl=cl, cd=cd, phi=phi)

        # Compute loss coefficients (if necessary)
        F = self.compute_loss(r=r, phi=phi)

        NUM = 4.0 * F * r * (np.sin(phi)) ** 2.0
        DEN = sigma * ct

        frac = NUM / DEN - 1.0

        a = 1.0 / frac

        return a

    def compute_angular_coefficient(self, r, phi, cl, cd, sigma, F):

        cn = self.compute_cn(cl=cl, cd=cd, phi=phi)

        # Compute loss coefficients (if necessary)
        F = self.compute_loss(r=r, phi=phi)

        NUM = 4.0 * F * r * np.sin(phi) * np.cos(phi)
        DEN = sigma * cn

        frac = NUM / DEN + 1.0

        b = 1.0 / frac

        return b

    def compute_inflow_angle(self, r, phi0=0.0):
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

        # Fix section
        func = partial(self._residual, r=r)

        try:
            phi = newton_krylov(func, phi0)
        except Exception as ex:
            print(ex)
            phi = np.array(np.nan)

        return phi.item()

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
            axi = self.compute_axial_coefficient(r, phi, self.beta)
            tng = self.compute_cn(r, phi, self.beta)

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
