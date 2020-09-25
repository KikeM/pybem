from functools import partial
from math import pi

import numpy as np
from scipy.integrate import simps as integrate
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
    propeller : Propeller-like object
    flight : FlightConditions-like object
    tip_loss : bool
        Use the tip correction factor? Default is `False`.
    hub_loss : bool
        Use the tip correction factor? Default is `False`.
    r_dist : list of floats
    phi : list of floats
    C_T : list of floats
    C_Q : list of floats
    N_SECTIONS : int
    EPSILON : float
    """

    N_SECTIONS = 100
    EPSILON = 1e-6

    def __init__(self, J, propeller, flight=None, tip_loss=False, hub_loss=False):

        self.J = J

        self.flight = flight
        self.propeller = propeller

        self.tip_loss = tip_loss
        self.hub_loss = hub_loss

        # Dimensionless thrust and torque distribution
        self.r_dist = []
        self.phi_dist = []
        self.a_dist = []
        self.b_dist = []
        self.F_dist = []
        self.dCTdr_dist = []
        self.dCQdr_dist = []

        self.CT = None
        self.CQ = None

    @staticmethod
    def compute_ct(cl, cd, phi):
        """Compute section tangential force coefficient.

        Parameters
        ----------
        cl : float
        cd : float
        phi : float
            Incidence angle in degrees.

        Returns
        -------
        ct : float
        """

        phi = np.deg2rad(phi)

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
            Incidence angle in degrees.

        Returns
        -------
        cn : float
        """

        phi = np.deg2rad(phi)

        cn = cl * np.sin(phi) + cd * np.cos(phi)

        return cn

    @staticmethod
    def force_coeff_integrand(r, F, a):
        """Compute the slope dCT/dr to integrate CT.

        Parameters
        ----------
        r : float
        F : float
        a : float

        Returns
        -------
        dCTdr : float
        """

        dCTdr = 8.0 * pi * r * (1.0 + a) * a * F

        return dCTdr

    @staticmethod
    def torque_coeff_integrand(r, F, a, b, J):
        """Compute the slope dCT/dr to integrate CT.

        Parameters
        ----------
        r : float
        F : float
        a : float
        b : float
        J : float

        Returns
        -------
        dCQdr : float
        """

        dCQdr = 8.0 * pi * (r ** 3.0) * (1.0 + a) * b * F / J

        return dCQdr

    def solve(self):
        """Solve the Blade Element Problem

        Returns
        -------
        r_dist : array-like of floats
            Dimensionless radii.
        phi : array-like of floats
            Incidence angles for each r.
        """

        r_min = self.propeller.radii[0]

        # Create dimensionless radius distribution
        # Take into account when losses are actived we cannot solve at the exact section
        start = r_min if self.hub_loss == False else r_min * (1.0 + self.EPSILON)
        stop = 1.0 if self.tip_loss == False else 1.0 - self.EPSILON

        r_dist = np.linspace(start=start, stop=stop, num=self.N_SECTIONS)
        self.r_dist = r_dist

        phi0 = self.propeller.compute_beta(r_min)
        phi = []
        a = []
        b = []
        dCTdr = []
        dCQdr = []
        F = []
        for r in r_dist:

            # (BEM loop)
            # Solve inflow angle.
            _phi = self.compute_inflow_angle(r=r, phi0=phi0)

            # Save station value
            phi.append(_phi)

            # Update starting point for the next iteration
            phi0 = _phi

            # (Performance calculation loop)
            # Compute differential force and torque slopes
            _F = self.compute_prandtl_loss(r=r, phi=_phi)
            _a, _b = self.compute_induction_coefficients(r=r, phi=_phi)

            _dCTdr = self.force_coeff_integrand(r=r, a=_a, F=_F)
            _dCQdr = self.torque_coeff_integrand(r=r, a=_a, F=_F, b=_b, J=self.J)

            # Save station values
            F.append(_F)
            dCTdr.append(_dCTdr)
            dCQdr.append(_dCQdr)
            a.append(_a)
            b.append(_b)

        self.phi_dist = phi
        self.F_dist = F
        self.dCTdr_dist = dCTdr
        self.dCQdr_dist = dCQdr
        self.a_dist = a
        self.b_dist = b

        return r_dist, phi

    def integrate_forces(self):
        """Integrate thrust and torque distributions.

        Returns
        -------
        CT : float
        CQ : float

        Notes
        -----
        Uses `self` variables:
            - r_dist
            - dCTdr_dist
            - dCQdr_dist
        """

        # Creates integrator based on the radius distribution
        _integrate = partial(integrate, x=self.r_dist)

        CT = _integrate(self.dCTdr_dist)
        CQ = _integrate(self.dCQdr_dist)

        self.CT = CT
        self.CQ = CQ

        return CT, CQ

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

        # Compute angle of attack
        a, b = self.compute_induction_coefficients(r=r, phi=phi)

        # Compute residual
        J = self.J

        phi = np.deg2rad(phi)

        SUM_1 = np.sin(phi) / (1.0 + a)
        SUM_2 = J * np.cos(phi) / (r * (1.0 - b))

        res = SUM_1 - SUM_2

        return res

    def compute_induction_coefficients(self, r, phi):
        """Compute the velocity induction coefficients.

        Parameters
        ----------
        r : float
            Dimensionless coefficients.
        phi : float
            Incidence angle in degrees.

        Returns
        -------
        a : float
            Axial induction coefficient.
        b : float
            Tangential induction coefficient.
        """

        propeller = self.propeller

        # Compute angle of attack
        beta = propeller.compute_beta(r)
        alpha = beta - phi

        # Compute lift and drag coefficients
        cl = propeller.compute_cl(r=r, alpha=alpha)
        cd = propeller.compute_cd(r=r, alpha=alpha)

        # Compute incidence coefficients
        sigma = propeller.compute_solidity(r)

        a = self.compute_axial_coefficient(r=r, phi=phi, cl=cl, cd=cd, sigma=sigma)
        b = self.compute_angular_coefficient(r=r, phi=phi, cl=cl, cd=cd, sigma=sigma)

        return a, b

    def compute_prandtl_loss(self, r, phi):
        """Compute tip and hub losses according to the Prandtl model.

        Parameters
        ----------
        r : float
        phi : float
            Incidence angle in degrees

        Returns
        -------
        F : float
            Correction factor
        """

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

    def compute_axial_coefficient(self, r, phi, cl, cd, sigma):

        ct = self.compute_ct(cl=cl, cd=cd, phi=phi)

        # Compute loss coefficients (if necessary)
        F = self.compute_prandtl_loss(r=r, phi=phi)

        # Convert to radians
        phi = np.deg2rad(phi)

        NUM = 4.0 * F * r * (np.sin(phi)) ** 2.0
        DEN = sigma * ct

        # Prevent division by zero warning
        if np.isclose(DEN, 0.0):
            frac = np.inf
        else:
            frac = NUM / DEN - 1.0

        a = 1.0 / frac

        return a

    def compute_angular_coefficient(self, r, phi, cl, cd, sigma):

        cn = self.compute_cn(cl=cl, cd=cd, phi=phi)

        # Compute loss coefficients (if necessary)
        F = self.compute_prandtl_loss(r=r, phi=phi)

        # Convert to radians
        phi = np.deg2rad(phi)

        NUM = 4.0 * F * r * np.sin(phi) * np.cos(phi)
        DEN = sigma * cn

        # Prevent division by zero warning
        if np.isclose(DEN, 0.0):
            frac = np.inf
        else:
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

        # try:
        phi = newton_krylov(func, phi0)
        # except Exception as ex:
        #     print(ex)
        #     phi = np.array(np.nan)

        return phi.item()
