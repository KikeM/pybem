from math import pi

import numpy as np
from scipy.interpolate import interp1d


class BaseAirfoil:
    """Analytical airfoil implementation.

    Parameters
    ----------
    cl_coeff : float
        Constant to multiply the analytical lift coefficient.
    cd_coeff : float
        Constant to multiply the analytical drag coefficient.

    Attributes
    ----------
    alpha : float
        Angle of attack in degrees.
    cl : float
        Lift coefficient.
    cd : float
        Drag coefficient.
    """

    def __init__(self, cl_coeff=1.0, cd_coeff=1.0):

        self.alpha = None
        self.cl = None
        self.cd = None

        self.cl_coeff = cl_coeff
        self.cd_coeff = cd_coeff

    def compute_cl(self, alpha):
        """Compute lift coefficient.

        Parameters
        ----------
        alpha : float
            Angle of attack in degrees.

        Returns
        -------
        cl : float
        """

        alpha = np.deg2rad(alpha)

        cl = 2.0 * pi * alpha * self.cl_coeff

        self.alpha = alpha
        self.cl = cl

        return cl

    def compute_cd(self, alpha):
        """Compute drag coefficient.

        Parameters
        ----------
        alpha : float
            Angle of attack in degrees.

        Returns
        -------
        cd : float
        """

        cl = self.compute_cl(alpha=alpha)

        cd = self.cd_coeff * (cl ** 2.0)

        self.alpha = alpha
        self.cl = cl
        self.cd = cd

        return cd


class Airfoil(BaseAirfoil):
    """Airfoil section aerodynamic coefficients.

    Attributes
    ----------
    polar_cl : array-like of floats
        [[alpha_0, cl_0], [alpha_1, cl_1], ..., [alpha_N, cl_N]]
    polar_cd : array-like of floats
        [[alpha_0, cd_0], [alpha_1, cd_1], ..., [alpha_N, cd_N]]
    interpolant_cl : scipy.interpolate.interp1-like object
    interpolant_cd : scipy.interpolate.interp1-like object
    """

    def __init__(self, polar_cl, polar_cd):
        """
        Parameters
        ----------
        polar_cl : array-like of floats
            Lift coefficients.
            [[alpha_0, cl_0], [alpha_1, cl_1], ..., [alpha_N, cl_N]]
        polar_cd : array-like of floats
            Drag coefficients.
            - [[alpha_0, cd_0], [alpha_1, cd_1], ..., [alpha_N, cd_N]]
            or
            - [[cl_0, cd_0], [cl_1, cd_1], ..., [cl_N, cd_N]]


        Notes
        -----
        The angles for interpolation for the Cl and the Cd do not have
        to be necessarily the same.
        """

        super().__init__()

        # Store the data
        self.polar_cl = polar_cl
        self.polar_cd = polar_cd

        # Create interpolants
        self.interpolant_cl = interp1d(polar_cl[:, 0], polar_cl[:, 1])
        self.interpolant_cd = interp1d(polar_cd[:, 0], polar_cd[:, 1])

    def compute_cl(self, alpha):
        """Compute lift coefficient.

        Parameters
        ----------
        alpha: float

        Returns
        -------
        cl: float
            Interpolation based on lift polar.
        """

        # Compute
        cl = self.interpolant_cl(alpha).item()

        # Update state
        self.alpha = alpha
        self.cl = cl

        return cl

    def compute_cd(self, cl=None, alpha=None):
        """Compute drag coefficient.

        Parameters
        ----------
        alpha : float
            Angle of attack.
        cl : float
            Lift coefficient.

        Returns
        -------
        cd: float
            Interpolation based on drag polar.
        """

        #  If the drag polar is given in terms of cl
        if alpha is None:
            cd = self.interpolant_cd(cl).item()
            self.cl = cl

        else:
            cd = self.interpolant_cd(alpha).item()
            self.alpha = alpha

        # Update state
        self.cd = cd

        return cd
