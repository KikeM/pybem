from math import pi

from scipy.interpolate import interp1d


class BaseAirfoil:
    """Analytical airfoil implementation.

    Attributes
    ----------
    alpha : float
        Angle of attack in radians.
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
            Angle of attack in radians.

        Returns
        -------
        cl : float
        """

        cl = 2.0 * pi * alpha * self.cl_coeff

        self.alpha = alpha
        self.cl = cl

        return cl

    def compute_cd(self, alpha):
        """Compute drag coefficient.

        Parameters
        ----------
        alpha : float
            Angle of attack in radians.

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

    Parameters
    ----------
    alpha : array-like of floats
        Angles of attack for `polar_cl` in radians.
    polar_cl : array-like of floats
        Lift coefficients.
    polar_cd : array-like of floats
        Drag coefficients.
    """

    def __init__(self, alpha, polar_cl, polar_cd):

        super().__init__()

        # Store the data
        self.polar_cl = polar_cl
        self.polar_cd = polar_cd

        # Create interpolants
        self.interpolant_cl = interp1d(alpha, polar_cl)
        self.interpolant_cd = interp1d(polar_cl, polar_cd)

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

        if alpha is None:
            cd = self.interpolant_cd(cl).item()
        else:
            cl = self.compute_cl(alpha=alpha)
            cd = self.interpolant_cd(cl).item()

        # Update state
        self.cl = cl
        self.cd = cd

        return cd
