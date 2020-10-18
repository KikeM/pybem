from bisect import bisect_left
from functools import partial
from math import pi

import numpy as np
from scipy.integrate import simps as integrate
from scipy.interpolate import interp1d


class Section:
    """Propeller blade section.

    This class is used by the `Propeller` class, to build the distribution of
    chord, twist and solidity across the blade.

    Attributes
    ----------
    r : float
        Dimensionless radius.
    beta : float
        Geometrical twist.
    chord : float
        Airfoil section chord.
    solidity : float
    airfoil : AirfoilBase-like object
    """

    def __init__(self, r, beta, chord, airfoil, name=None):

        self.r = r
        self.beta = beta
        self.chord = chord
        self.airfoil = airfoil
        self.name = name

    def __repr__(self):
        return (
            f"Section {self.name} (r = {self.r}, c = {self.chord}, beta = {self.beta})"
        )

    def cl(self, alpha):
        """Section lift coefficient.

        Parameters
        ----------
        alpha : float

        Returns
        -------
        cl : float
        """

        _cl = self.airfoil.compute_cl(alpha=alpha)

        return _cl

    def cd(self, alpha):
        """Section drag coefficient.

        Parameters
        ----------
        alpha : float

        Returns
        -------
        cd : float
        """

        _cd = self.airfoil.compute_cd(alpha=alpha)

        return _cd

    def solidity(self, B, R):
        """Section solidity.

        Parameters
        ----------
        B : int
            Number of blades in the propeller.
        R : float
            Propeller radius.

        Returns
        -------
        sigma : float
        """

        NUM = B * self.chord
        DEN = 2.0 * pi * R

        sigma = NUM / DEN

        return sigma


class Propeller:
    """Propeller definition.

    Parameters
    ----------
    B : int
        Number of blades.
    sections : list of Sections-like objects

    Attributes
    ----------
    R : float
    R_hub : float
    sections : sorted list of Sections-like objects
    radii : list-like of floats
        Distribution of dimensionless radius distribution
    """

    def __init__(self, B, sections):

        self.B = B

        # Get blade bounds
        self.R = max(section.r for section in sections)
        self.R_hub = min(section.r for section in sections)

        # Sort sections by ascending order of radius
        func_selector = lambda section: section.r
        _sections = sorted(sections, key=func_selector)

        self.sections = []
        dist_r = []
        dist_beta = []
        dist_chord = []
        for section in _sections:

            # Make radius dimensionless
            section.r = section.r / self.R

            dist_r.append(section.r)
            dist_beta.append(section.beta)
            dist_chord.append(section.chord)

            # Save
            self.sections.append(section)

        self._interpolant_chord = interp1d(dist_r, dist_chord)
        self._interpolant_beta = interp1d(dist_r, dist_beta)

        self.radii = dist_r

    def _find_bracket(self, r):
        """For a given dimensionless location `r`, find the
        two surrounding airfoil Sections.

        Parameters
        ----------
        r : float
            Dimensionless radius.

        Returns
        -------
        left : Section-like object.
        right : Section-like object.
        """

        _radius = [section.r for section in self.sections]

        idx = bisect_left(_radius, r)

        left = self.sections[idx - 1]
        right = self.sections[idx]

        return left, right

    def compute_chord(self, r):
        """Compute section chord via interpolation between sections.

        Parameters
        ----------
        r : float

        Returns
        -------
        chord : float
        """
        chord = self._interpolant_chord(r).item()

        return chord

    def compute_beta(self, r):
        """Compute section beta via interpolation between sections.

        Parameters
        ----------
        r : float

        Returns
        -------
        beta : float
        """
        beta = self._interpolant_beta(r).item()

        return beta

    def compute_solidity(self, r):
        """Compute section solidity via interpolation between sections.

        Parameters
        ----------
        r : float

        Returns
        -------
        sgima : float
        """

        left, right = self._find_bracket(r)

        radii = [left.r, right.r]
        solidities = [
            left.solidity(B=self.B, R=self.R),
            right.solidity(B=self.B, R=self.R),
        ]
        f = interp1d(radii, solidities)

        sigma = f(r).item()

        return sigma

    def compute_cl(self, r, alpha):
        """Compute section lift coefficient via interpolation between sections.

        Parameters
        ----------
        r : float
        alpha : float

        Returns
        -------
        cl : float
        """

        left, right = self._find_bracket(r)

        radii = [left.r, right.r]
        coeffs = [
            left.cl(alpha=alpha),
            right.cl(alpha=alpha),
        ]

        f = interp1d(radii, coeffs)

        cl = f(r).item()

        return cl

    def compute_cd(self, r, alpha):
        """Compute section drag coefficient via interpolation between sections.

        Parameters
        ----------
        r : float
        alpha : float

        Returns
        -------
        cd : float
        """

        left, right = self._find_bracket(r)

        radii = [left.r, right.r]
        coeffs = [
            left.cd(alpha=alpha),
            right.cd(alpha=alpha),
        ]

        f = interp1d(radii, coeffs)

        cd = f(r).item()

        return cd

    def _dSdr(self, r):

        R = self.R

        return R * self.compute_chord(r)

    @property
    def area(self):

        space = np.linspace(self.R_hub, self.R, num=1000)
        space /= self.R

        _integrate = partial(integrate, x=space)

        dSdr = list(map(self._dSdr, space))

        S = _integrate(dSdr)

        return S

    @property
    def aspect_ratio(self):

        S = self.area
        AR = self.R ** 2.0 / S

        return AR