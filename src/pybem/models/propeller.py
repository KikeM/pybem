from scipy.interpolate import interp1d


class Propeller:
    """Propeller definition.

    Parameters
    ----------
    r_hub : float
    r_tip : float
    r_dist : array-like of floats
    beta_dist : array-like of floats
    chord_dist : array-like of floats
    """

    def __init__(self, r_hub, r_tip, r_dist, beta_dist, chord_dist):

        # Get blade bounds
        self.r_hub = r_hub
        self.r_tip = r_tip

        # Get airfoil properties
        self.dist_twist = beta_dist
        self.dist_chord = chord_dist
        self.dist_r = r_dist

        self.interpolant_twist = interp1d(r_dist, beta_dist)
        self.interpolant_chord = interp1d(r_dist, chord_dist)

        self.r = None
        self._beta = None
        self._chord = None

    def get_beta(self, r):
        """Airfoil twist at location r.

        Parameters
        ----------
        r: float
            Location based on r_loc.

        Returns
        -------
        beta: float
        """
        self.r = r

        _beta = self.interpolant_twist(r).item()

        self._beta = _beta

        return _beta

    def get_chord(self, r):
        """Airfoil chord at location r.

        Parameters
        ----------
        r: float
            Location based on r_loc.

        Returns
        -------
        chord: float
        """
        self.r = r

        _chord = self.interpolant_chord(r).item()

        self._chord = _chord

        return _chord
