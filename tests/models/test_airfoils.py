from math import pi

import pytest
from numpy.testing import assert_allclose, assert_array_almost_equal
import numpy as np

from pybem.models.airfoils import Airfoil, BaseAirfoil


@pytest.fixture
def polars():

    analytical = BaseAirfoil()

    alphas = np.linspace(0, 1, num=1000)
    polar_cl = []
    polar_cd = []

    for alpha in alphas:

        _cl = analytical.compute_cl(alpha)
        _cd = analytical.compute_cd(alpha=alpha)

        polar_cl.append(_cl)
        polar_cd.append(_cd)

    return alphas, np.array(polar_cl), np.array(polar_cd)


@pytest.mark.parametrize(
    argnames="alpha, expected_cl",
    argvalues=[(0.0, 0.0), (1, 2.0 * pi)],
    ids=["Zero angle of attack", "Unitary angle of attack"],
)
def test_analytical_airfoil_cl(alpha, expected_cl):

    airfoil = BaseAirfoil()

    result_cl = airfoil.compute_cl(alpha)

    deg2rad = pi / 180.0

    expected_cl = 2.0 * pi * (alpha * deg2rad)

    assert_allclose(expected_cl, result_cl)


def test_airfoil_analytical_interpolation(polars):

    alphas, expected_cl, expected_cd = polars

    polar_cl = np.column_stack((alphas, expected_cl))
    polar_cd = np.column_stack((alphas, expected_cd))

    empirical = Airfoil(polar_cl=polar_cl, polar_cd=polar_cd)
    result_cl = []
    result_cd = []

    for alpha in alphas:

        _cl = empirical.compute_cl(alpha=alpha)
        _cd = empirical.compute_cd(alpha=alpha)

        result_cl.append(_cl)
        result_cd.append(_cd)

    assert_array_almost_equal(expected_cl, result_cl)
    assert_array_almost_equal(expected_cd, result_cd)


def test_airfoil_unseen_interpolation(polars):

    alphas, polar_cl, polar_cd = polars

    polar_cl = np.column_stack((alphas, polar_cl))
    polar_cd = np.column_stack((alphas, polar_cd))

    empirical = Airfoil(polar_cl=polar_cl, polar_cd=polar_cd)

    alpha0 = 0.5
    result_cl = empirical.compute_cl(alpha0)
    result_cd = empirical.compute_cd(alpha=alpha0)

    deg2rad = pi / 180.0

    expected_cl = 2.0 * pi * (alpha0 * deg2rad)
    expected_cd = expected_cl ** 2.0

    assert_allclose(expected_cl, result_cl)
    assert_allclose(expected_cd, result_cd, rtol=1e-5)
