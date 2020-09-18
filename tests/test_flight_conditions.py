import pytest

from pybem.core import FlightConditions

from numpy.testing import assert_allclose


def test_airspeed():

    conditions = FlightConditions(airspeed=1000, omega=1, altitude=1000)

    expected = 277.777777
    result = conditions.v

    assert_allclose(expected, result)


def test_omega():

    conditions = FlightConditions(airspeed=1000, omega=1, altitude=1000)

    expected = 0.1047197551
    result = conditions.omega

    assert_allclose(expected, result)


def test_density():

    conditions = FlightConditions(airspeed=1000, omega=1, altitude=1000)

    expected = 1.1116589850558272
    result = conditions.rho

    assert_allclose(expected, result)
