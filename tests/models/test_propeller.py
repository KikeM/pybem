from math import pi

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pybem.models.airfoils import BaseAirfoil
from pybem.models.propeller import Propeller, Section

TOL = 1e-4


@pytest.fixture
def propeller():

    sections = [
        Section(
            name="Hub",
            r=0.3,
            beta=60,
            chord=0.2,
            airfoil=BaseAirfoil(cl_coeff=1.0, cd_coeff=1.0),
        ),
        Section(
            name="Middle",
            r=0.6,
            beta=45,
            chord=0.35,
            airfoil=BaseAirfoil(cl_coeff=2.0, cd_coeff=2.0),
        ),
        Section(
            name="Tip",
            r=1.2,
            beta=30,
            chord=0.5,
            airfoil=BaseAirfoil(cl_coeff=3.0, cd_coeff=3.0),
        ),
    ]

    B = 6

    propeller = Propeller(B=B, sections=sections)

    return propeller


@pytest.fixture
def section():

    return Section(r=1.0, beta=0, chord=1.0, airfoil=BaseAirfoil())


def test_section_solidity(section):

    result = section.solidity(1, 1.0)

    expected = 1.0 / (2.0 * pi)

    assert_allclose(expected, result)


def test_section_cl(section):

    alpha = 1.0
    deg2rad = pi / 180.0

    result = section.cl(alpha=alpha)

    expected = 2.0 * pi * alpha * deg2rad

    assert_allclose(expected, result)


def test_section_cd(section):

    alpha = 1.0
    deg2rad = pi / 180.0

    result = section.cd(alpha=1.0)
    expected = (2.0 * pi * alpha * deg2rad) ** 2.0

    assert_allclose(expected, result)


def test_propeller_init_radius(propeller):

    result = [propeller.R_hub, propeller.R]
    expected = [0.3, 1.2]

    assert_allclose(expected, result)


def test_propeller_compute_chord(propeller):

    result = propeller.compute_chord(0.9)
    expected = 0.47

    assert_allclose(expected, result)


def test_propeller_compute_beta(propeller):

    result = propeller.compute_beta(0.9)
    expected = 33.0

    assert_allclose(expected, result)


def test_propeller_find_bracket_left(propeller):

    result_sections = propeller._find_bracket(0.3)

    expected_sections = (propeller.sections[0], propeller.sections[1])

    assert expected_sections == result_sections


def test_propeller_find_bracket_right(propeller):

    result_sections = propeller._find_bracket(0.9)

    expected_sections = (propeller.sections[1], propeller.sections[2])

    assert expected_sections == result_sections


def test_propeller_solidity(propeller):

    result_sigma = propeller.compute_solidity(0.9)

    expected_sigma = 0.37401411626595404

    assert_allclose(expected_sigma, result_sigma)


def test_propeller_compute_cl(propeller):

    alpha = 1.0
    deg2rad = pi / 180.0

    result_cl = propeller.compute_cl(r=0.75, alpha=alpha)
    expected_cl = 2.0 * pi * (alpha * deg2rad) * (3.0 + 2.0) / 2.0

    assert_allclose(expected_cl, result_cl)


def test_propeller_compute_cd(propeller):

    alpha = 1.0
    deg2rad = pi / 180.0

    result_cd = propeller.compute_cd(r=0.75, alpha=alpha)

    expected_cd = (
        2.0
        * pi ** 2.0
        * (alpha * deg2rad) ** 2.0
        * (3.0 * 3.0 ** 2.0 + 2.0 * 2.0 ** 2.0)
    )

    assert_allclose(expected_cd, result_cd)


@pytest.fixture
def R0_c0():
    return 0.3, 0.1


@pytest.fixture
def R_c():
    return 1.2, 0.5


def linear_chord(r, R, R0, c, c0):
    return r * (c - c0) / (R - R0) + (R * c0 - R0 * c) / (R - R0)


@pytest.fixture
def propeller_linear_chord(R_c, R0_c0):

    R, c = R_c
    R0, c0 = R0_c0

    N = 10

    sections = []
    for idx, r in enumerate(np.linspace(R0, R, num=N)):

        name = "S" + str(idx)
        chord = linear_chord(r=r, R=R, R0=R0, c0=c0, c=c)

        _section = Section(
            name=name,
            r=r,
            beta=0,
            chord=chord,
            airfoil=BaseAirfoil(),
        )

        sections.append(_section)

    B = 2

    propeller = Propeller(B=B, sections=sections)

    return propeller


@pytest.fixture
def expected_area(R_c, R0_c0):

    R, c = R_c
    R0, c0 = R0_c0

    expected_area = (1.0 / 2.0) * (-R + R0) ** 2.0 * (c + c0) / (R - R0)

    return expected_area


@pytest.fixture
def expected_aspect_ratio(expected_area, R_c):

    R, c = R_c

    AR = R ** 2.0 / expected_area

    return AR


def test_area(expected_area, propeller_linear_chord):

    result_area = propeller_linear_chord.area

    assert_allclose(expected_area, result_area, atol=TOL, rtol=TOL)


def test_aspect_ratio(expected_aspect_ratio, propeller_linear_chord):

    result_aspect_ratio = propeller_linear_chord.aspect_ratio

    assert_allclose(expected_aspect_ratio, result_aspect_ratio, atol=TOL, rtol=TOL)
