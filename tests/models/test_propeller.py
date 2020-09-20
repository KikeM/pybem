from math import pi

import pytest
from numpy.testing import assert_allclose

from pybem.models.propeller import Propeller, Section


@pytest.fixture
def propeller():
    sections = [
        Section(name="Hub", r=0.3, beta=60, chord=0.2, airfoil=None),
        Section(name="Middle", r=0.6, beta=45, chord=0.35, airfoil=None),
        Section(name="Tip", r=1.2, beta=30, chord=0.5, airfoil=None),
    ]

    B = 6

    propeller = Propeller(B=B, sections=sections)

    return propeller


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


def test_section_solidity():

    section = Section(r=1.0, beta=0, chord=1.0, airfoil=None)

    result = section.solidity(1, 1.0)
    expected = 1.0 / (2.0 * pi)

    assert_allclose(expected, result)
