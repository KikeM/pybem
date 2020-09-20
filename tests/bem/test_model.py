import pytest

from pybem.models import Propeller, Section, BaseAirfoil
from pybem.bem import BladeElementMethod

from numpy.testing import assert_allclose


@pytest.fixture
def propeller():

    sections = [
        Section(
            name="Hub",
            r=0.3,
            beta=60,
            chord=0.4,
            airfoil=BaseAirfoil(cl_coeff=1.0, cd_coeff=1e-2),
        ),
        Section(
            name="Middle",
            r=0.6,
            beta=45,
            chord=0.35,
            airfoil=BaseAirfoil(cl_coeff=0.85, cd_coeff=1e-3),
        ),
        Section(
            name="Tip",
            r=1.2,
            beta=30,
            chord=0.2,
            airfoil=BaseAirfoil(cl_coeff=0.5, cd_coeff=1e-3),
        ),
    ]

    B = 6

    propeller = Propeller(B=B, sections=sections)

    return propeller


def test_run_bem_no_losses(propeller):

    bem = BladeElementMethod(J=1, propeller=propeller)

    bem.solve()

    pass  # What could I assert?


def test_run_bem_with_losses(propeller):

    bem = BladeElementMethod(J=1, propeller=propeller, tip_loss=True, hub_loss=True)

    bem.solve()

    pass  # What could I assert?


def test_reproduce_bem_no_losses(propeller):
    """This is not a proper test, but it gives you some guarantees.

    With the actual state of the code, we are getting a result that seems
    to be correct for the type of phenomena we are solving.

    Perhaps if we linearise the equations we could validate the solver, but
    I haven't had time for the moment. So I check at least that for the same inputs,
    I am getting the same outputs from here on.

    If a bug is found in the future, this test would fail, but that is not
    necessarilly wrong. If the bug is proved to be correct, then the expected values
    need to be updated.
    """

    bem = BladeElementMethod(J=1, propeller=propeller)

    bem.N_SECTIONS = 10

    bem.solve()

    result_phi = bem.phi

    expected_phi = [
        65.02457613402427,
        61.53158664709961,
        58.06287232558151,
        54.680737977050924,
        52.70649811115369,
        50.648517749244384,
        48.57169404768725,
        46.52659452619687,
        44.549823308619374,
        42.66534592922986,
    ]

    assert_allclose(expected_phi, result_phi)


def test_reproduce_bem_with_losses(propeller):
    """This is not a proper test, but it gives you some guarantees.

    With the actual state of the code, we are getting a result that seems
    to be correct for the type of phenomena we are solving.

    Perhaps if we linearise the equations we could validate the solver, but
    I haven't had time for the moment. So I check at least that for the same inputs,
    I am getting the same outputs from here on.

    If a bug is found in the future, this test would fail, but that is not
    necessarilly wrong. If the bug is proved to be correct, then the expected values
    need to be updated.
    """

    bem = BladeElementMethod(J=1, propeller=propeller, tip_loss=True, hub_loss=True)

    bem.solve()

    bem.N_SECTIONS = 10

    bem.solve()

    result_phi = bem.phi

    expected_phi = [
        60.011868623160844,
        60.2642555093087,
        57.29588901630371,
        54.09221438071569,
        52.18007160086735,
        50.04596910839645,
        47.76712196831369,
        45.351112313287835,
        42.5933251947877,
        30.13349622416542,
    ]

    assert_allclose(expected_phi, result_phi)