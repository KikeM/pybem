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


class TestRun:
    def test_no_losses(self, propeller):

        bem = BladeElementMethod(J=1, propeller=propeller)

        bem.solve()

        pass  # What could I assert?

    def test_with_losses(self, propeller):

        bem = BladeElementMethod(J=1, propeller=propeller, tip_loss=True, hub_loss=True)

        bem.solve()

        pass  # What could I assert?


class TestReproduceOutputs:
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

    def test_solve_no_losses(self, propeller):

        bem = BladeElementMethod(J=1, propeller=propeller)

        bem.N_SECTIONS = 10

        bem.solve()

        result_phi = bem.phi_dist

        expected_phi = [
            65.02458658223601,
            61.531591303050206,
            58.062871180704974,
            54.68073133672458,
            52.70649067486621,
            50.64850667164986,
            48.571679570213234,
            46.5265770548006,
            44.549803327627394,
            42.66532393117762,
        ]

        assert_allclose(expected_phi, result_phi)

    def test_solve_with_losses(self, propeller):

        bem = BladeElementMethod(J=1, propeller=propeller, tip_loss=True, hub_loss=True)

        bem.N_SECTIONS = 10

        bem.solve()

        result_phi = bem.phi_dist

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

    def test_forces_no_losses(self, propeller):

        bem = BladeElementMethod(
            J=1, propeller=propeller, tip_loss=False, hub_loss=False
        )

        bem.solve()

        CT, CQ = bem.integrate_forces()

        result_forces = [CT, CQ]
        expected_forces = [-0.7608500036338367, -0.597588130498909]

        assert_allclose(expected_forces, result_forces)

    def test_forces_with_losses(self, propeller):

        bem = BladeElementMethod(J=1, propeller=propeller, tip_loss=True, hub_loss=True)

        bem.solve()

        CT, CQ = bem.integrate_forces()

        result_forces = [CT, CQ]
        expected_forces = [-0.6787816798433173, -0.513027309014405]

        assert_allclose(expected_forces, result_forces)