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

        bem = BladeElementMethod(_lambda=1, propeller=propeller)

        bem.solve()

        pass  # What could I assert?

    def test_with_losses(self, propeller):

        bem = BladeElementMethod(
            _lambda=1, propeller=propeller, tip_loss=True, hub_loss=True
        )

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

    TOL = 1e-3

    _lambda = 0.2

    def test_solve_no_losses(self, propeller):

        bem = BladeElementMethod(_lambda=self._lambda, propeller=propeller)

        bem.N_SECTIONS = 10

        bem.solve()

        result_phi = bem.phi_dist

        expected_phi = [
            53.7559900627035,
            46.51360450433435,
            40.24566201204451,
            34.82356412209791,
            31.020411856568106,
            27.60175003439307,
            24.517814495030052,
            21.727980253100842,
            19.200977399649037,
            16.913981930219638,
        ]

        assert_allclose(expected_phi, result_phi, rtol=self.TOL, atol=self.TOL)

    def test_solve_with_losses(self, propeller):

        bem = BladeElementMethod(
            _lambda=self._lambda, propeller=propeller, tip_loss=True, hub_loss=True
        )

        bem.N_SECTIONS = 10

        bem.solve()

        result_phi = bem.phi_dist

        expected_phi = [
            59.98375776054978,
            47.80654198333023,
            40.79851063596392,
            35.07657290878467,
            31.187956608585395,
            27.78122466482189,
            24.817388223831337,
            22.331258522446035,
            20.59650318210961,
            29.83229695496505,
        ]

        assert_allclose(expected_phi, result_phi, rtol=self.TOL, atol=self.TOL)

    def test_forces_no_losses(self, propeller):

        bem = BladeElementMethod(
            _lambda=self._lambda, propeller=propeller, tip_loss=False, hub_loss=False
        )

        bem.solve()

        CT, CQ = bem.integrate_forces()

        result_forces = [CT, CQ]
        expected_forces = [0.06015230682040654, 0.02035913185617931]

        assert_allclose(expected_forces, result_forces, atol=self.TOL, rtol=self.TOL)

    def test_forces_with_losses(self, propeller):

        bem = BladeElementMethod(
            _lambda=self._lambda, propeller=propeller, tip_loss=True, hub_loss=True
        )

        bem.solve()

        CT, CQ = bem.integrate_forces()

        result_forces = [CT, CQ]
        expected_forces = [0.05496400990860528, 0.019365623278049134]

        assert_allclose(expected_forces, result_forces, atol=self.TOL, rtol=self.TOL)