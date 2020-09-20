import pytest

from pybem.models import Propeller, Section, BaseAirfoil
from pybem.bem import BladeElementMethod


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


def test_no_losses(propeller):

    bem = BladeElementMethod(J=1, propeller=propeller)

    bem.solve()

    pass