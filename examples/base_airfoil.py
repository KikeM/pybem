from pybem.models import Propeller, Section, BaseAirfoil
from pybem.bem import BladeElementMethod

# Define known sections
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

# Define propeller
B = 6
propeller = Propeller(B=B, sections=sections)

# Define flow conditions and BEM method
J = 0.2
bem = BladeElementMethod(J=J, propeller=propeller, tip_loss=False, hub_loss=False)

# Solve
bem.solve()

# Compute forces
CT, CQ = bem.integrate_forces()
