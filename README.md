# pybem

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.org/KikeM/pybem.svg?branch=master)](https://travis-ci.org/KikeM/pybem)

Blade Element Method implementation for propeller calculations.

## Installation

To run it as a user, simply invoke `pip`:
```bash
pip install .
```

### For developers

If you want to contribute to the library, or tweak it to your own needs, install it in developer mode, including the development libraries:
```bash
pip install -e . --requirement requirements-dev.txt
```

## Quickstart

Running the code consists of easy and uncoupled steps:
1. Declare the airfoil sections with their corresponding geometrical definition.
1. Create a propeller by putting together the sections and the number of blades.
1. Create a solver by putting together the propeller and a advance ratio.
1. Solve the flow.
1. Compute the force and torque coefficients.

Here is an example with an airfoil defined by an analytical lift and drag polars. 

```python
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
```