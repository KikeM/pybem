from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Blade Element Method
from pybem.bem.model import BladeElementMethod

# Core
from pybem.models.airfoils import Airfoil, BaseAirfoil
from pybem.models.flight_conditions import FlightConditions
from pybem.models.propeller import Propeller, Section
