from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Core
from pybem.models.propeller import Propeller
from pybem.models.airfoils import Airfoil
from pybem.models.flight_conditions import FlightConditions


# Blade Element Method
from pybem.bem.model import BladeElementMethod
