from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Core
from pybem.core import Propeller
from pybem.models.airfoils import Airfoil
from pybem.core import FlightConditions


# Blade Element Method
from pybem.bem import BladeElementMethod
