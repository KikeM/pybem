import pytest

from pybem import Propeller


@pytest.mark.xfail(
    raises=NotImplementedError, reason="Propeller definition needs to be updated."
)
def test_propeller():
    raise NotImplementedError
