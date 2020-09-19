from math import pi

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pybem.bem.loss import compute_hub_loss, compute_tip_loss, func


def test_auxiliary_function():

    phi = 30.0  # deg
    phi = np.deg2rad(phi)
    x = 1.0
    a = 1.0 / 2.0

    result = func(phi=phi, x=x, a=a)

    expected = -1.0 / 2.0

    assert_allclose(expected, result)


@pytest.mark.parametrize(
    argnames="B", argvalues=[(1), (2)], ids=["One blade", "Multiple blades"]
)
def test_tip_loss_factor(B):

    phi = 30.0  # deg
    phi = np.deg2rad(phi)
    x = 1.0 / 2.0

    result = compute_tip_loss(phi=phi, r=x, B=B)

    expected = 2.0 / pi * np.arccos(np.exp(-B))

    assert_allclose(expected, result)


@pytest.mark.parametrize(
    argnames="B", argvalues=[(1), (2)], ids=["One blade", "Multiple blades"]
)
def test_hub_loss_factor(B):

    phi = 30.0  # deg
    phi = np.deg2rad(phi)
    x = 1.0 / 2.0
    r_hubR = 1.0 / 5.0

    result = compute_hub_loss(phi=phi, r=x, B=B, r_hubR=r_hubR)

    expected = 2.0 / pi * np.arccos(np.exp(-B * 3.0 / 5.0))

    assert_allclose(expected, result)
