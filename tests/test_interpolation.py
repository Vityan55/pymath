import pytest
from pymath import lagrange_interpolation, spline_interpolation


def test_lagrange_interpolation():
    x = [0, 1, 2]
    y = [0, 1, 4]
    assert lagrange_interpolation(x, y, 1.5) == pytest.approx(2.25)


def test_spline_interpolation_boundary():
    x = [0, 1, 2, 3]
    y = [0, 1, 4, 9]
    assert spline_interpolation(x, y, 3) == pytest.approx(9, abs=1e-2)


def test_spline_interpolation_midpoint():
    x = [0, 2, 4]
    y = [1, 4, 9]
    assert spline_interpolation(x, y, 2) == pytest.approx(4, abs=1e-2)
