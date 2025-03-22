import pytest
from pymath import find_extrema, integrate, solve_ode


def test_find_extrema_quadratic():
    f = lambda x: x ** 2 - 2 * x + 1
    assert find_extrema(f, 0, 3) == pytest.approx(1.0, abs=1e-3)


def test_integrate_polynomial():
    result = integrate(lambda x: x ** 2, 0, 2)
    assert result == pytest.approx(8 / 3, rel=1e-3)


def test_ode_exponential():
    f = lambda t, y: y
    solution = solve_ode(f, 1.0, 0.0, 1.0, h=0.01)
    assert solution[-1][1] == pytest.approx(2.718, rel=1e-2)
