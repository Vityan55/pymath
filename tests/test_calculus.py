import pytest
from pymath import find_extrema, integrate, solve_ode

def test_find_extrema_quadratic():
    f = lambda x: x**2 - 4*x + 3
    extremum = find_extrema(f, 0.0, 3.0, epsilon=1e-5) 
    assert pytest.approx(extremum, abs=1e-5) == 2.0
def test_integrate_polynomial():
    f = lambda x: x**2
    result = integrate(f, 0, 2, n=1000)
    assert pytest.approx(result, rel=1e-3) == 8/3

def test_ode_exponential():
    f = lambda x, y: y 
    solution = solve_ode(f, y0=1.0, x0=0.0, x_end=1.0, h=0.01) 
    x_last, y_last = solution[-1]
    assert pytest.approx(y_last, rel=1e-2) == 2.718