import pytest
from pymath import monte_carlo_integral


def test_monte_carlo():
    result = monte_carlo_integral(lambda x: x ** 2, 0, 2, n=100000)
    assert result == pytest.approx(8 / 3, rel=0.1)
