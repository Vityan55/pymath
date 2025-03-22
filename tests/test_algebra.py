import pytest
from pymath import solve_slae


def test_solve_slae_basic():
    A = [[2, 1], [1, 3]]
    b = [5, 10]
    assert solve_slae(A, b) == pytest.approx([1.0, 3.0], rel=1e-3)


def test_solve_slae_singular():
    with pytest.raises(ValueError):
        solve_slae([[1, 1], [1, 1]], [2, 2])
