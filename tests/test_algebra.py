import pytest
from pymath import solve_slae

def test_solve_slae_basic():
    A = [[2, 1], [1, 3]]
    b = [5, 10]
    expected = [1.0, 3.0]
    result = solve_slae(A, b)
    assert pytest.approx(result, rel=1e-5) == expected

def test_solve_slae_singular():
    A = [[1, 2], [2, 4]]
    b = [3, 6]
    with pytest.raises(ValueError) as exc_info: 
        solve_slae(A, b)
    assert "singular matrix" in str(exc_info.value).lower()