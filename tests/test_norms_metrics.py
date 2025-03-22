import pytest
from pymath import vector_l2_norm, matrix_frobenius_norm


def test_vector_norms():
    assert vector_l2_norm([3, 4]) == 5
    assert matrix_frobenius_norm([[1, 0], [0, 1]]) == pytest.approx(2 ** 0.5)
