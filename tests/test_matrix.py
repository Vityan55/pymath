import pytest
from pymath import determinant, matrix_multiply, matrix_inverse


def test_determinant():
    assert determinant([[1, 2], [3, 4]]) == -2
    assert determinant([[2]]) == 2


def test_matrix_multiplication():
    A = [[1, 2], [3, 4]]
    B = [[5, 6], [7, 8]]
    result = matrix_multiply(A, B)
    assert result[0] == [19, 22]
    assert result[1] == [43, 50]


def test_matrix_inverse():
    matrix = [[1, 2], [3, 4]]
    inv = matrix_inverse(matrix)
    product = matrix_multiply(matrix, inv)
    assert product[0] == pytest.approx([1.0, 0.0], abs=1e-3)
