from pymath import determinant, matrix_multiply, matrix_inverse
import pytest

def test_determinant():
    matrix = [[1.0, 2.0], [3.0, 4.0]]
    result = determinant(matrix)
    
    assert isinstance(result, complex)
    assert result.real == pytest.approx(-2.0)
    assert result.imag == pytest.approx(0.0)

def test_matrix_multiplication():
    A = [[1.0, 2.0], [3.0, 4.0]]
    B = [[5.0, 6.0], [7.0, 8.0]]
    expected = [[19.0, 22.0], [43.0, 50.0]]
    
    result = matrix_multiply(A, B)
    
    for i in range(2):
        for j in range(2):
            assert isinstance(result[i][j], complex)
            assert result[i][j].real == pytest.approx(expected[i][j])
            assert result[i][j].imag == 0.0

def test_matrix_inverse():
    matrix = [[1.0, 2.0], [3.0, 4.0]]
    inverse = matrix_inverse(matrix)
    expected = [[-2.0, 1.0], [1.5, -0.5]]
    
    for i in range(2):
        for j in range(2):
            assert isinstance(inverse[i][j], complex)
            assert inverse[i][j].real == pytest.approx(expected[i][j], abs=1e-6)
            assert inverse[i][j].imag == pytest.approx(0.0, abs=1e-6)