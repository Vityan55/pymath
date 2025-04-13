#pragma once
#include <vector>
#include <stdexcept>
#include <complex>

namespace matrix {

    using Complex = std::complex<double>;
    using Matrix = std::vector<std::vector<Complex>>;

    void check_square(const Matrix& matrix);

    Matrix multiply_scalar(const Matrix& matrix, const Complex& scalar);

    Complex determinant(const Matrix& matrix);
    Matrix hermitian_conjugate(const Matrix& matrix);
    Matrix multiply(const Matrix& A, const Matrix& B);
    Matrix exponential(const Matrix& matrix, int iterations = 10);
    Matrix power(const Matrix& matrix, int power);
    Matrix inverse(const Matrix& matrix);
    
    Matrix identity_matrix(size_t n);
    Matrix matrix_minor(const Matrix& matrix, size_t row, size_t col);
}