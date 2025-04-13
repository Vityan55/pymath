#include "matrix.hpp"
#include <cmath>
#include <algorithm>

namespace matrix {

    void check_square(const Matrix& matrix) {
        if (matrix.empty() || matrix.size() != matrix[0].size()) {
            throw std::invalid_argument("Matrix must be square");
        }
    }

    Matrix multiply_scalar(const Matrix& matrix, const Complex& scalar) {
        Matrix result = matrix;
        for (auto& row : result) {
            for (auto& element : row) {
                element *= scalar;
            }
        }
        return result;
    }

    Matrix identity_matrix(size_t n) {
        Matrix I(n, std::vector<Complex>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            I[i][i] = 1.0;
        }
        return I;
    }

    Matrix matrix_minor(const Matrix& matrix, size_t row, size_t col) {
        Matrix minor;
        for (size_t i = 0; i < matrix.size(); ++i) {
            if (i == row) continue;
            std::vector<Complex> new_row;
            for (size_t j = 0; j < matrix[i].size(); ++j) {
                if (j != col) new_row.push_back(matrix[i][j]);
            }
            minor.push_back(new_row);
        }
        return minor;
    }

    Complex determinant(const Matrix& matrix) {
        check_square(matrix);
        if (matrix.size() == 1) return matrix[0][0];
        
        Complex det = 0.0;
        for (size_t col = 0; col < matrix[0].size(); ++col) {
            Complex sign = (col % 2 == 0) ? 1.0 : -1.0;
            det += matrix[0][col] * sign * determinant(matrix_minor(matrix, 0, col));
        }
        return det;
    }

    Matrix hermitian_conjugate(const Matrix& matrix) {
        Matrix result(matrix[0].size(), std::vector<Complex>(matrix.size()));
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[i].size(); ++j) {
                result[j][i] = std::conj(matrix[i][j]);
            }
        }
        return result;
    }

    Matrix multiply(const Matrix& A, const Matrix& B) {
        check_square(A);
        check_square(B);
        size_t n = A.size();
        Matrix result(n, std::vector<Complex>(n, 0.0));
        
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }

    Matrix exponential(const Matrix& matrix, int iterations) {
        check_square(matrix);
        size_t n = matrix.size();
        Matrix result = identity_matrix(n);
        Matrix current_term = result;
        
        for (int k = 1; k <= iterations; ++k) {
            current_term = multiply(current_term, matrix);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    result[i][j] += current_term[i][j] / static_cast<double>(k);
                }
            }
        }
        return result;
    }

    Matrix power(const Matrix& matrix, int power) {
        check_square(matrix);
        if (power == 0) return identity_matrix(matrix.size());
        
        Matrix result = identity_matrix(matrix.size());
        Matrix base = matrix;
        int exponent = abs(power);
        
        while (exponent > 0) {
            if (exponent % 2 == 1) {
                result = multiply(result, base);
            }
            base = multiply(base, base);
            exponent /= 2;
        }
        return (power < 0) ? inverse(result) : result;
    }

    Matrix inverse(const Matrix& matrix) {
        check_square(matrix);
        size_t n = matrix.size();
        Complex det = determinant(matrix);
        if (std::abs(det) < 1e-10) {
            throw std::invalid_argument("Matrix is singular");
        }
        
        Matrix adjugate(n, std::vector<Complex>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                Complex sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
                adjugate[j][i] = sign * determinant(matrix_minor(matrix, i, j));
            }
        }
        return multiply_scalar(adjugate, Complex(1.0)/det);
    }
}