#include "algebra.hpp"
#include <cmath> 

namespace algebra {
    std::vector<double> solve_slae(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
        size_t n = A.size();
        std::vector<std::vector<double>> matrix = A;
        std::vector<double> result = b;

        for (size_t i = 0; i < n; ++i) {
            size_t max_row = i;
            for (size_t k = i + 1; k < n; ++k) {
                if (std::abs(matrix[k][i]) > std::abs(matrix[max_row][i])) {
                    max_row = k;
                }
            }

            if (std::abs(matrix[max_row][i]) < 1e-10) {
                throw std::invalid_argument("Singular matrix detected");
            }

            std::swap(matrix[i], matrix[max_row]);
            std::swap(result[i], result[max_row]);

            for (size_t k = i + 1; k < n; ++k) {
                double factor = matrix[k][i] / matrix[i][i];
                result[k] -= factor * result[i];
                for (size_t j = i; j < n; ++j) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }

        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (size_t j = i + 1; j < n; ++j) {
                sum += matrix[i][j] * x[j];
            }
            x[i] = (result[i] - sum) / matrix[i][i];
        }

        return x;
    }
}