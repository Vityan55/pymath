#include "norms_metrics.hpp"
#include "math_utils.hpp"

namespace norms_metrics {

    double vector_l1_norm(const std::vector<double>& v) {
        double sum = 0;
        for (double x : v) {
            sum += std::abs(x);
        }
        return sum;
    }

    double vector_l2_norm(const std::vector<double>& v) {
        double sum = 0;
        for (double x : v) {
            sum += x * x;
        }
        return math_utils::newton_sqrt(sum);
    }

    double vector_linf_norm(const std::vector<double>& v) {
        return *std::max_element(v.begin(), v.end(), [](double a, double b) {
            return std::abs(a) < std::abs(b);
        });
    }

    double matrix_frobenius_norm(const std::vector<std::vector<double>>& matrix) {
        double sum = 0;
        for (const auto& row : matrix) {
            for (double x : row) {
                sum += x * x;
            }
        }
        return math_utils::newton_sqrt(sum);
    }

    double matrix_l1_norm(const std::vector<std::vector<double>>& matrix) {
        double max_sum = 0;
        size_t n_cols = matrix[0].size();
        for (size_t i = 0; i < n_cols; ++i) {
            double col_sum = 0;
            for (const auto& row : matrix) {
                col_sum += std::abs(row[i]);
            }
            max_sum = std::max(max_sum, col_sum);
        }
        return max_sum;
    }

    double matrix_linf_norm(const std::vector<std::vector<double>>& matrix) {
        double max_sum = 0;
        for (const auto& row : matrix) {
            double row_sum = 0;
            for (double x : row) {
                row_sum += std::abs(x);
            }
            max_sum = std::max(max_sum, row_sum);
        }
        return max_sum;
    }

    double euclidean_distance(const std::vector<double>& v1, const std::vector<double>& v2) {
        double sum = 0;
        for (size_t i = 0; i < v1.size(); ++i) {
            sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
        }
        return math_utils::newton_sqrt(sum, 1e-10);
    }

    double manhattan_distance(const std::vector<double>& v1, const std::vector<double>& v2) {
        double sum = 0;
        for (size_t i = 0; i < v1.size(); ++i) {
            sum += std::abs(v1[i] - v2[i]);
        }
        return sum;
    }

    double chebyshev_distance(const std::vector<double>& v1, const std::vector<double>& v2) {
        double max_diff = 0;
        for (size_t i = 0; i < v1.size(); ++i) {
            max_diff = std::max(max_diff, std::abs(v1[i] - v2[i]));
        }
        return max_diff;
    }

}
