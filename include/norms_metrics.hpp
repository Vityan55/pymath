#pragma once
#include <vector>

namespace norms_metrics {
    double vector_l1_norm(const std::vector<double>& v);
    double vector_l2_norm(const std::vector<double>& v);
    double vector_linf_norm(const std::vector<double>& v);
    
    double matrix_frobenius_norm(const std::vector<std::vector<double>>& matrix);
    double matrix_l1_norm(const std::vector<std::vector<double>>& matrix);
    double matrix_linf_norm(const std::vector<std::vector<double>>& matrix);
    
    double euclidean_distance(
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );
    
    double manhattan_distance(
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );
    
    double chebyshev_distance(
        const std::vector<double>& v1,
        const std::vector<double>& v2
    );
}