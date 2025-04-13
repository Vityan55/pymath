#pragma once
#include <vector>
#include <stdexcept>

namespace algebra {
    std::vector<double> solve_slae(
        const std::vector<std::vector<double>>& A, 
        const std::vector<double>& b
    );
}