#pragma once
#include <vector>

namespace interpolation {
    double lagrange_interpolation(
        const std::vector<double>& x_vals,
        const std::vector<double>& y_vals,
        double x
    );
    
    double newton_interpolation(
        const std::vector<double>& x_vals,
        const std::vector<double>& y_vals,
        double x
    );
    
    double spline_interpolation(
        const std::vector<double>& x,
        const std::vector<double>& y, 
        double xi
    );
}