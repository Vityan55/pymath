#include "visualization.hpp"
#include <iostream>
#include <iomanip>  

namespace visualization {
void plot_solution(const std::vector<std::pair<double, double>>& solution) {
    if (solution.empty()) {
        return;
    }

    double min_y = solution[0].second, max_y = solution[0].second;
    for (const auto& [t, y] : solution) {
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }

    int scale = 20;
    for (const auto& [t, y] : solution) {
        int pos = (max_y > min_y) ? int(scale * (y - min_y) / (max_y - min_y)) : 0;
        std::cout << std::fixed << std::setprecision(2) << t << " |" + std::string(pos, ' ') + "*" << std::endl;
    }
}
}