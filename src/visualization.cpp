#include "visualization.hpp"
#include <vector>
#include <cstdio>
#include <algorithm>
#include <string>

namespace visualization {
// Plots a solution as an ASCII graph with time (t) and value (y)
// Parameters:
//   solution - Vector of (t, y) pairs representing the solution
void plot_solution(const std::vector<std::pair<double, double>>& solution) {
    if (solution.empty()) {
        return; // Silently return for empty solution
    }

    // Find min and max y for scaling
    double min_y = solution[0].second;
    double max_y = solution[0].second;
    for (const auto& [t, y] : solution) {
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
    }

    // Plot each point
    const int scale = 20; // Maximum width for the plot
    for (const auto& [t, y] : solution) {
        // Calculate position of *
        int pos = 0;
        if (max_y > min_y + 1e-10) { // Avoid division by zero
            pos = static_cast<int>(scale * (y - min_y) / (max_y - min_y));
            pos = std::max(0, std::min(pos, scale));
        }

        // Create spaces string
        std::string spaces(pos, ' ');

        // Output using printf with Windows-compatible newline
        std::printf("%.2f |%s *\r\n", t, spaces.c_str());
        std::fflush(stdout); // Ensure output is flushed
    }
}
} // namespace visualization