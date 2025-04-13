#include "montecarlo.hpp"
#include <random>

namespace montecarlo {

    double monte_carlo_integral(std::function<double(double)> f, double a, double b, int n) {
        double total = 0;
        unsigned seed = 1; 
        for (int i = 0; i < n; ++i) {
            seed = (seed * 1103515245 + 12345) % (2U << 30);
            double x = a + (b - a) * (static_cast<double>(seed) / (2U << 30));
            total += f(x);
        }
        return (b - a) * total / n;
    }

}
