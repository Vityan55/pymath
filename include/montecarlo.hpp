#pragma once
#include <functional>

namespace montecarlo {

    double monte_carlo_integral(std::function<double(double)> f, double a, double b, int n = 1000);

}

