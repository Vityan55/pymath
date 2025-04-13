#include "math_utils.hpp"
#include <stdexcept>
#include <cmath>

namespace math_utils {

double newton_sqrt(double x, double epsilon) {
    if (x < 0.0) {
        throw std::domain_error("Negative number under the square root");
    }
    double guess = x;
    while (std::abs(guess * guess - x) > epsilon) {
        guess = (guess + x / guess) / 2.0;
    }
    return guess;
}

double taylor_sin(double x, int terms) {
    x = std::fmod(x, 2 * PI);
    double result = 0.0;
    double power = x;
    int sign = 1;
    for (int n = 1; n <= terms; ++n) {
        result += sign * power;
        power *= x * x / ((2 * n) * (2 * n + 1));
        sign *= -1;
    }
    return result;
}

double taylor_cos(double x, int terms) {
    x = std::fmod(x, 2 * PI);
    double result = 0.0;
    double power = 1.0;
    int sign = 1;
    for (int n = 0; n < terms; ++n) {
        result += sign * power;
        power *= x * x / ((2 * n + 1) * (2 * n + 2));
        sign *= -1;
    }
    return result;
}

double taylor_exp(double x, int terms) {
    double result = 1.0;
    double term = 1.0;
    for (int n = 1; n <= terms; ++n) {
        term *= x / n;
        result += term;
    }
    return result;
}

std::pair<double, double> complex_exp(std::pair<double, double> z) {
    double a = z.first;
    double b = z.second;
    double exp_a = taylor_exp(a);
    double real = taylor_cos(b) * exp_a;
    double imag = taylor_sin(b) * exp_a;
    return {real, imag};
}

}
