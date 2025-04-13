#pragma once
#include <utility>

namespace math_utils {

constexpr double PI = 3.141592653589793;

double newton_sqrt(double x, double epsilon = 1e-10);
double taylor_sin(double x, int terms = 15);
double taylor_cos(double x, int terms = 15);
double taylor_exp(double x, int terms = 20);
std::pair<double, double> complex_exp(std::pair<double, double> z);

}
