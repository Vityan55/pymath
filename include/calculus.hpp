#pragma once
#include <functional>
#include <vector>
#include <utility>

namespace calculus {

double find_extrema(const std::function<double(double)>& f, double a, double b, double step = 1e-5);

double integrate(const std::function<double(double)>& f, double a, double b, int n = 1000);

double integrate_double(const std::function<double(double, double)>& f,
                        double ax, double bx, double ay, double by,
                        int nx = 100, int ny = 100);

double integrate_triple(const std::function<double(double, double, double)>& f,
                        double ax, double bx, double ay, double by, double az, double bz,
                        int nx = 20, int ny = 20, int nz = 20);

std::vector<std::pair<double, double>> solve_ode(const std::function<double(double, double)>& f,
                                                 double y0, double t0, double t1, double h = 0.01);

}
