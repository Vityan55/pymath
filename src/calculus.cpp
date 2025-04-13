#include "calculus.hpp"
namespace calculus {
    double find_extrema(const std::function<double(double)>& f, double a, double b, double step) {
        double prev = a;
        double curr = a + step;
        while (curr <= b) {
            if ((f(curr) - f(prev)) * (f(curr + step) - f(curr)) < 0)
                return curr;
            prev = curr;
            curr += step;
        }
        return -1.0;
    }

    double integrate(const std::function<double(double)>& f, double a, double b, int n) {
        double h = (b - a) / n;
        double result = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; ++i) {
            result += f(a + i * h);
        }
        return result * h;
    }

    double integrate_double(const std::function<double(double, double)>& f,
                            double ax, double bx, double ay, double by,
                            int nx, int ny) {
        double hx = (bx - ax) / nx;
        double hy = (by - ay) / ny;
        double result = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x = ax + i * hx;
                double y = ay + j * hy;
                result += f(x, y) * hx * hy;
            }
        }
        return result;
    }

    double integrate_triple(const std::function<double(double, double, double)>& f,
                            double ax, double bx, double ay, double by, double az, double bz,
                            int nx, int ny, int nz) {
        double hx = (bx - ax) / nx;
        double hy = (by - ay) / ny;
        double hz = (bz - az) / nz;
        double result = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    double x = ax + i * hx;
                    double y = ay + j * hy;
                    double z = az + k * hz;
                    result += f(x, y, z) * hx * hy * hz;
                }
            }
        }
        return result;
    }

    std::vector<std::pair<double, double>> solve_ode(const std::function<double(double, double)>& f,
                                                    double y0, double t0, double t1, double h) {
        double t = t0;
        double y = y0;
        std::vector<std::pair<double, double>> solution;
        solution.emplace_back(t, y);
        while (t < t1) {
            y += h * f(t, y);
            t += h;
            solution.emplace_back(t, y);
        }
        return solution;
    }
}