#include "interpolation.hpp"

namespace interpolation {

double lagrange_interpolation(const std::vector<double>& x_vals, const std::vector<double>& y_vals, double x) {
    double result = 0.0;
    size_t n = x_vals.size();
    for (size_t i = 0; i < n; ++i) {
        double term = y_vals[i];
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j]);
            }
        }
        result += term;
    }
    return result;
}

double newton_interpolation(const std::vector<double>& x_vals, const std::vector<double>& y_vals, double x) {
    size_t n = x_vals.size();
    std::vector<double> coef = y_vals;
    for (size_t j = 1; j < n; ++j) {
        for (size_t i = n - 1; i >= j; --i) {
            coef[i] = (coef[i] - coef[i - 1]) / (x_vals[i] - x_vals[i - j]);
        }
    }
    double result = coef[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x - x_vals[i]) + coef[i];
    }
    return result;
}

double spline_interpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    size_t n = x.size();
    std::vector<double> h(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    std::vector<double> alpha(n, 0.0);
    for (size_t i = 1; i < n - 1; ++i) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    std::vector<double> l(n), mu(n), z(n), c(n), b(n - 1), d(n - 1);
    l[0] = 1.0;
    for (size_t i = 1; i < n - 1; ++i) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.0;

    for (int j = n - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    size_t i = n - 2;
    for (size_t j = 0; j < n - 1; ++j) {
        if (x[j] <= xi && xi <= x[j + 1]) {
            i = j;
            break;
        }
    }

    double dx = xi - x[i];
    return y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}

} 
