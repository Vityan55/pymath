#pragma once
#include <vector>
#include <complex>

namespace transforms {
    std::vector<std::complex<double>> dft(
        const std::vector<std::complex<double>>& signal
    );
    
    std::vector<double> idft(
        const std::vector<std::complex<double>>& transform
    );
    
    std::vector<std::complex<double>> fft(
        const std::vector<std::complex<double>>& signal
    );
    
    std::pair<double, double> laplace_transform(
        double (*f)(double),
        double s_real,
        double s_imag,
        double t_max,
        double dt
    );
}