#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "algebra.hpp"
#include "calculus.hpp"
#include "interpolation.hpp"
#include "math_utils.hpp"
#include "matrix.hpp"
#include "montecarlo.hpp"
#include "norms_metrics.hpp"
#include "transforms.hpp"
#include "visualization.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pymath, m) {

    m.def("solve_slae", &algebra::solve_slae, 
        "Solve system of linear equations",
        py::arg("A"), 
        py::arg("b")
    );

    m.def("find_extrema", &calculus::find_extrema, 
        "Find extrema using numerical differentiation",
        py::arg("f"), py::arg("a"), py::arg("b"), py::arg("epsilon") = 1e-5
    );
    m.def("integrate", &calculus::integrate, 
        "Compute definite integral (trapezoidal rule)",
        py::arg("f"), py::arg("a"), py::arg("b"), py::arg("n") = 1000
    );
    m.def("integrate_double", &calculus::integrate_double, 
        "Compute double integral",
        py::arg("f"), py::arg("ax"), py::arg("bx"), py::arg("ay"), py::arg("by"), 
        py::arg("nx") = 100, py::arg("ny") = 100
    );
    m.def("integrate_triple", &calculus::integrate_triple, 
        "Compute triple integral",
        py::arg("f"), py::arg("ax"), py::arg("bx"), py::arg("ay"), py::arg("by"), 
        py::arg("az"), py::arg("bz"), py::arg("nx") = 20, py::arg("ny") = 20, 
        py::arg("nz") = 20
    );
    m.def("solve_ode", &calculus::solve_ode, 
        "Solve ODE using Euler's method",
        py::arg("f"), py::arg("y0"), py::arg("x0"), py::arg("x_end"), py::arg("h") = 0.01
    );
    
    m.def("newton_sqrt", &math_utils::newton_sqrt, "Square root using Newton's method",
        py::arg("x"), py::arg("epsilon") = 1e-10);
    m.def("taylor_sin", &math_utils::taylor_sin, "Sine via Taylor series",
        py::arg("x"), py::arg("terms") = 15);
    m.def("taylor_cos", &math_utils::taylor_cos, "Cosine via Taylor series",
        py::arg("x"), py::arg("terms") = 15);
    m.def("taylor_exp", &math_utils::taylor_exp, "Exponential via Taylor series",
        py::arg("x"), py::arg("terms") = 20);

    m.def("determinant", &matrix::determinant, 
        "Compute determinant of a complex matrix",
        py::arg("matrix"));
    m.def("hermitian_conjugate", &matrix::hermitian_conjugate,
        "Compute Hermitian conjugate of a matrix",
        py::arg("matrix"));
    m.def("multiply_scalar", &matrix::multiply_scalar,
        "Multiply matrix by scalar",
        py::arg("matrix"), py::arg("scalar"));
    m.def("matrix_multiply", &matrix::multiply,
        "Multiply two complex matrices",
        py::arg("A"), py::arg("B"));
    m.def("matrix_exponential", &matrix::exponential,
        "Compute matrix exponential",
        py::arg("matrix"), py::arg("iterations") = 10);
    m.def("matrix_power", &matrix::power,
        "Raise matrix to integer power",
        py::arg("matrix"), py::arg("power"));
    m.def("matrix_inverse", &matrix::inverse,
        "Compute matrix inverse",
        py::arg("matrix"));

    m.def("monte_carlo_integral", &montecarlo::monte_carlo_integral, 
        "Estimate integral using Monte Carlo method", 
        py::arg("f"), py::arg("a"), py::arg("b"), py::arg("n") = 1000);

    m.def("lagrange_interpolation", &interpolation::lagrange_interpolation,
        "Lagrange interpolation",
        py::arg("x_vals"), py::arg("y_vals"), py::arg("x"));
    m.def("newton_interpolation", &interpolation::newton_interpolation,
        "Newton interpolation",
        py::arg("x_vals"), py::arg("y_vals"), py::arg("x"));
    m.def("spline_interpolation", &interpolation::spline_interpolation,
        "Cubic spline interpolation",
        py::arg("x"), py::arg("y"), py::arg("xi"));

    m.def("vector_l1_norm", &norms_metrics::vector_l1_norm,
        "Compute L1 norm of vector");
    m.def("vector_l2_norm", &norms_metrics::vector_l2_norm,
        "Compute L2 norm of vector");
    m.def("vector_linf_norm", &norms_metrics::vector_linf_norm,
        "Compute Linf norm of vector");
    m.def("matrix_frobenius_norm", &norms_metrics::matrix_frobenius_norm,
        "Compute Frobenius norm of matrix");
    m.def("matrix_l1_norm", &norms_metrics::matrix_l1_norm,
        "Compute L1 norm of matrix");
    m.def("matrix_linf_norm", &norms_metrics::matrix_linf_norm,
        "Compute Linf norm of matrix");
    m.def("euclidean_distance", &norms_metrics::euclidean_distance,
        "Euclidean distance between vectors",
        py::arg("v1"), py::arg("v2"));
    m.def("manhattan_distance", &norms_metrics::manhattan_distance,
        "Manhattan distance between vectors",
        py::arg("v1"), py::arg("v2"));
    m.def("chebyshev_distance", &norms_metrics::chebyshev_distance,
        "Chebyshev distance between vectors",
        py::arg("v1"), py::arg("v2"));

    m.def("dft", &transforms::dft,
        "Discrete Fourier Transform",
        py::arg("signal"));
    m.def("idft", &transforms::idft,
        "Inverse Discrete Fourier Transform",
        py::arg("transform"));
    m.def("fft", &transforms::fft,
        "Fast Fourier Transform",
        py::arg("signal"));
    m.def("laplace_transform", &transforms::laplace_transform,
        "Laplace Transform",
        py::arg("f"), py::arg("s_real"), py::arg("s_imag"),
        py::arg("t_max"), py::arg("dt") = 0.01);

    m.def("plot_solution", &visualization::plot_solution,
        "Plot ODE solution",
        py::arg("solution"));
}