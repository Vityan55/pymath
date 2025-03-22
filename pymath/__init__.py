"""
Pymath: A mathematical package for various computations.
"""

from .calculus import find_extrema, integrate, integrate_double, integrate_triple, solve_ode
from .algebra import solve_slae
from .interpolation import (lagrange_interpolation, newton_interpolation, spline_interpolation,
                            lagrange_approximation, newton_approximation, spline_approximation)
from .matrix import determinant, matrix_multiply, matrix_exponential, hermitian_conjugate, matrix_power, matrix_inverse
from .montecarlo import monte_carlo_integral
from .visualization import plot_solution
from .math_utils import newton_sqrt, taylor_exp, taylor_sin, taylor_cos, complex_exp
from .norms_metrics import (vector_l1_norm, vector_linf_norm, vector_l2_norm, matrix_linf_norm, matrix_frobenius_norm,
                            matrix_l1_norm, manhattan_distance, euclidean_distance, chebyshev_distance)
from .transforms import dft, idft, fft, laplace_transform
