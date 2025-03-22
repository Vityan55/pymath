from .math_utils import PI, taylor_sin, taylor_cos, taylor_exp  # Added a dot before the module name


def dft(signal):
    """
    Computes the Discrete Fourier Transform (DFT) without external dependencies.

    :param signal: List of real or complex numbers representing the input signal.
    :return: List of tuples (real, imaginary) representing the DFT result.
    """
    N = len(signal)
    transform = []
    for k in range(N):
        real = 0.0
        imag = 0.0
        for n in range(N):
            angle = -2 * PI * k * n / N
            c = taylor_cos(angle)
            s = taylor_sin(angle)
            real += signal[n] * c
            imag += signal[n] * s
        transform.append((real, imag))
    return transform


def idft(transform):
    """
    Computes the Inverse Discrete Fourier Transform (IDFT).

    :param transform: List of tuples (real, imaginary) representing the frequency domain signal.
    :return: List of real numbers representing the reconstructed time-domain signal.
    """
    N = len(transform)
    signal = []
    for n in range(N):
        real = 0.0
        imag = 0.0  # Initialize imaginary part
        for k in range(N):
            angle = 2 * PI * k * n / N
            c = taylor_cos(angle)
            s = taylor_sin(angle)
            # Correct calculation of both components
            real += transform[k][0] * c - transform[k][1] * s
            imag += transform[k][0] * s + transform[k][1] * c
        # Normalize and return only the real part
        signal.append(real / N)
    return signal


def fft(signal):
    """
    Recursive implementation of the Fast Fourier Transform (Cooley-Tukey algorithm).
    Requires the signal length to be a power of 2.

    :param signal: List of real or complex numbers representing the input signal.
    :return: List of tuples (real, imaginary) representing the FFT result.
    """
    N = len(signal)

    # Base case for recursion
    if N <= 1:
        return signal

    # Ensure the length is a power of 2
    if (N & (N - 1)) != 0:
        raise ValueError("Signal length must be a power of 2")

    # Split into even and odd elements
    even = fft(signal[::2])
    odd = fft(signal[1::2])

    # Combine results
    result = [0] * N
    for k in range(N // 2):
        twiddle_factor = (-2 * PI * k / N)
        c = taylor_cos(twiddle_factor)
        s = taylor_sin(twiddle_factor)

        odd_real = odd[k][0] * c - odd[k][1] * s
        odd_imag = odd[k][0] * s + odd[k][1] * c

        result[k] = (
            even[k][0] + odd_real,
            even[k][1] + odd_imag
        )

        result[k + N // 2] = (
            even[k][0] - odd_real,
            even[k][1] - odd_imag
        )

    return result


def laplace_transform(f, s_real=1.0, s_imag=0.0, t_max=10.0, dt=0.01):
    """
    Numerical computation of the Laplace Transform.

    :param f: Function of time t -> float.
    :param s_real: Real part of the complex frequency s.
    :param s_imag: Imaginary part of the complex frequency s.
    :param t_max: Upper limit of integration.
    :param dt: Time step for numerical integration.
    :return: Tuple (real, imaginary) representing the Laplace transform result.
    """
    integral_real = 0.0
    integral_imag = 0.0
    t = 0.0

    while t <= t_max:
        # Compute the exponent term e^(-st) = e^(-σt) * [cos(ωt) - j sin(ωt)]
        st_real = s_real * t
        st_imag = s_imag * t
        exp_real = taylor_exp(-st_real) * taylor_cos(st_imag)
        exp_imag = -taylor_exp(-st_real) * taylor_sin(st_imag)

        # Multiply by f(t)
        ft = f(t)
        product_real = ft * exp_real
        product_imag = ft * exp_imag

        # Numerical integration
        integral_real += product_real * dt
        integral_imag += product_imag * dt

        t += dt

    return integral_real, integral_imag
