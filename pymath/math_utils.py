PI = 3.141592653589793


def newton_sqrt(x, epsilon=1e-10):
    """Compute the square root using Newton's method."""
    if x < 0:
        raise ValueError("Negative number under the square root")
    guess = x
    while abs(guess * guess - x) > epsilon:
        guess = (guess + x / guess) / 2
    return guess


def taylor_sin(x, terms=15):
    """Compute sine using the Taylor series expansion."""
    x = x % (2 * PI)
    result = 0
    power = x
    sign = 1
    for n in range(1, terms + 1):
        result += sign * power
        power *= x * x / ((2 * n) * (2 * n + 1))
        sign *= -1
    return result


def taylor_cos(x, terms=15):
    """Compute cosine using the Taylor series expansion."""
    x = x % (2 * PI)
    result = 0
    power = 1
    sign = 1
    for n in range(terms):
        result += sign * power
        power *= x * x / ((2 * n + 1) * (2 * n + 2))
        sign *= -1
    return result


def complex_exp(z):
    """Compute the exponential of a complex number (z = (a, b))."""
    a, b = z
    real = taylor_cos(b) * taylor_exp(a)
    image = taylor_sin(b) * taylor_exp(a)
    return real, image


def taylor_exp(x, terms=20):
    """Compute the exponential function using the Taylor series expansion."""
    result = 1.0
    term = 1.0
    for n in range(1, terms + 1):
        term *= x / n
        result += term
    return result
