def lagrange_interpolation(x_vals, y_vals, x):
    """Perform Lagrange interpolation."""
    result = 0
    n = len(x_vals)
    for i in range(n):
        term = y_vals[i]
        for j in range(n):
            if i != j:
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
        result += term
    return result


def newton_interpolation(x_vals, y_vals, x):
    """Perform Newton interpolation."""
    n = len(x_vals)
    coef = [y for y in y_vals]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (x_vals[i] - x_vals[i - j])
    result = coef[-1]
    for i in range(n - 2, -1, -1):
        result = result * (x - x_vals[i]) + coef[i]
    return result


def spline_interpolation(x, y, xi):
    """Perform cubic spline interpolation using natural boundary conditions."""
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]

    # Step 1: Compute the coefficients alpha
    alpha = [0] * n
    for i in range(1, n - 1):
        alpha[i] = (1 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1])

    # Step 2: Solve the tridiagonal system using forward elimination
    l = [1] + [0] * (n - 1)
    mu = [0] * n
    z = [0] * n

    for i in range(1, n - 1):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    # Step 3: Compute coefficients c, b, d
    c = [0] * n
    b = [0] * (n - 1)
    d = [0] * (n - 1)

    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    # Step 4: Find the correct spline segment
    i = max([j for j in range(n - 1) if x[j] <= xi])  # Binary search for efficiency
    dx = xi - x[i]
    return y[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3


def lagrange_approximation(x_vals, y_vals):
    """Perform Lagrange polynomial approximation."""
    return lambda x: lagrange_interpolation(x_vals, y_vals, x)


def newton_approximation(x_vals, y_vals):
    """Perform Newton polynomial approximation."""
    return lambda x: newton_interpolation(x_vals, y_vals, x)


def spline_approximation(x_vals, y_vals):
    """Perform cubic spline approximation."""
    return lambda x: spline_interpolation(x_vals, y_vals, x)
