def find_extrema(f, a, b, step=1e-5):
    """Find extrema of a function using numerical differentiation."""
    prev, curr = a, a + step
    while curr <= b:
        if (f(curr) - f(prev)) * (f(curr + step) - f(curr)) < 0:
            return curr
        prev, curr = curr, curr + step
    return None


def integrate(f, a, b, n=1000):
    """Compute definite integral using the trapezoidal rule."""
    h = (b - a) / n
    result = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        result += f(a + i * h)
    return result * h


def integrate_double(f, ax, bx, ay, by, nx=100, ny=100):
    """Compute double integral using the rectangular method."""
    hx = (bx - ax) / nx
    hy = (by - ay) / ny
    result = 0
    for i in range(nx):
        for j in range(ny):
            x = ax + i * hx
            y = ay + j * hy
            result += f(x, y) * hx * hy
    return result


def integrate_triple(f, ax, bx, ay, by, az, bz, nx=20, ny=20, nz=20):
    """Compute triple integral using the rectangular method."""
    hx = (bx - ax) / nx
    hy = (by - ay) / ny
    hz = (bz - az) / nz
    result = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x = ax + i * hx
                y = ay + j * hy
                z = az + k * hz
                result += f(x, y, z) * hx * hy * hz
    return result


def solve_ode(f, y0, t0, t1, h=0.01):
    """Solve first-order ODE using Euler's method."""
    t, y = t0, y0
    solution = [(t, y)]
    while t < t1:
        y += h * f(t, y)
        t += h
        solution.append((t, y))
    return solution
