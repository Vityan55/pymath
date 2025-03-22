def monte_carlo_integral(f, a, b, n=1000):
    """Estimate integral using Monte Carlo method"""
    total = 0
    seed = 1
    for _ in range(n):
        seed = (seed * 1103515245 + 12345) % (2 ** 31)
        x = a + (b - a) * (seed / (2 ** 31))
        total += f(x)
    return (b - a) * total / n
