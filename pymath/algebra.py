def solve_slae(A, b):
    """Solve a system of linear equations using Gaussian elimination with partial pivoting."""
    n = len(A)
    for i in range(n):
        # Partial pivoting: find the maximum element in the current column
        max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
        if A[max_row][i] == 0:
            raise ValueError("Singular matrix: system has no unique solution")

        # Swap rows if needed
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
            b[i], b[max_row] = b[max_row], b[i]

        # Gaussian elimination
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]

    # Back substitution
    x = [0] * n
    for i in range(n - 1, -1, -1):
        if A[i][i] == 0:
            raise ValueError("Singular matrix detected during back substitution")
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    return x


