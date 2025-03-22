def determinant(matrix):
    """Compute the determinant of a square matrix using recursion."""
    if len(matrix) == 1:
        return matrix[0][0]

    det = 0
    for col in range(len(matrix)):
        # Compute the minor matrix by removing the first row and current column
        minor = [row[:col] + row[col + 1:] for row in matrix[1:]]
        # Apply Laplace expansion
        det += (-1) ** col * matrix[0][col] * determinant(minor)

    return det


def hermitian_conjugate(matrix):
    """Compute the Hermitian conjugate (conjugate transpose) of a complex matrix."""
    return [[(x[0], -x[1]) for x in col] for col in zip(*matrix)]


def matrix_multiply(A, B):
    """Multiply two matrices A and B."""
    return [
        [sum(a * b for a, b in zip(A_row, B_col)) for B_col in zip(*B)]
        for A_row in A
    ]


def matrix_exponential(matrix, iterations=10):
    """Compute the matrix exponential using the power series expansion."""
    n = len(matrix)

    # Initialize the result as an identity matrix
    result = [[(1.0 if i == j else 0.0, 0.0) for j in range(n)] for i in range(n)]
    current_term = [[(x[0], x[1]) for x in row] for row in result]

    for k in range(1, iterations + 1):
        # Compute the next term in the power series
        current_term = matrix_multiply(current_term, matrix)
        current_term = [[(x[0] / k, x[1] / k) for x in row] for row in current_term]

        # Add the term to the result
        for i in range(n):
            for j in range(n):
                result[i][j] = (result[i][j][0] + current_term[i][j][0],
                                result[i][j][1] + current_term[i][j][1])

    return result


def matrix_power(matrix, power):
    """Raise a square matrix to an integer power using exponentiation by squaring."""
    if not matrix:
        return []

    n = len(matrix)

    # Ensure the matrix is square
    if any(len(row) != n for row in matrix):
        raise ValueError("Matrix must be square")

    # Initialize the identity matrix
    result = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    current = [row.copy() for row in matrix]
    exponent = abs(power)

    while exponent > 0:
        if exponent % 2 == 1:
            result = matrix_multiply(result, current)
        current = matrix_multiply(current, current)
        exponent //= 2

    # Handle negative exponent by computing the inverse
    if power < 0:
        result = matrix_inverse(result)

    return result


def matrix_inverse(matrix):
    """Compute the inverse of a square matrix using the Gauss-Jordan elimination method."""
    n = len(matrix)

    # Augment the matrix with the identity matrix
    augmented = [row + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(matrix)]

    for i in range(n):
        # Normalize the pivot row
        pivot = augmented[i][i]
        if pivot == 0:
            raise ValueError("Matrix is singular and cannot be inverted")

        for j in range(2 * n):
            augmented[i][j] /= pivot

        # Eliminate other entries in the current column
        for k in range(n):
            if k != i and augmented[k][i] != 0:
                factor = augmented[k][i]
                for j in range(2 * n):
                    augmented[k][j] -= factor * augmented[i][j]

    # Extract the inverse matrix
    return [row[n:] for row in augmented]
