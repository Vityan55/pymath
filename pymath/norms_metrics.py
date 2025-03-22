from .math_utils import newton_sqrt


def vector_l1_norm(v):
    """
    Compute the L1 norm (Manhattan norm) of a vector.

    Args:
        v (iterable): Input vector.

    Returns:
        float: L1 norm of the vector.
    """
    return sum(abs(x) for x in v)


def vector_l2_norm(v):
    """
    Compute the L2 norm (Euclidean norm) of a vector.

    Args:
        v (iterable): Input vector.

    Returns:
        float: L2 norm of the vector.
    """
    return newton_sqrt(sum(x ** 2 for x in v))


def vector_linf_norm(v):
    """
    Compute the L-infinity norm (maximum absolute value) of a vector.

    Args:
        v (iterable): Input vector.

    Returns:
        float: L-infinity norm of the vector.
    """
    return max(abs(x) for x in v)


def matrix_frobenius_norm(matrix):
    """
    Compute the Frobenius norm of a matrix.

    Args:
        matrix (iterable of iterables): Input matrix.

    Returns:
        float: Frobenius norm of the matrix.
    """
    return newton_sqrt(sum(sum(x ** 2 for x in row) for row in matrix))


def matrix_l1_norm(matrix):
    """
    Compute the L1 norm (maximum column sum) of a matrix.

    Args:
        matrix (iterable of iterables): Input matrix.

    Returns:
        float: L1 norm of the matrix.
    """
    return max(sum(abs(row[i]) for row in matrix) for i in range(len(matrix[0])))


def matrix_linf_norm(matrix):
    """
    Compute the L-infinity norm (maximum row sum) of a matrix.

    Args:
        matrix (iterable of iterables): Input matrix.

    Returns:
        float: L-infinity norm of the matrix.
    """
    return max(sum(abs(x) for x in row) for row in matrix)


def euclidean_distance(v1, v2):
    """
    Compute the Euclidean distance between two vectors.

    Args:
        v1 (iterable): First vector.
        v2 (iterable): Second vector.

    Returns:
        float: Euclidean distance between v1 and v2.
    """
    return newton_sqrt(sum((x - y) ** 2 for x, y in zip(v1, v2)))


def manhattan_distance(v1, v2):
    """
    Compute the Manhattan distance between two vectors.

    Args:
        v1 (iterable): First vector.
        v2 (iterable): Second vector.

    Returns:
        float: Manhattan distance between v1 and v2.
    """
    return sum(abs(x - y) for x, y in zip(v1, v2))


def chebyshev_distance(v1, v2):
    """
    Compute the Chebyshev distance between two vectors.

    Args:
        v1 (iterable): First vector.
        v2 (iterable): Second vector.

    Returns:
        float: Chebyshev distance between v1 and v2.
    """
    return max(abs(x - y) for x, y in zip(v1, v2))