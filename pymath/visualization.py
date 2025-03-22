def plot_solution(solution):
    """Visualize the solution of an ODE using ASCII plotting."""
    min_y = min(y for _, y in solution)
    max_y = max(y for _, y in solution)
    scale = 20
    for t, y in solution:
        pos = int(scale * (y - min_y) / (max_y - min_y)) if max_y > min_y else 0
        print(f"{t: .2f} |" + " " * pos + "*")
