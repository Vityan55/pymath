from pymath import plot_solution


def test_plot_solution():
    solution = [(0, 1), (1, 2), (2, 3)]
    try:
        plot_solution(solution)
        assert True
    except:
        assert False
