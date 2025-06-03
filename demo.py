from pymath import (
    matrix_inverse,
    monte_carlo_integral,
    lagrange_interpolation,
    find_extrema,
    plot_solution,
    solve_ode
)

def main():
    # Пример правой части ОДУ: dy/dt = -y
    def ode_func(t, y):
        return -y

    # Параметры задачи
    y0 = 1.0      # начальное значение y
    t0 = 0.0      # начальное время
    t1 = 5.0      # конечное время
    h = 0.2       # шаг

    # Вычисление численного решения
    solution = solve_ode(ode_func, y0, t0, t1, h)

    # Отображение результата
    plot_solution(solution)
    
    # # Demonstration of matrix_inverse
    # A = [[1.0, 2.0], [3.0, 4.0]]
    # inv = matrix_inverse(A)
    # print("Inverse matrix of [[1, 2], [3, 4]]:")
    # for row in inv:
    #     print(row)

    # # Demonstration of Monte Carlo integration
    # integral = monte_carlo_integral(lambda x: x**2, 0, 1, 10000)
    # print(f"\nIntegral of x^2 on [0, 1] ≈ {integral:.4f} (expected ~0.3333)")

    # # Lagrange interpolation
    # x_vals = [0, 1, 2]
    # y_vals = [1, 3, 2]
    # value = lagrange_interpolation(x_vals, y_vals, 1.5)
    # print(f"\nValue of the Lagrange interpolation polynomial at 1.5: {value:.4f}")

    # # Extremum finding
    # f = lambda x: -(x - 2)**2 + 4
    # x_ext = find_extrema(f, 0, 4)
    # print(f"\nExtremum of the function -(x - 2)^2 + 4 found at x = {x_ext:.4f}")

if __name__ == "__main__":
    main()

