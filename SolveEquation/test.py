import solve_equation
from scipy.optimize import fsolve
import math

eq = solve_equation.SolveEquation()

# fun = lambda x: (x + 2 / x) / 2
# fun2 = lambda x: x**3 + x - 1
# dirf = lambda x: 3 * x**2 + 1

# fun3 = lambda x: x**2
# dirf3 = lambda x: 2 * x
# xc = eq.newton(fun3, dirf3, 1, 10, m=1)
# print(xc)
# print(fsolve(fun3, 0))

# xc = eq.fpi(fun, 1, 8)
# print("FPI:", xc)
# xc = eq.bisect(fun2, -1, 0)
# print("Binsect:", xc)
# xc = fsolve(fun2, 0)
# print("Fsolve:", xc[0])

fun = lambda x: math.sin(x) + x**2 * math.cos(x) - x**2 - x
dirf = lambda x: math.cos(x) + 2 * x * math.cos(x) - x**2 * math.sin(x) - 2 * x - 1

fun1 = lambda x: x**3 + x - 1

xc = eq.newton(fun, dirf, 1, 5, m=3)
print(xc)

xc = eq.secant(fun1, 0, 1, 9)
print(xc)