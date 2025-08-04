import sympy as sp
import numpy as np


x, y = sp.symbols("x,y")

f = x**3 * y + sp.sin(x * y)

dfdx = sp.diff(f, x)
dfdy = sp.diff(f, y)
dfdxdy = sp.diff(dfdx, y)
dfdydx = sp.diff(dfdy, x)

print(dfdx)
print(dfdy)
print(dfdxdy)
print(dfdydx)

print(True if sp.simplify(dfdxdy - dfdydx) == 0 else False)