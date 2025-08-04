from matplotlib import pyplot as plt
import numpy as np
import sympy as sp
import scipy as scp


"""
Find the roots of the equation xcosx - sinx in the range [4,5].
"""


def f(x) -> float:
    return x**3 - 6*x**2 + 11*x - 6


x_vals = np.linspace(4, 5, 1000)
f_vals = f(x_vals)


def bisector_root_finder(function, start, end, epsilon=0.00001):
    iteration = 0

    while True:
        c = (start + end) / 2
        fc = function(c)

        if fc == 0 or abs(fc) < epsilon:
            print(f"It took {iteration} iterations to get to root within {epsilon}.")
            return c

        if fc * function(start) < 0:
            end = c
        else:
            start = c

        iteration += 1


def newton_raphson(function, start, end, epsilon=0.0001):
    diff = 1.0
    x0 = (start + end) / 2

    def df(f, x, h=1e-5):
        return (f(x + h) - f(x - h)) / (2 * h)
    while diff > epsilon:
        fxn = function(x0)
        dfxn = df(function, x0)
        if abs(dfxn) < 1e-10:
            raise ZeroDivisionError("Derivate too small -- possible divergence")
        x2 = x0 - (fxn / dfxn)
        diff = abs(x2 - x0)
        x0 = x2

    return x2


root_finders = {
    "newton": newton_raphson,
    "bisection": bisector_root_finder,
}

def find_root(f, method, start, end, N):
    if method not in root_finders:
        raise ValueError(f"Unknown method '{method}'")
    h = (end - start) / N
    x0 = start
    intervals = []
    for _ in range(N):
        x1 = x0 + h
        if f(x0) * f(x1) < 0:
            intervals.append((x0, x1))
        x0 = x1

    roots = []
    for i, j in intervals:
        roots.append(root_finders[method](f, i, j))

    if len(roots) == 1:
        return roots[0]
    else:
        return roots


print(find_root(f, "bisection", 0, 4, N=150))
print(find_root(f, "newton", 0, 4, N=150))
