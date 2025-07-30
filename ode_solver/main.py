"""
Solve dx/dt = -x
"""
import numpy as np
from matplotlib import pyplot as plt


def f(t, x):
    return -x


t0 = 0
x0 = 1
t_end = 5
h = 0.05
N = int((t_end - t0) / h)

# ===============
# Euler Method
# ===============

"""
To find the value at the next step you do the value at the previous step
plus the gradient times the step value at the previous step.
"""


def euler(f, t0, x0, h, N):
    t_values = [t0]
    x_values = [x0]
    x = x0
    t=t0

    for _ in range(N):
        x = x + h * f(t, x)
        t += h
        x_values.append(x)
        t_values.append(t)

    return t_values, x_values

t_exact = np.linspace(t0, t_end, 1000)
x_exact = np.exp(-t_exact)

t_euler, x_euler = euler(f, t0, x0, h, N)

plt.plot(t_exact, x_exact, 'k--', label='Exact')
plt.plot(t_euler, x_euler, label='Euler')
plt.legend()
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("ODE: dx/dt = -x")
plt.grid(True)
plt.show()
