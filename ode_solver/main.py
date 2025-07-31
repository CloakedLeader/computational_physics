"""
Solve dx/dt = -x
"""
import numpy as np
from matplotlib import pyplot as plt
import math


def f(t, x):
    return -1 * math.sin(t)


t0 = 0
x0 = 0
t_end = 3
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
    t = t0

    for _ in range(N):
        x = x + h * f(t, x)
        t += h
        x_values.append(x)
        t_values.append(t)

    return t_values, x_values

# ===============================
# Improved Euler/ Heun's Method
# ===============================


"""
Find the slope at the beginning and end of the step. The value of y at
the next point is the value at the previous point plus the average of
the two slopes times the step height.
"""


def heuns(f, t0, x0, h, N):
    t_vals = [t0]
    x_vals = [x0]
    t = t0
    x = x0

    for _ in range(N):
        k1 = f(t, x)
        t = t + h
        k2 = f(t, x + h * k1)
        x = x + h / 2 * (k1 + k2)
        t_vals.append(t)
        x_vals.append(x)

    return t_vals, x_vals


# =======================
# Runge-Kutta 4th Order
# =======================

"""
Involves finding 4 slopes, beginning, two in the middle and the end,
then averaging these to get a good approximation for the function in the
next interval.
"""


def rk4(f, t0, x0, h, N):
    t_vals = [t0]
    x_vals = [x0]
    t = t0
    x = x0
    half = h / 2

    for _ in range(N):
        k1 = f(t, x)
        k2 = f(t + half, x + half * k1)
        k3 = f(t + half, x + half * k2)
        t += h
        k4 = f(t, x + h * k3)
        x += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t_vals.append(t)
        x_vals.append(x)

    return t_vals, x_vals

# =======================
# RK4 with 2nd Order ODE
# =======================


def harmonic_oscillator(t, y, omega=2*np.pi):
    x, v = y
    dxdt = v
    dvdt = -omega**2 * x
    return np.array([dxdt, dvdt])


def second_rk4(f, t0, y0, h, N):
    t_vals = [t0]
    y_vals = [y0]
    t = t0
    y = np.array(y0)

    for _ in range(N):
        k1 = f(t, y)
        k2 = f(t + h/2, y + h/2 * k1)
        k3 = f(t + h/2, y + h/2 * k2)
        k4 = f(t + h, y + h * k3)

        y += h/6 * (k1 + 2*k2 + 2*k3 + k4)
        t += h
        t_vals.append(t)
        y_vals.append(y.copy())

    return np.array(t_vals), np.array(y_vals)


def damped_pendulum(t, y, gamma=0.1, g=9.81, L=1.0):
    theta, omega = y
    dtheta_dt = omega
    domega_dt = -gamma*omega - (g/L)*np.sin(theta)
    return np.array([dtheta_dt, domega_dt])


def duffing_oscillator(t, y, delta=0.1, alpha=-1.0, beta=1.0):
    x, v = y
    dx_dt = v
    dv_dt = -delta * v - alpha * x - beta * x**3
    return np.array([dx_dt, dv_dt])


t_vals, y_vals = second_rk4(duffing_oscillator, 0.0, [0, 2.0], 0.01, 5000)

x_vals = y_vals[:, 0]
v_vals = y_vals[:, 1]

plt.plot(t_vals, x_vals, label="x(t) - Position")
plt.plot(t_vals, v_vals, label="v(t) - Velocity")
plt.title("Duffing Oscillator (RK4)")
plt.xlabel("Time")
plt.ylabel("Value")
plt.legend()
plt.grid(True)
plt.show()
