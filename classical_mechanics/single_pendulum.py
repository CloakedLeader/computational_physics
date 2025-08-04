from main_functions import generate_paths
from matplotlib import pyplot as plt
import numpy as np
from pde_solver.ode_solver import second_rk4


N = 1000
t = np.linspace(0, 1, N)
dt = t[1] - t[0]


def compute_pendulum_action(theta, dt, m=1.0, g=9.81, L=1.0):
    dtheta_dt = np.diff(theta) / dt
    kinetic = 0.5 * L**2 * m * dtheta_dt**2
    potential = m * g * L * (1 - np.cos(theta[:-1]))
    L = kinetic - potential
    S = np.sum(L) * dt
    return S


paths = []
actions = []

for _ in range(500):
    theta_path = generate_paths(np.pi/2, 0, t, 10, 0.2)
    S = compute_pendulum_action(theta_path, dt)
    paths.append(theta_path)
    actions.append(S)

min_idx = np.argmin(actions)
best_path = paths[min_idx]


def pendulum_ode(t, y, g=9.81, L=1.0):
    theta, omega = y
    dtheta_dt = omega
    domega_dt = -g*L*np.sin(theta)
    return np.array([dtheta_dt, domega_dt])


t_vals, y_vals = second_rk4(pendulum_ode, 0, [np.pi/2, 0.0], 0.01, 1000)
theta_vals = y_vals[:, 0]
omega_vals = y_vals[:, 1]

for i in range(10):
    plt.plot(t, paths[i], alpha=0.4, label=f"S={actions[i]:.3f}")

plt.plot(t, best_path, 'k', linewidth=2,
         label=f'Min Action={actions[min_idx]:.3f}')
plt.plot(t_vals, theta_vals, label="Solution to ODE")
plt.legend()
plt.xlabel('Time')
plt.ylabel('x(t)')
plt.title('Sample Paths Between (0,0) and (1,1)')
plt.grid(True)
plt.show()
