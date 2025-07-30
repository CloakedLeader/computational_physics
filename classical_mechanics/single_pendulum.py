from main_functions import generate_paths
from matplotlib import pyplot as plt
import numpy as np


N = 100
t = np.linspace(0, 1, N)
dt = t[1] - t[0]


def compute_pendulum_action(theta, dt, m=1.0, g=9.81, l=1.0):
    dtheta_dt = np.diff(theta) / dt
    kinetic = 0.5 * l**2 * m * dtheta_dt**2
    potential = m * g * l * (1 - np.cos(theta[:-1]))
    L = kinetic - potential
    S = np.sum(L) * dt
    return S


paths = []
actions = []

for _ in range(500):
    theta_path = generate_paths(0, np.pi/2, t, 5, 0.2)
    S = compute_pendulum_action(theta_path, dt)
    paths.append(theta_path)
    actions.append(S)

min_idx = np.argmin(actions)
best_path = paths[min_idx]

for i in range(10):
    plt.plot(t, paths[i], alpha=0.4, label=f"S={actions[i]:.3f}")

plt.plot(t, best_path, 'k', linewidth=2, label='Min Action')
plt.plot(t, t, '--', label='Classical path')
plt.legend()
plt.xlabel('Time')
plt.ylabel('x(t)')
plt.title('Sample Paths Between (0,0) and (1,1)')
plt.grid(True)
plt.show()
