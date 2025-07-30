"""
Between 2 fixed endpoints x(0) = 0 and x(T) = 1
find the path which minimises the action.
Use the free particle Lagrangian:
L = 1/2 * m * v^2
S = int_0^1 L dt
where L is a single variable function of t.

Paths to do by hand:
1) Straight line x(t) = t
2) Quadratic x(t) = t^2
3) Cubic x(t) = 3t^2 - 2t^3
4) Sine pertubation x(t) = t + a sin(pi t)
5) Step function
6) Exponential x(t) = (e^t - 1)/(e - 1)
7) Wiggle x(t) = t + a sin(2pi t)

For each compute:
v(t), L, S and use m = 1 to simplify.
"""

import numpy as np
from matplotlib import pyplot as plt
from main_functions import generate_paths


N = 100
t = np.linspace(0, 1, N)
dt = t[1] - t[0]
num_paths = 100
actions = []
paths = []


def compute_action(x, dt, m=1.0):
    dxdt = np.diff(x) / dt
    L = 0.5 * m * dxdt**2
    S = np.sum(L) * dt
    return S


for _ in range(num_paths):
    x_path = generate_paths(t, num_modes=5, epsilon=0.2)
    S = compute_action(x_path, dt)
    paths.append(x_path)
    actions.append(S)

min_idx = np.argmin(actions)
best_path = paths[min_idx]

for i in range(5):
    plt.plot(t, paths[i], alpha=0.4, label=f"S={actions[i]:.3f}")

plt.plot(t, best_path, 'k', linewidth=2, label='Min Action')
plt.plot(t, t, '--', label='Classical path')
plt.legend()
plt.xlabel('Time')
plt.ylabel('x(t)')
plt.title('Sample Paths Between (0,0) and (1,1)')
plt.grid(True)
plt.show()
