import numpy as np
from matplotlib import pyplot as plt


def intial_cond(x):
    return np.exp(-100 * (x - 0.5)**2)


def solve_heat_eq(L=1.0, alpha=0.05, T=1.0, N=20, conv=0.5):
    delta_x = L / N
    delta_t = conv * (delta_x**2) / alpha
    num_time_points = max(2, int(T / delta_t))

    grid = np.zeros((N, num_time_points))
    x = np.linspace(0, L, N)
    grid[:, 0] = np.sin(np.pi * intial_cond(x))

    grid[0, :] = 0
    grid[1, :] = 0

    for n in range(0, num_time_points - 1):
        for i in range(1, N - 1):
            grid[i, n+1] = grid[i, n] + conv * (grid[i-1, n] - 2*grid[i, n] + grid[i+1, n])

    return x, grid, delta_t


x, grid, dt = solve_heat_eq()
time_indices = [0, int(0.1/dt), int(0.2/dt), int(0.3/dt), int(0.5/dt)]

plt.figure(figsize=(10, 6))
# for t_idx in time_indices:
#     print(grid[:, t_idx])
#     plt.plot(x, grid[:, t_idx], label=f"t = {round(t_idx * dt, 3)}s")

plt.plot(x, grid[:, 0], label="Intial")
plt.plot(x, grid[:, -1], label="Final")
plt.title("Heat Diffusion Over Time")
plt.xlabel("Position along the rod")
plt.ylabel("Temperature")
plt.legend()
plt.grid(True)
plt.show()
