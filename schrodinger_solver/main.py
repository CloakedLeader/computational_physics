import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd


"""
Finite Difference Code
"""

N = 500
T = 1000

x_min, x_max = -10, 10
x = np.linspace(x_min, x_max, N)
dx = x[1] - x[0]
dt = 0.1 * dx**2
print(dt)

x0 = -5
sigma = 1
k0 = 5
psi0 = np.exp(-(x - x0)**2 / (2 * sigma**2)) * np.exp(1j * k0 * x)

psi0 /= np.sqrt(np.sum(np.abs(psi0)**2) * dx)

psi_current = psi0.copy()
psi_next = np.zeros_like(psi0)

V = np.zeros(N)
V[220:300] = 10

psi_all = np.zeros((T, N), dtype=complex)
psi_all[0, :] = psi_current
for t in range(1, T):
    for j in range(1, N-1):
        laplacian = (psi_current[j+1] - 2 * psi_current[j] + psi_current[j-1]) / (2 * dx**2)
        psi_next[j] = psi_current[j] + 1j * dt * (laplacian - V[j] * psi_current[j])
    
    psi_next[0] = psi_next[-1] = 0

    psi_current, psi_next = psi_next, psi_current
    psi_all[t, :] = psi_current

plt.plot(x, np.abs(psi_all[0])**2, label="t=0")
plt.plot(x, np.abs(psi_all[T//4])**2, label="t=1/4 T")
plt.plot(x, np.abs(psi_all[T//2])**2, label="t=1/2 T")
plt.plot(x, np.abs(psi_all[-1])**2, label="t=T")
plt.legend()

plt.show()