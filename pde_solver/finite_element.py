import numpy as np
from matplotlib import pyplot as plt


"""
Start with 1 Dimension

Solve -u''(x) = f(x), where f(x) = 1
"""


def f(x):
    return -20*np.sin(5 * np.pi * x) + 30*np.exp(-10*x)


def finite_element_method(f, num_nodes, start, end):
    h = 1 / num_nodes
    nodes = np.linspace(start, end, num_nodes+1)
    A = np.zeros((num_nodes + 1, num_nodes + 1))
    b = np.zeros(num_nodes + 1)
    for e in range(num_nodes):
        x0 = nodes[e]
        x1 = nodes[e + 1]
        h = x1 - x0

        A_e = (1 / h) * np.array([[1, -1], [-1, 1]])
        A[e:e+2, e:e+2] += A_e

        xm = (x0 + x1) / 2
        f_mid = f(xm)
        b_local = f_mid * h * 0.5 * np.array([1, 1])
        b[e: e+2] += b_local

    A_red = A[1:-1, 1:-1]
    b_red = b[1:-1]

    U_inner = np.linalg.solve(A_red, b_red)
    U = np.zeros(num_nodes + 1)
    U[1:-1] = U_inner

    return nodes, U


nodes, U = finite_element_method(f, num_nodes=100, start=0, end=1)

plt.plot(nodes, U, "o-", label="Finite Element Method")
plt.xlabel("x")
plt.ylabel("y")
plt.title("FEM Solution of -u\'\' = 1 with Dirichlet BC's")
plt.grid(True)
plt.legend()
plt.show()
