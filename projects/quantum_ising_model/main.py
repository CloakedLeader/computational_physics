import numpy as np
from matplotlib import pyplot as plt


class Hamiltonian:

    def __init__(self, num_particles, h, J):
        self.n = num_particles
        self.h = h
        self.J = J
        self.sz = np.array([[1, 0], [0, -1]])
        self.sx = np.array([[0, 1], [1, 0]])
        self.I = np.eye(2)

    @staticmethod
    def tensor_product(*args):
        result = np.array([[1]])
        for op in args:
            result = np.kron(result, op)
        return result

    def flipping_term(self, target):
        product = [self.I] * self.n
        product[target] = self.sx
        return self.tensor_product(*product)

    def interaction_term(self, target1, target2):
        product = [self.I] * self.n
        product[target1] = self.sz
        product[target2] = self.sz
        return self.tensor_product(*product)

    def build(self):
        dim = 2 ** self.n
        H = np.zeros((dim, dim))

        for i in range(self.n):
            H -= self.h * self.flipping_term(i)

        for i in range(self.n - 1):
            H -= self.J * self.interaction_term(i, (i + 1) % self.n)

        return H


def find_eigenstates(n, h, J):
    H = Hamiltonian(num_particles=n, h=h, J=J).build()
    eigenvals, eigenvecs = np.linalg.eigh(H)
    return eigenvals


Js = np.linspace(0.1, 2.0, 20)
fig, ax = plt.subplots()
for i in Js:
    vals = find_eigenstates(6, 0.1, i)
    ax.plot(vals, label=f"J = {i}")

plt.legend()
plt.show()
