import numpy as np
from projects.quantum_ising_model.main import Hamiltonian


def test_is_hermitian():
    qh = Hamiltonian(num_particles=3, h=1.0, J=0.5)
    H = qh.build()
    assert np.allclose(H, H.conj().T), "Hamiltonian is not Hermitian"
