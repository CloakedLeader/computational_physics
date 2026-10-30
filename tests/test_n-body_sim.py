import numpy as np
from projects.n_body_simulation.main import Bodies, Simulation
from matplotlib import pyplot as plt
import math

def test_energy_conserved_circle():
    sun   = Bodies({'x': 0, 'y': 0, 'vx': 0, 'vy': 0, 'mass': 1.0, 'name': 'Sun'})
    earth = Bodies({'x': 1, 'y': 0, 'vx': 0, 'vy': 2*np.pi, 'mass': 3e-6, 'name': 'Earth'})
    sim = Simulation([sun, earth])
    sim.run(dt=0.0001, steps=1000)
    E = np.array(sim.total_energy_hist)
    error = np.max(np.abs(E - E[0])) / abs(E[0])
    print(f"Error: {error}")
    assert error < 1e-8

def test_energy_conservation_circle():
    sun   = Bodies({'x': 0, 'y': 0, 'vx': 0, 'vy': 0, 'mass': 1.0, 'name': 'Sun'})
    earth = Bodies({'x': 1, 'y': 0, 'vx': 0, 'vy': 2*np.pi, 'mass': 3e-6, 'name': 'Earth'})
    sim = Simulation([sun, earth])
    step_sizes = np.array([0.004, 0.002, 0.001, 0.0005])
    max_errors = []
    for size in step_sizes:
        T = 20
        sim.run(dt=size, steps=round(T / size))
        E = np.array(sim.total_energy_hist)
        error = np.max(np.abs(E - E[0])) / abs(E[0])
        max_errors.append(error)
    
    max_errors = np.array(max_errors)

    for i in range(len(max_errors)):
        assert max_errors[i] < 1e-6, \
            f"Energy drift {max_errors[i]} for step size {step_sizes[i]} is too large."
    
    
def test_energy_convergence_ellipse():
    sun   = Bodies({'x': 0, 'y': 0, 'vx': 0, 'vy': 0, 'mass': 1.0, 'name': 'Sun'})
    earth = Bodies({'x': 1, 'y': 0, 'vx': 0, 'vy': 0.8*2*np.pi, 'mass': 3e-6, 'name': 'Earth'})
    sim = Simulation([sun, earth])
    step_sizes = np.array([0.008, 0.004, 0.002, 0.001, 0.0005])
    max_errors = []
    for size in step_sizes:
        T = 20
        sim.run(dt=size, steps=round(T / size))
        E = np.array(sim.total_energy_hist)
        error = np.max(np.abs(E - E[0])) / abs(E[0])
        max_errors.append(error)

    max_errors = np.array(max_errors)

    gradient, _ = np.polyfit(np.log10(step_sizes), np.log10(max_errors), 1)
    print(gradient)
    plt.plot(np.log10(step_sizes), np.log10(max_errors))
    plt.show()

    assert math.isclose(gradient, 2.0, rel_tol=0.1), \
        f"Gradient {gradient} is not close to {2.0}"
