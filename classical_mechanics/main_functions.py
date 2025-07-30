import numpy as np


def generate_paths(start, end, t, num_modes, epsilon):
    """
    Creates a random path by altering the points in a linear distribution
    by a random coefficent and sinusoidal variation.
    """
    x = start + (end - start) * t
    for n in range(1, num_modes + 1):
        amplitude = epsilon * np.random.randn()
        x += amplitude * np.sin(n * np.pi * t)

    x[0] = start
    x[-1] = end
    return x
