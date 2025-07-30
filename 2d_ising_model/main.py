import numpy as np
import matplotlib.pyplot as plt
import random
import math


L = 25


class IsingLattice:
    def __init__(self):
        self.lattice = self.create_base_lattice()

    def create_base_lattice(self):
        a = np.random.randint(0, 2, (L, L))
        a = (a*2) - 1
        return a

    def find_neighbours(self, position):
        i, j = position
        north = (i, (j-1) % L)
        south = (i, (j+1) % L)
        west = ((i-1) % L, j)
        east = ((i+1) % L, j)
        return north, east, south, west

    def random_pos(self):
        pos = (random.randint(0, L-1), random.randint(0, L-1))
        return pos

    def compute_energy_change(self, pos):
        neigbour_locations = self.find_neighbours(pos)
        neigbour_values_sum = 0
        for coord in neigbour_locations:
            neigbour_values_sum += self.lattice[coord]
        deltaE = 2*self.lattice[pos]*neigbour_values_sum
        return deltaE

    def flip_spin(self, pos):
        self.lattice[pos] = self.lattice[pos]*-1

    def find_boltzman_prob(self, energy_change, temp):
        power = -1*energy_change/(temp)
        P = math.exp(power)
        return P

    def measure_magnetisation(self):
        total_magnetisation = (np.sum(self.lattice))
        return abs(total_magnetisation)

    def measure_energy(self):
        energy = 0
        for i in range(L):
            for j in range(L):
                spin_state = self.lattice[i, j]
                all_neighbours = self.find_neighbours((i, j))
                right = self.lattice[all_neighbours[1]]
                down = self.lattice[all_neighbours[2]]
                neigbours_sum = right + down
                energy += -spin_state * neigbours_sum
        return energy

    def one_sweep(self, T):
        for _ in range(L**2):
            position = self.random_pos()
            energy_change = self.compute_energy_change(position)
            if energy_change < 0:
                self.flip_spin(position)
            else:
                probability = self.find_boltzman_prob(energy_change, T)
                rand_float = random.random()
                if rand_float <= probability:
                    self.flip_spin(position)

    def find_avgs(self, mag_energ):
        mag, energ = mag_energ
        total_spins = self.lattice.size
        avg_M = sum(mag) / len(mag)
        avg_E = sum(energ) / len(energ)
        mag_per_spin = avg_M / total_spins
        energ_per_spin = avg_E / total_spins
        return mag_per_spin, energ_per_spin

    def run(self, measure_sweeps=500, equil_sweeps=500, T=2.0):
        mag = []
        energ = []
        for sweep in range(equil_sweeps):
            self.one_sweep(T=T)
        for sweep in range(measure_sweeps):
            self.one_sweep(T=T)
            energy = self.measure_energy()
            magnetisation = self.measure_magnetisation()
            mag.append(magnetisation)
            energ.append(energy)
        return mag, energ


temperatures = np.linspace(1, 4, num=25)
mag_vals = []
energ_vals = []
for i in temperatures:
    dummy = IsingLattice()
    info = dummy.run(measure_sweeps=100, equil_sweeps=100, T=i)
    mag, energ = dummy.find_avgs(info)
    mag_vals.append(mag)
    energ_vals.append(energ)

print(f"Length temperatures: {len(temperatures)}")
print(f"Length mag_vals : {len(mag_vals)}")
print(f"Length energ_vals: {len(energ_vals)}")

plt.scatter(temperatures, mag_vals, color='blue',
            marker='x', label="|M| per spin")
plt.scatter(temperatures, energ_vals, color='red',
            marker='x', label="E per spin")

plt.xlabel("Temperture (T)")
plt.ylabel("Magnetisation per spin")
plt.title("Ising Model Magnetisation and Energy vs Temperature")
plt.legend()
plt.grid(True)

plt.show()
