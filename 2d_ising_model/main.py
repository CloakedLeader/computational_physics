import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import random
import math

L = 50
class IsingLattice:
    def __init__(self):
        self.lattice = self.create_base_lattice()
        

    def create_base_lattice(self):
        a = np.random.randint(0, 2, (L,L))
        a = (a*2) - 1
        return a
   
    def find_neighbours(self, position):
        i, j = position
        north = (i, (j-1)%L)
        south = (i, (j+1)%L)
        west = ((i-1)%L, j)
        east = ((i+1)%L, j)
        return north, east, south, west

    def random_spin(self):
        pos = (random.randint(0,L-1), random.randint(0,L-1))
        return pos
    
    def compute_energy_change(self, pos):
        neigbour_locations = self.find_neighbours(pos)
        neigbour_values_sum = 0
        for coord in neigbour_locations:
            neigbour_values_sum += self.lattice[coord]
        deltaE = 2*self.lattice[self.random_spin()]*neigbour_values_sum
        return deltaE
    
    def flip_spin(self, pos):
        self.lattice[pos] = self.lattice[pos]*-1
    
    def find_boltzman_prob(self, energy_change, temp):
        power = -1*energy_change/(temp)
        P = math.exp(power)
        return P


    def run(self, sweeps=1000, T=2.0):
        plt.ion()
        fig, ax = plt.subplots()
        im = ax.imshow(self.lattice, cmap='coolwarm', interpolation='nearest')
        ax.set_title("Ising Model: Sweep 0")
        ax.axis('off')
        for sweep in range(sweeps):
            for _ in range(L**2):
                i = 0
                while i < 500:
                    i += 1
                    position = self.random_spin()
                    energy_change = self.compute_energy_change(position)
                    if energy_change < 0:
                        self.flip_spin(position)
                    else:
                        probability = self.find_boltzman_prob(energy_change, T)
                        rand_float = random.random()
                        if rand_float <= probability:
                            self.flip_spin(position)
            im.set_data(self.lattice)
            ax.set_title(f"Ising Model: Sweep {sweep + 1}")
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.1)
        plt.ioff()
        plt.show()

tester = IsingLattice()
tester.run()
