import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

def random_num() -> float:
    return np.random.random()


class PercolationLattice:
    def __init__(self, size: int, probability: float) -> None:
        self.lattice = np.zeros((size, size))
        self.labels = np.zeros((size, size), dtype=object)
        self.prob = probability
        self.size = size
        self.all_labels = None
        self.num_clusters = None
        for i in range(size):
            for j in range(size):
                if random_num() < probability:
                    self.lattice[i][j] = 1
                else:
                    continue

    def __repr__(self) -> str:
        first_row = np.array2string(self.lattice[0], separator=", ")
        return "First row: " + first_row
    
    def label_clusters(self) -> None:
        label = 1
        for i in range(self.size):
            for j in range(self.size):
                if self.lattice[i][j] == 1 and self.labels[i][j] == 0:
                    self.identify_clusters((i, j), label)
                    label += 1

        self.num_clusters = label

    def check_neighbours(self, position: tuple[int, int]) -> list[tuple[int, int]]:
        positives = []
        x, y = position
        for nx, ny in [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]:
            if 0 <= nx < self.size and 0 <= ny < self.size:
                if self.lattice[nx][ny] == 1:
                    positives.append((nx, ny))
        return positives

    def identify_clusters(self, position: tuple[int, int], current_label) -> None:
        x, y = position
        if self.labels[x][y] != 0:
            return
        self.labels[x][y] = current_label
        for nx, ny in self.check_neighbours(position):
            if self.labels[nx][ny] == 0:
                self.identify_clusters((nx, ny), current_label)

    def is_spanning_cluster(self) -> bool:
        updownlabels = []
        for i in range(self.size):
            if self.labels[0][i] != 0:
                label = self.labels[0][i]
                updownlabels.append(label)
        if len(updownlabels) != 0:
            for i in range(self.size):
                if self.labels[self.size -1][i] in updownlabels:
                    return True
        
        leftrightlabels = []
        for i in range(self.size):
            if self.labels[i][0] != 0:
                label = self.labels[i][0]
                leftrightlabels.append(label)
        if len(leftrightlabels) != 0:
            for i in range(self.size):
                if self.labels[i][self.size - 1] in leftrightlabels:
                    return True
        
        return False
    
    def run(self) -> bool:
        self.label_clusters()
        return self.is_spanning_cluster()


def average_over_p(num_of_p: int, lattice_size: int, trials_per_p: int) -> list[tuple[float, float]]:
    ps = np.linspace(start=0.05, stop=0.95, num=num_of_p)
    results: list[tuple[float, float]] = []
    for prob in ps:
        spannings = 0
        for _ in range(trials_per_p):
            percol = PercolationLattice(lattice_size, prob)
            spannings += 1 if percol.run() else 0
        results.append((prob, spannings/trials_per_p))
    
    return results

res = average_over_p(20, 30, 100)
ps, spanning_probs = zip(*res)
f = interp1d(spanning_probs, ps)
p_c = float(f(0.5))
print("Estimated critical probability: ", p_c)
plt.plot(ps, spanning_probs, 'o-')
plt.xlabel('p')
plt.ylabel('Spanning Probability')
plt.show()

