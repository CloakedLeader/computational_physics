#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <array>

/* TODO: Implement an Enum instead of numbers for the classification of points on the lattice.*/

double random_num() {
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

class ForestLattice {
private:
    int size;
    double fire_prob;
    double grow_prob;
    int steps;
    int subLatticeSize;

    std::vector<std::vector<int>> current_lattice;
    std::vector<std::vector<int>> next_lattice;

    struct MetricsData{
        int newIgnitions = 0;
        int totalBurned = 0;
    };

    MetricsData Metrics;

public:
    ForestLattice(int size_, int steps_, double f, double p)
        : size(size_), steps(steps_), fire_prob(f), grow_prob(p)

    {
        current_lattice.resize(size, std::vector<int>(size, 0));
        next_lattice.resize(size, std::vector<int>(size, 0));
        subLatticeSize = static_cast<int>(size * 0.1);
    }

    void lightNeighbourTrees(int x, int y) {
        const int xOffsets[4] = {-1, 1, 0, 0};
        const int yOffsets[4] = {0, 0, -1, 1};
        for (int k = 0; k < 4; k++) {
            int nx = x + xOffsets[k];
            int ny = y + yOffsets[k];
            if (current_lattice[nx][ny] == 1 &&
                0 <= nx && nx < size && 0 <= ny && ny < size) {
                    next_lattice[nx][ny] = 2;
                }
        }

    }

    std::array<std::vector<int>, 2> chooseRandomPoints () {
        std::vector<int> xs(subLatticeSize);
        std::vector<int> ys(subLatticeSize);
        for (int i = 0; i < subLatticeSize; i++) {
            xs[i] = static_cast<int>(random_num() * (subLatticeSize - 1));
            ys[i] = static_cast<int>(random_num() * (subLatticeSize - 1));
        }

        return {xs, ys};
    }

    void advanceOneStep() {
        std::swap(current_lattice, next_lattice);
        std::array<std::vector<int>, 2> randomPoints = chooseRandomPoints();
        for (int i = 0; i < subLatticeSize; i++) {
            int coords[2] = {randomPoints[0][i], randomPoints[1][i]};
            if (current_lattice[coords[0]][coords[1]] == 0) {
                    if (random_num() < grow_prob) {
                        next_lattice[coords[0]][coords[1]] = 1;
                    }
                }
                else if (current_lattice[coords[0]][coords[1]] == 1) {
                    if (random_num() < fire_prob) {
                        Metrics.newIgnitions++;
                        next_lattice[coords[0]][coords[1]] = 2;
                    }
                }
                else if (current_lattice[coords[0]][coords[1]] == 2) {
                    lightNeighbourTrees(coords[0], coords[1]);
                    Metrics.totalBurned++;
                    next_lattice[coords[0]][coords[1]] = 0;

                }

        }
        
        /*
        for (int i = 0; i < size; i++) {
            for (int j = 0; j< size; j++) {
                if (current_lattice[i][j] == 0) {
                    if (random_num() < grow_prob) {
                        next_lattice[i][j] = 1;
                    }
                }
                else if (current_lattice[i][j] == 1) {
                    if (random_num() < fire_prob) {
                        Metrics.newIgnitions++;
                        next_lattice[i][j] = 2;
                    }
                }
                else if (current_lattice[i][j] == 2) {
                    lightNeighbourTrees(i, j);
                    Metrics.totalBurned++;
                    next_lattice[i][j] = 0;

                }
            }
        }
        */
    }

    std::array<int, 3> countInstMetrics () {

        std::array<int, 3> counts = {0, 0, 0}; 
        /*Empty, tree, burning*/ 
        for (int i = 0; i < size; i++) {
            for (int j = 0; j< size; j++) {
                if (next_lattice[i][j] == 0)
                    counts[2]++;    
                else if (next_lattice[i][j] == 1)
                    counts[0]++;
                else if (next_lattice[i][j] == 2)
                    counts[1]++;
            }
        }

        return counts;
    }


    void run() {
        std::ofstream file("results.csv");
        file << "#step_number, empty, tree, burning" << '\n';
        for (int t = 0; t < steps; t++) {
            advanceOneStep();
            std::array<int, 3> inst = countInstMetrics();
            file << t << "," << inst[0] << "," << inst[1] << "," << inst[2] << '\n';
        }
        std::cout << "Data written to results.csv" << '\n';
    }
};

int main() {
    ForestLattice forest(100, 500, 0.001, 0.2);
    forest.run();

    return 0;
}