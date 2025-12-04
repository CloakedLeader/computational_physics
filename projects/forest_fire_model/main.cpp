#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <array>

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

    std::vector<std::vector<int>> current_lattice;
    std::vector<std::vector<int>> next_lattice;

    struct {
        int newIgnitions = 0;
        int totalBurned = 0;
    } Metrics;


public:
    ForestLattice(int size_, int steps_, double f, double p)
        : size(size_), steps(steps_), fire_prob(f), grow_prob(p)

    {
        current_lattice.resize(size, std::vector<int>(size, 0));
        next_lattice.resize(size, std::vector<int>(size, 0));
    }

    void lightNeighbourTrees(int x, int y) {
        const int xOffsets[4] = {-1, 1, 0, 0};
        const int yOffsets[4] = {0, 0, -1, 1};
        for (int k = 0; k < 4; k++) {
            int nx = x + xOffsets[k];
            int ny = y + yOffsets[k];
            if (current_lattice[nx][ny] == 1 &&
                0 <= nx < size && 0 <= ny < size) {
                    next_lattice[nx][ny] = 2;
                }
        }

    }

    void advanceOneStep() {

        next_lattice = current_lattice;
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
        file << "#step_number, empty, tree, burning";
        for (int t = 0; t < steps; t++) {
            advanceOneStep();
            std::array<int, 3> inst = countInstMetrics();
            file << t << "," << inst[0] << "," << inst[1] << "," << inst[2] << '\n';
        }
        file.close();
        std::cout << "Data written to results.csv" << '\n';
    }
};

void main () {
    ForestLattice forest(100, 500, 0.005, 0.2);
    forest.run();
}