#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>

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

    std::vector<std::vector<double,double>> current_lattice;
    std::vector<std::vector<double,double>> next_lattice;

public:
    ForestLattice(int size_, double f, double p)
        : size(size_), fire_prob(f), grow_prob(p)

    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                current_lattice[i][j] = 0;
            }
        }
    }

    void advanceOneStep() {
        
    }




};