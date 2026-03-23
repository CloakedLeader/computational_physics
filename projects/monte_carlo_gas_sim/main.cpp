#define _USE_MATHS_DEFINES

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


class GasLattice{
private:
    int number;
    int steps;
    float size;
    float delta = 0.01;
    
    std::vector<std::vector<float>> current_positions;
    std::vector<std::vector<float>> next_positions;


public:
    GasLattice(int number_, int steps_, float size_) : number(number_), steps(steps_), size(size_)
    {
        current_positions.resize(number, std::vector<float>(3, 0));
        next_positions.resize(number, std::vector<float>(3, 0));
    }

    std::vector<float> generate_small_vector() {
        double phi = random_num() * 2 * M_PI;
        double theta = random_num() * M_PI;
        std::vector<float> unit_vec = {
            static_cast<float>(sin(theta)*cos(phi)),
            static_cast<float>(sin(phi)*sin(theta)),
            static_cast<float>(cos(theta))};
        for (int i = 0; i < 3; i++) {
            unit_vec[i] *= delta;
        };
        return unit_vec;
    }

    void monte_carlo_step() {
        int particle_num = static_cast<int>(random_num()*number);
        std::vector<float> particle_pos = current_positions[particle_num];
        std::vector<float> tiny_displacement = generate_small_vector();
        std::vector<float> new_pos;
        new_pos.resize(3);
        for (int i = 0; i < 3; i++) {
            float dist = particle_pos[i] + tiny_displacement[i];
            if (dist > size) {
                new_pos[i] = dist - size;
            }
            else if (dist < 0) {
                new_pos[i] = size + dist;
            }
            else {
                new_pos[i] = dist;
            };
        };
        for (int i = 0; i < number; i++) {
            if (current_positions[i] == new_pos) {
                next_positions[particle_num] = new_pos;
                break;
                /* Cannot use equality since the floating points will never be directly equal.
                Must use a epsilon i.e a radius in which we can consider the particles to "touching".
                Implement that here and then steps should be valid.
                Ideally find a way to reduce check time by not considering all other positions but only those
                that are within a quite large range, otherwise the code will slow down signinficantly. Especially
                when many bodies are being simulated.*/
            };
        };
    }

    void monte_carlo_sweep() {
        for (int i = 0; i < number; i++) {
            monte_carlo_step();
        }
    }

};