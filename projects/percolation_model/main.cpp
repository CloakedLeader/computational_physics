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

class PercolationLattice {
private:
    int size;
    double prob;
    int num_clusters;

    std::vector<std::vector<int>> lattice;
    std::vector<std::vector<int>> labels;

public:
    PercolationLattice(int size_, double probability)
        : size(size_), prob(probability), num_clusters(0)
    {
        lattice.resize(size, std::vector<int>(size, 0));
        labels.resize(size, std::vector<int>(size, 0));

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (random_num() < prob) {
                    lattice[i][j] = 1;
                }
            }
        }
    }

    void label_clusters() {
        int label = 1;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (lattice[i][j] == 1 && labels[i][j] == 0) {
                    identify_clusters(i, j, label);
                    label++;
                }
            }
        }
        num_clusters = label - 1;
    }

    void identify_clusters(int x0, int y0, int current_label) {
        std::queue<std::pair<int, int>> q;
        q.emplace(x0, y0);

        const int dx[4] = {1, -1, 0, 0};
        const int dy[4] = {0, 0, 1, -1};

        while (!q.empty()) {
            auto [x, y] = q.front();
            q.pop();

            if (labels[x][y] != 0) {
                continue;
            }

            labels[x][y] = current_label;
            for (int k = 0; k < 4; k++) {
                int nx = x + dx[k];
                int ny = y + dy[k];
                if (nx >= 0 && nx < size &&
                    ny >= 0 && ny < size &&
                    lattice[nx][ny] == 1 &&
                    labels[nx][ny] == 0)
                {
                    q.emplace(nx, ny);
                }
            }
        }

    }

    bool is_spanning_cluster() {
        std::vector<int> top_labels;
        for (int i = 0; i < size; i++) {
            if (labels[0][i] != 0) {
                top_labels.push_back(labels[0][i]);
            }
        }
        if (!top_labels.empty()) {
            for (int i = 0; i < size; i++) {
                int label = labels[size - 1][i];
                if (std::find(top_labels.begin(), top_labels.end(), label) != top_labels.end()) {
                    return true;
                }
            }
        }

        std::vector<int> left_labels;
        for (int i = 0; i < size; i++) {
            if (labels[i][0] != 0)
                left_labels.push_back(labels[i][0]);
        }
        if (!left_labels.empty()) {
            for (int i = 0; i < size; i++) {
                int label = labels[i][size - 1];
                if (std::find(left_labels.begin(), left_labels.end(), label) != left_labels.end())
                    return true;
            }
        }

        return false;
    }

    bool run() {
        label_clusters();
        return is_spanning_cluster();
    }

};

std::vector<std::pair<double, double>> average_over_p(int num_of_p, int lattice_size, int trials) {
    std::vector<std::pair<double, double>> results;
    std::reverse(results.begin(), results.end());

    double p_start = 0.05;
    double p_end = 0.95;

    for (int i = 0; i < num_of_p; i++) {
        double p = p_start + (p_end - p_start) * (double(i) / (num_of_p - 1));

        int spannings = 0;
        for (int t = 0; t < trials; t++) {
            PercolationLattice percol(lattice_size, p);
            if (percol.run())
                spannings++;
        }

        results.emplace_back(p, double(spannings) / trials);
    }

    return results;

}

int main() {
    
    auto results = average_over_p(30, 200, 100);

    std::ofstream file("results.csv");
    file << "# p, spanning_probability" << '\n';
    for (auto& [p, sp] : results) {
        file << p << "," << sp << '\n';
    }
    file.close();

    std::cout << "Wrote results" << '\n';
}