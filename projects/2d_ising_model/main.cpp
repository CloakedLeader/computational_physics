#include <Eigen/Dense> 
#include <random>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>

struct Neighbours{
    std::array<Eigen::Array2i, 4> coords;

    Neighbours(const Eigen::Array2i& centre, int latticeSize) {
        int x = centre.x();
        int y = centre.y();

        coords[0] = {x, (y - 1 + latticeSize) % latticeSize}; // North
        coords[1] = {(x + 1) % latticeSize, y};               // East
        coords[2] = {x, (y + 1) % latticeSize};               // South
        coords[3] = {(x - 1 + latticeSize) % latticeSize, y}; // West
    }

    const Eigen::Array2i& north() const { return coords[0]; }
    const Eigen::Array2i& east()  const { return coords[1]; }
    const Eigen::Array2i& south() const { return coords[2]; }
    const Eigen::Array2i& west()  const { return coords[3]; }

    auto begin() const { return coords.begin(); }
    auto end()   const { return coords.end(); }
    constexpr std::size_t size() const { return 4; }
    const Eigen::Array2i& operator[](std::size_t i) const { return coords[i]; }

};


class SpinLattice{

    public:

    SpinLattice() : L(0) {}
    
    SpinLattice(Eigen::ArrayXXi grid) : lattice(grid) {
        L = grid.rows();
    };
    
    const Eigen::ArrayXXi getLattice() const {
        return lattice;
    }

    int& operator()(const Eigen::Array2i& site) {
        return lattice(site.x(), site.y());
    }

    const int& operator()(const Eigen::Array2i& site) const {
        return lattice(site.x(), site.y());
    }

    void flip(const Eigen::Array2i& site) {
        (*this)(site) *= -1;
    }

    unsigned int magnetisation() {
        return std::abs(lattice.sum());
    }

    double energy(double J = 1.0, double h = 0.0) const {
        int horizBulk = (lattice.leftCols(L - 1) * lattice.rightCols(L - 1)).sum();
        int horizWrap = (lattice.col(L - 1) * lattice.col(0)).sum();

        int vertBulk = (lattice.topRows(L - 1) * lattice.bottomRows(L - 1)).sum();
        int vertWrap = (lattice.row(L - 1) * lattice.row(0)).sum();

        int interactionSum = horizBulk + horizWrap + vertBulk + vertWrap;

        return -J * interactionSum - h * lattice.sum();
    }

    private:
        Eigen::ArrayXXi lattice;
        int L;
};

class IsingSim{

    public:
    IsingSim(int size, float upProb = 0.5,  unsigned int seed) : latticeSize(size), rng(seed), dist(upProb), intDist(0, size - 1), doubDist(0, 1.0) {
        Eigen::ArrayXXi grid = Eigen::ArrayXXi::NullaryExpr(
            latticeSize, latticeSize,
            [this]() { return dist(rng) ? 1 : -1; }
        );
        lattice = SpinLattice(grid);
    };
    

    Neighbours findNeighbours(const Eigen::Array2i& centre) const {
        return Neighbours(centre, latticeSize);
    }

    Eigen::Array2i randomPos() {
        Eigen::Array2i pos(intDist(rng), intDist(rng));
        return pos;
    }

    double computeEnergyChange(Eigen::Array2i pos) {
        Neighbours adjacent = findNeighbours(pos);
        int adjacentSum = 0;
        for (const auto& nbr : adjacent){
            adjacentSum += lattice(nbr);
        }
        double dE = 2 * lattice(pos);
        return dE;
    }

    void flipSpin(Eigen::Array2i pos) {
        lattice.flip(pos);
    }

    double boltzmannProb(double dE, double temp) {
        double power = -1 * dE / temp;
        return std::exp(power);
    }

    unsigned int measureMag() {
        return lattice.magnetisation();
    }

    double measureEnergy() {
        return lattice.energy();
    }

    void oneSweep(double T) {
        for (int i = 0; i < latticeSize * latticeSize; i++) {
            Eigen::Array2i pos = randomPos();
            float dE = computeEnergyChange(pos);
            if (dE < 0) {
                flipSpin(pos);
            }
            else {
                double prob = boltzmannProb(dE, T);
                double rand = doubDist(rng);
                if (rand <= prob) {
                    flipSpin(pos);
                }
            }
        }
    }

    std::tuple<Eigen::ArrayXi, Eigen::ArrayXd> run(int equilSweeps, int measureSweeps, double T) {
        int N = measureSweeps;
        Eigen::ArrayXi magnetisation(N);
        Eigen::ArrayXd energy(N);

        for (int i = 0; i < equilSweeps; i++) { 
            oneSweep(T);
        }
        for (int i = 0; i < measureSweeps; i++) {
            oneSweep(T);

            magnetisation[i] = measureMag();
            energy[i] = measureEnergy();
        }

        return {magnetisation, energy};
    }

    private:
        int latticeSize;
        SpinLattice lattice;
        std::mt19937 rng;
        std::bernoulli_distribution dist;
        std::uniform_int_distribution<int> intDist;
        std::uniform_real_distribution<double> doubDist;

};


int main() {
    
    IsingSim sim = IsingSim(100, 0.5, 1309200426);
    auto observs = sim.run(500, 500, 2.0);
    auto mag = std::get<0>(observs);
    auto energy= std::get<1>(observs);

    double avgEnergy = energy.mean();
    std::cout << "Average Energy in One Sweep: " << avgEnergy << std::endl;
    return 0;
};