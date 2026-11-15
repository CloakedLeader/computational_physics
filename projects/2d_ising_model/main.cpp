#include <Eigen/Dense> 
#include <random>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

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
    
    SpinLattice(Eigen::ArrayXXi&& grid) : lattice(std::move(grid)), L(lattice.rows()) {};
    
    const Eigen::ArrayXXi getLattice() const {
        return lattice;
    }

    int size() {
        return L;
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
    IsingSim(int size, unsigned int seed, double upProb = 0.5) : latticeSize(size), rng(seed), dist(upProb), intDist(0, size - 1), doubDist(0, 1.0) {
        Eigen::ArrayXXi grid = Eigen::ArrayXXi::NullaryExpr(
            latticeSize, latticeSize,
            [this]() { return dist(rng) ? 1 : -1; }
        );
        lattice = SpinLattice(std::move(grid));
    };

    IsingSim(SpinLattice grid, unsigned int seed) : lattice(grid), rng(seed) {};

    int size() {
        return lattice.size();
    }
    
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
        double dE = 2 * lattice(pos) * adjacentSum;
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
            double dE = computeEnergyChange(pos);
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

double getExactAvgEnergy(int size, double T) {
    double Z = 0.0;
    double weightedEnergy = 0.0;
    int N = size * size;
    for (std::uint32_t state = 0; state < (1u << N); ++state) {

        Eigen::ArrayXXi lattice(size, size);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {

                int bit = i * size + j;

                lattice(i, j) =
                    (state & (1u << bit)) ? 1 : -1;
            }
        }
        
        // double M = std::abs(lattice.sum());
        double E = SpinLattice(std::move(lattice)).energy();
        double weight = std::exp(-E / T);

        Z += weight;
        weightedEnergy += E * weight;
    }
    return weightedEnergy / Z;

}

std::tuple<double, double> compareExactSmallLattice(int size, double T, int sweeps) {
    IsingSim sim = IsingSim(size, 1309200426, 0.5);
    auto observs = sim.run(sweeps, sweeps, T);
    auto energies = std::get<1>(observs);
    double avgEnergy = energies.mean();

    double exactAverageEnergy = getExactAvgEnergy(size, T);

    std::cout << "Metroplis Average Energy: " << avgEnergy << '\n';
    std::cout << "Exact Average Energy: " << exactAverageEnergy << std::endl;

    return {avgEnergy, exactAverageEnergy};
};

void testSweepConvergence(int size, double T) {
    std::vector<int> sweepCounts = {
        100, 500, 1000, 5000, 10000, 50000, 100000
    };
    
    constexpr int numRuns = 20;
    std::vector<unsigned int> seeds;
    std::mt19937 rng(1309200426);
   
    for (int i = 0; i < numRuns; i++) {
        seeds.push_back(rng());
    }
    double exactE = getExactAvgEnergy(size, T);
    std::cout << std::setprecision(10);
    std::cout << "Sweeps, Mean, StdDev, StdError, AbsError, z \n";

    for (int sweeps : sweepCounts) {
        Eigen::ArrayXd avgEnergies(numRuns);
        
        for (int i = 0; i < numRuns; ++i) {
            IsingSim sim(size, seeds[i], 0.5);
            auto observs = sim.run(10000, sweeps, T);
            auto energy = std::get<1>(observs);
            avgEnergies[i] = energy.mean();
        }
        double avgEnergy = avgEnergies.mean();
        double variance =
            (avgEnergies - avgEnergy).square().sum()
            / (numRuns - 1);
        double stdEnergy = std::sqrt(variance);
        double stdError = stdEnergy / std::sqrt(numRuns);
        double absError = std::abs(avgEnergy - exactE);
        double z = std::abs(avgEnergy - exactE) / stdError;
        
        std::cout << sweeps << ", "
                  << avgEnergy << ", "
                  << stdEnergy << ", "
                  << stdError << ", "
                  << absError << ", "
                  << z << '\n';
    }
}

void runIsingSim(int size, double T, int sweeps) {
    IsingSim sim(size, 1309200426, 0.5);
    auto observs = sim.run(sweeps, sweeps, T);
    auto energy = std::get<1>(observs);

    double avgEnergyPerSpin = energy.mean() / (sim.size() * sim.size());
    std::cout << "Average Energy Per Sweep in One Sweep: " << avgEnergyPerSpin << std::endl;
}

void testOnsager(int size) {
    std::vector<double> temps = {
        1.5, 1.7, 1.9, 2.0, 2.5,
        2.7, 3.0, 3.5, 4.0
    };

    constexpr int numRuns = 20;
    std::vector<unsigned int> seeds;
    std::mt19937 rng(1309200426);
   
    for (int i = 0; i < numRuns; i++) {
        seeds.push_back(rng());
    }

    constexpr int equilSweeps = 10000;
    constexpr int measureSweeps = 10000;

    std::ofstream file("onsager_results.csv");
    if (!file) {
        throw std::runtime_error("Could not open .csv file");
    }

    file << "T, EnergyPerSpin, StdError\n";
    file << std::setprecision(12);
    for (double T : temps) {
        std::cout << "Starting T = " << T << std::endl;
        Eigen::ArrayXd avgEnergies(numRuns);
        for (int i = 0; i < numRuns; i++){
            IsingSim sim(size, seeds[i], 0.5);
            auto observs = sim.run(equilSweeps, measureSweeps, T);
            auto energy = std::get<1>(observs);

            double avgEnergyPerSpin = energy.mean() / (size * size);
            avgEnergies[i] = avgEnergyPerSpin;
        }
            double energyPerSpin = avgEnergies.mean();
            double variance =
                (avgEnergies - energyPerSpin).square().sum()
                / (numRuns - 1);
            double stdEnergy = std::sqrt(variance);
            double stdError = stdEnergy / std::sqrt(numRuns);

            file << T << ", " << energyPerSpin << ", " << stdError << '\n';
            file.flush();
            std::cout << "Finished T = " << T << ", Energy/spin = " << energyPerSpin << std::endl;
    }
        std::cout << "All simulations are complete\n";
}

int main() {
    testOnsager(100);
    return 0;
};