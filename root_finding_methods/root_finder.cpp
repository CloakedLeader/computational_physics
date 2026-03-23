#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <functional>


class RootFinder {
private:
    std::function<double(double)> func;
    double start;
    double end;
    double epsilon;

    double epsilon = 0.0001; 

public:
    RootFinder(std::function<double(double)> function, double start_, double end_)
        : func(function), start(start_), end(end_)

    double bisectorRootFinder(int iterations) {
        int iter = 0;

        while (iter < iterations) {
            double mid = (start + end) / 2;
            double fc = func(mid);

            if (fc == 0 or (abs(fc) < epsilon)) {
                std::cout << "It took " << iter << " iterations to get to root within " << epsilon;
                return mid; 
            }
            if ((fc * func(start)) < 0 ) {
                end = mid;
            }
            else {
                start = mid;
            }
            iter += 1
        };

        std::cout << "Could not find a root in the specified number of attempts: " << 

    }
      
}