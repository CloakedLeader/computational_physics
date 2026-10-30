# N-Body Gravitation Simulation

This project simulates a system of particles interacting under Newtonian gravity.
The goal of this project was to explore **numerical stability** and **energy conservation** in integrator methods.

## Methods Used
- Euler Method
- Verlet Integrator
- Heun's Method

I originally used RK4 since this I assumed that the highest order method will be the most accurate but it led to the system not conserving energy and it slowly creeping up over time. So I did some research and came across sympletic integrators which are made to conserve the volume in phase space and hence the energy of the system. I then implemented a velocity-Verlet method which led to much more physical results.

## Results

As a measure of numerical error I use the maximum relative energy drift.

For a circular orbit the error for different step sizes (0.004, 0.002, 0.001 and 0.0005) all are under $10^{-6}.$



Example trajectory of the solar system:

![nbody demo](demo.gif)
