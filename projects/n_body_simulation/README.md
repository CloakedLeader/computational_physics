# N-Body Gravitation Simulation

This project simulates a system of particles interacting under Newtonian gravity.
The goal of this project was to explore **numerical stability** and **energy conservation** in integrator methods.

## Methods Used
- Euler Method
- Heun's Method
- Runge-Kutta 4th Order

In the end I decided on using RK4 as it conserves the total energy of the system very well, even though it is slower than the other methods.

Example trajectory of the solar system:

![nbody demo](demo.gif)
