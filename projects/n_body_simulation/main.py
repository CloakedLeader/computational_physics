import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import math
from typing import List
import tkinter as tk
from tkinter import messagebox
import random

"""
Structure of dictionary to add bodies:
'x': ...
'y': ...
'vx': ...
'vy': ...
'mass': ...
"""

pi = math.pi


def csv_to_listofdicts(path: str):
    df = pd.read_csv(f'{path}', dtype={'name': str}, comment='#',
                     header=None, names=['x', 'y', 'vx', 'vy', 'mass', 'name'])
    return df.to_dict('records')


# --- Body Class ---
class Bodies:

    G = 4 * pi**2

    body_counter = 1

    _instances: List["Bodies"] = []

    def __init__(self, data: dict) -> None:
        self.pos = np.array([data['x'], data['y']], dtype=float)
        self.vel = np.array([data['vx'], data['vy']], dtype=float)
        self.name = str(data['name'])
        self.force = np.zeros(2)
        self.mass = data['mass']
        self.identifier = Bodies.body_counter
        Bodies.body_counter += 1
        self.history = []
        Bodies._instances.append(self)

    def __repr__(self) -> str:
        return f"Body (x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy}, mass={self.mass})"

    def all_instances(cls) -> List["Bodies"]:
        return cls._instances

    def dist(self, other: "Bodies") -> tuple[float, np.ndarray]:
        dist = other.pos - self.pos
        r = np.linalg.norm(dist)
        return r, dist

    def reset_force(self):
        self.force = np.zeros(2)

    def net_grav_force(self, others: List["Bodies"]) -> None:
        self.reset_force()
        for other in others:
            if other is self:
                continue
            r_mag, r_vec = self.dist(other)
            if r_mag == 0:
                continue
            r_hat = r_vec / r_mag
            force_mag = Bodies.G * self.mass * other.mass / r_mag ** 2
            self.force += force_mag * r_hat

    @property
    def accel(self) -> float:
        return self.force / self.mass

    def verlet_pos_update(self, dt: float):
        self.pos += self.vel * dt + 0.5 * self.accel * dt**2

    def update(self, others: List["Bodies"], dt: float):
        old_accel = self.accel
        self.verlet_pos_update(dt)
        self.net_grav_force(others)
        new_accel = self.accel
        self.vel += 0.5 * (old_accel + new_accel) * dt

    def kin_energy(self):
        vel_mag = np.linalg.norm(self.vel)
        out = float(0.5 * self.mass * vel_mag**2)
        return out

    def pot_energy(self, others: List["Bodies"]):
        total = 0
        for other in others:
            if other is self:
                continue
            r_mag = self.dist(other)[0]
            if r_mag == 0:
                continue
            total += - Bodies.G * self.mass * other.mass / r_mag
            return total


# --- Simulation Class ---
class Simulation:

    def __init__(self, bodies: List["Bodies"]):
        self.bodies = bodies
        self.kin = 0
        self.pot = 0
        self.total = 0

    def list_of_names(self):
        dummy_list = []
        for body in self.bodies:
            dummy_list.append(body.name)
        return dummy_list

    def total_kin_energy(self):
        total = 0
        for body in self.bodies:
            total += body.kin_energy()
        return total

    def total_pot_energy(self):
        total = 0
        for body in self.bodies:
            total += body.pot_energy(self.bodies)
        return total

    def tot_energy(self):
        self.total = self.kin + self.pot

    def step(self, dt):
        for body in self.bodies:
            body.net_grav_force(self.bodies) 

        for body in self.bodies:
            body.update(self.bodies, dt)

        self.kin = self.total_kin_energy()
        self.pot = self.total_pot_energy()
        self.total = self.tot_energy()

    def run(self, dt: float, steps: int):
        self.kin_energy_hist = []
        self.pot_energy_hist = []
        self.total_energy_hist = []
        self.time_history = []
        self.kin = self.total_kin_energy()
        self.pot = self.total_pot_energy()
        self.total = self.tot_energy()
        self.kin_energy_hist.append(self.kin)
        self.pot_energy_hist.append(self.pot)
        self.total_energy_hist.append(self.kin + self.pot)
        self.time_history.append(0)

        for step in range(1, steps + 1):
            for body in self.bodies:
                body.net_grav_force(self.bodies)
            for body in self.bodies:
                body.update(self.bodies, dt)
                body.history.append(body.pos.copy())

            self.kin = self.total_kin_energy()
            self.pot = self.total_pot_energy()
            self.total = abs(self.kin) + abs(self.pot)
            self.kin_energy_hist.append(self.kin)
            self.pot_energy_hist.append(self.pot)
            self.total_energy_hist.append(self.total)
            self.time_history.append(step * dt)


def initialise_many_bodies(input: list) -> List:
    body_list = []
    for i in input:
        body_list.append(Bodies(i))
    return body_list


class NBodyGUI:

    def __init__(self, root):
        self.root = root
        self.root.title("N-Body Simulation Setup")

        tk.Label(root, text="Timestep (dt): ").grid(row=0, column=0)
        self.dt_entry = tk.Entry(root)
        self.dt_entry.grid(row=0, column=1)

        tk.Label(root, text="Number of steps: ").grid(row=1, column=0)
        self.steps_entry = tk.Entry(root)
        self.steps_entry.grid(row=1, column=1)

        tk.Label(root, text="Initial Configuration: ").grid(row=3, column=0)
        self.config_var = tk.StringVar(root)
        self.config_var.set("Simple Solar System")
        options = ["Random", "Simple Solar System", "Three Body Problem", "Binary System"]
        self.dropdown = tk.OptionMenu(root, self.config_var, *options)
        self.dropdown.grid(row=3, column=1)

        run_button = tk.Button(root, text="Run Simulation", command=self.run_sim)
        run_button.grid(row=4, column=0, columnspan=2, pady=10)

    config_files = {
        "Simple Solar System" : r"D:\computational_physics\projects\n_body_simulation\data_files\solar_system.csv",
        "Three Body Problem": r"D:\computational_physics\projects\n_body_simulation\data_files\three_body.csv",
        "Binary System": r"D:\computational_physics\projects\n_body_simulation\data_files\binary_system.csv"
    }

    def run_sim(self):
        try:
            dt = float(self.dt_entry.get())
            steps = int(self.steps_entry.get())
            config = self.config_var.get()

            if config == "Random":
                filename = random.choice(list(NBodyGUI.config_files.values()))

            filename = NBodyGUI.config_files.get(config)
            if not filename:
                raise ValueError("Invalid configuration selected.")
        except ValueError:
            messagebox.showerror("Invalid input", "Please enter valid numbers.")
            return

        dummy = csv_to_listofdicts(filename)
        list_of_bodies = initialise_many_bodies(dummy)
        sim = Simulation(list_of_bodies)
        sim.run(dt, steps)

        bodies_history = []
        for body in list_of_bodies:
            bodies_history.append(body.history)

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)

        scatters = [ax.plot([], [], 'o')[0] for _ in list_of_bodies]
        labels = [ax.text(0, 0, name, fontsize=8, ha='left', va='bottom') for
                  name in sim.list_of_names()]
        lines = [ax.plot([], [], lw=1)[0] for _ in list_of_bodies]

        def init():
            for scatter, line in zip(scatters, lines):
                scatter.set_data([], [])
                line.set_data([], [])
            for label in labels:
                label.set_position((0, 0))
            return scatters + lines + labels

        def update(frame):
            print(f"Frame {frame}:")
            for i, body in enumerate(list_of_bodies):
                history = np.array(body.history)
                lines[i].set_data(history[:frame+1, 0], history[:frame+1, 1])
                x, y = history[frame]
                print(f"Body {i} position: ({x}, {y})")
                scatters[i].set_data([x], [y])
                labels[i].set_position((x+0.05, y+0.05))
            return scatters + lines + labels

        ani = FuncAnimation(
            fig, update, frames=len(list_of_bodies[0].history),
            init_func=init, blit=False, interval=20
        )

        plt.figure(figsize=(10, 6))
        plt.plot(sim.time_history, sim.kin_energy_hist, label="Kinetic Energy", color="blue")
        plt.plot(sim.time_history, sim.pot_energy_hist, label="Potential Energy", color="green")
        plt.plot(sim.time_history, sim.total_energy_hist, label="Total Energy", color="red", linestyle= "--")
        plt.ylabel = "Energies"
        plt.xlabel = "Time"
        plt.title = "Energy vs Time"
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        plt.show()
        print(sim.list_of_names())


if __name__ == "__main__":
    root = tk.Tk()
    app = NBodyGUI(root)
    root.mainloop()

# filename = r"D:\computational_physics\n_body_simulation\solar_system.csv"
# dummy = csv_to_listofdicts(filename)
# print( dummy[1] )
# list_of_bodies = initialise_many_bodies(dummy)
# sim = Simulation(list_of_bodies)
