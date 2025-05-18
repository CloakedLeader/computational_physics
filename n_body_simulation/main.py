import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import scipy as sp
from typing import List

"""
Structure of dictionary to add bodies:
'x': ...
'y': ...
'vx': ...
'vy': ...
'mass': ...
"""

def csv_to_listofdicts( path: str ):
    df = pd.read_csv(f'{path}', comment='#', header=None, names=['x', 'y', 'vx', 'vy', 'mass'])
    return df.to_dict('records')

# --- Body Class ---
class Bodies:

    G = sp.constants.gravitational_constant

    body_counter = 1

    _instances: List["Bodies"] = []


    def __init__( self, data: dict) -> None:
        self.pos = np.array( [ data['x'], data['y'] ], dtype=float )
        self.vel = np.array( [ data['vx'], data['vy'] ], dtype=float )
        self.force = np.zeros(2)
        self.mass = data['mass']
        self.identifier = Bodies.body_counter
        Bodies.body_counter += 1
        self.history = []
        Bodies._instances.append(self)

    def __repr__( self ) -> str:
        return f"Body (x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy}, mass={self.mass})"
    
    def all_instances(cls) -> List["Bodies"]:
        return cls._instances
    
    def dist( self, other: "Bodies" ) -> tuple[ float, np.ndarray ]:
        dist = other.pos - self.pos
        r = np.linalg.norm( dist )
        return r, dist
    
    def reset_force( self ):
        self.force = np.zeros(2)
    
    def net_grav_force( self, others: List[ "Bodies" ] ) -> None:
        self.reset_force()
        for other in others:
            if other is self:
                continue
            r_mag, r_vec = self.dist( other )
            if r_mag == 0:
                continue
            r_hat = r_vec / r_mag
            force_mag = Bodies.G * self.mass * other.mass / r_mag ** 2
            self.force += force_mag * r_hat

    
    def update( self, dt: float ):
        accel = self.force / self.mass
        self.vel += accel * dt
        self.pos += self.vel * dt
        self.history.append(self.pos.copy())

# --- Simulation Class ---
class Simulation:

    def __init__( self, bodies: List["Bodies"] ):
        self.bodies = bodies

    def step( self, dt ):
        for body in self.bodies:
            body.net_grav_force(self.bodies)

            for body in self.bodies:
                body.update(dt)

    def run( self, dt: float, steps: int ):
        for _ in range(steps):
            for body in self.bodies:
                body.net_grav_force(self.bodies)
            for body in self.bodies:
                body.update(dt)

            
def initialise_many_bodies( input: list ) -> List:
    body_list = []
    for i in input:
        body_list.append(Bodies(i))
    return body_list



dummy = csv_to_listofdicts(r"D:\computational_physics\n_body_simulation\bodies.csv")
list_of_bodies = initialise_many_bodies(dummy)
sim = Simulation(list_of_bodies)

dt = 10000
steps = 3154
sim.run(dt, steps)

bodies_history = []
for body in list_of_bodies:
    bodies_history.append(body.history)


fig, ax = plt.subplots()
ax.set_aspect( 'equal' )
ax.set_xlim( -5e12, 5e12 )
ax.set_ylim( -5e12, 5e12 )

scatters = [ax.plot([], [], 'o')[0] for _ in list_of_bodies]

def init():
    for scatter in scatters:
        scatter.set_data([], [])
    return scatters

def update(frame):
    print(f"Frame {frame}:")
    for i, body in enumerate(list_of_bodies):
        x, y = body.history[frame]
        print(f"Body {i} position: ({x}, {y})")
        scatters[i].set_data([x], [y])
    return scatters

ani = FuncAnimation(
    fig, update, frames=len(bodies_history[0]),
    init_func=init, blit=False, interval=20
)

plt.show()

print(f"Steps: {steps}")
for i, body in enumerate(list_of_bodies):
    print(f"Body {i} has {len(body.history)} positions.")

