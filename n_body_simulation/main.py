import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import math
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


def get_time_constraint() -> tuple:
    user_string = input("Enter time step and total time values seperated by a comma, no whitespace: ")
    user_string = user_string.split(',')
    return user_string[0], user_string[1]


class Bodies:

    G = sp.constants.gravitational_constant


    _instances: List["Bodies"] = []


    def __init__( self, data: dict, name: str ) -> None:
        self.pos = np.array( [ data['x'], data['y'] ], dtype=float )
        self.vel = np.array( [ data['vx'], data['vy'] ], dtype=float )
        self.force = np.zeros(2)
        self.mass = data['mass']
        self.name = name
        Bodies._instances.append(self)

    def __repr__( self ) -> str:
        return f"Body (x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy}, mass={self.mass})"
    
    def all_instances(cls) -> List["Bodies"]:
        return cls._instances
    
    def dist( self, other: "Bodies" ) -> tuple[ float, np.ndarray ]:
        dist = other.pos - self.pos
        r = np.linalg.norn( dist )
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
    

    
    #def net_force_vec( self, others: List[ "Bodies" ] ) -> tuple[ float, float ]:
        #self.force_vec = self.net_grav_force( others )[1:]
        #return self.force_vec
    
    
    #def net_force_mag( self, others: List[ "Bodies" ] ) -> float:
        #self.force_mag = self.net_grav_force( others )[0]
        #return self.force_mag 
    
    def update( self, dt: float ):
        accel = self.force / self.mass
        self.vel += accel * dt
        self.pos += self.vel * dt


class Simulation:

    def __init__( self, bodies: List["Bodies"] ):
        self.bodies = bodies

    def step( self, dt ):
        for body in self.bodies:
            body.net_grav_force(self.bodies)

            for body in self.bodies:
                body.update(dt)

    def run( self, dt: float, steps):
        for  _ in range(steps):
            self.steps(dt)

            
    
            










def initialise_many_bodies( input: list ) -> List:
    body_dict = []
    for i in input:
        input()
        body_dict.append(Bodies(i))
    return body_dict

def create_list_of_dt ( input: tuple ) -> np.ndarray:
    timesteps = np.arange(0, input[1] + input[0], input[0])
    if input[1] % input[0] == 0:
        return timesteps
    else:
        return timesteps[timesteps <= input[1]]
    
def update_bodies( bodies: List["Bodies"], dt: float ) -> None:
    for body in bodies:
        accel = body.force / body.mass
        body.vel += accel * dt
        body.pos += body.vel * dt

def run_simulation( bodies: List["Bodies"], dt: float, steps: float ):

    for step in range(steps):
        pass
        



