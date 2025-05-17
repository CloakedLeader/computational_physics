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

    G = 6.67430e-11

    _instances = []


    def __init__( self, data: dict, name: str ) -> None:
        self.x = data['x']
        self.y = data['y']
        self.vx = data['vx']
        self.vy = data['vy']
        self.mass = data['mass']
        self.name = name
        Bodies._instances.append(self)

    def __repr__( self ) -> str:
        return f"Body (x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy}, mass={self.mass})"
    
    def all_instances(cls):
        return cls._instances
    
    def dist( self, other: "Bodies" ) -> tuple[ float, float, float ]:
        dx = other.x - self.x
        dy = other.y - self.y
        return math.sqrt( dx ** 2 + dy ** 2), dx, dy
    
    
    def net_grav_force( self, others: List[ "Bodies" ] ) -> tuple[ float, float, float ]:
        fx_total = 0.0
        fy_total = 0.0
        for other in others:
            if other is self:
                continue
            r, dx, dy = self.dist( other )
            if r == 0:
                continue
            force_mag = Bodies.G * self.mass * other.mass / r ** 2
            fx_total += force_mag * dx / r
            fy_total += force_mag * dy / r
            self.force_vec = (fx_total, fy_total)
            self.force_mag = force_mag
            return force_mag, fx_total, fy_total
    
    
    def net_force_vec( self, others: List[ "Bodies" ] ) -> tuple[ float, float ]:
        self.force_vec = self.net_grav_force( others )[1:]
        return self.force_vec
    
    
    def net_force_mag( self, others: List[ "Bodies" ] ) -> float:
        self.force_mag = self.net_grav_force( others )[0]
        return self.force_mag 
    
    def accel( self ) -> tuple[ float, float ]:
        pass



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
    


