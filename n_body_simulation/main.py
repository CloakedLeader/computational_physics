import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def csv_to_listofdicts( path: str ):
    df = pd.read_csv(f'{path}', comment='#', header=None, names=['x', 'y', 'vx', 'vy', 'mass'])
    return df.to_dict('records')


def get_time_constraint() -> tuple:
    user_string = input("Enter time step and total time values seperated by a comma, no whitespace: ")
    user_string = user_string.split(',')
    return user_string[0], user_string[1]


class Bodies:

    def __init__( self, data: dict ):
        self.x = data['x']
        self.y = data['y']
        self.vx = data['vx']
        self.vy = data['vy']
        self.mass = data['mass']

    def __repr__( self ):
        return f"Body (x={self.x}, y={self.y}, vx={self.vx}, vy={self.vy}, mass={self.mass})"



def initialise_many_bodies( input: list ) -> None:
    body_dict = []
    for i in input:
        body_dict.append(Bodies(i))
    return body_dict

def create_list_of_dt ( input: tuple ) -> np.ndarray:
    timesteps = np.arange(0, input[1] + input[0], input[0])
    if input[1] % input[0] == 0:
        return timesteps
    else:
        return timesteps[timesteps <= input[1]]
    


list_of_instances = initialise_many_bodies(csv_to_listofdicts('D:\\computational_physics\\n_body_simulation\\bodies.csv'))


