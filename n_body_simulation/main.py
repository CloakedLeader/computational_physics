import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def csv_to_listofdicts( path: str ):
    df = pd.read_csv(f'{path}', comment='#', header=None, names=['x', 'y', 'vx', 'vy', 'mass'])
    return df.to_dict('records')


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

list_of_instances = initialise_many_bodies(csv_to_listofdicts('D:\\computational_physics\\n_body_simulation\\bodies.csv'))


