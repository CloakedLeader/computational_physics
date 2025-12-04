import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/workspaces/computational_physics/projects/forest_fire_model/results.csv", names=["timestep", "empty", "tree", "burning"], comment="#")

t = df["timestep"].to_numpy()
e = df["empty"].to_numpy()
trees = df["tree"].to_numpy()
b = df["tree"].to_numpy()

plt.plot(t, e, "r")
plt.plot(t, trees, "g")
plt.plot(t, b, "b")
plt.grid()
plt.xlabel("# of Timesteps")
plt.title("Simple Plot")
plt.savefig("plot.png")