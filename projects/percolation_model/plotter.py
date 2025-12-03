import pandas
from matplotlib import pyplot as plt

df = pandas.read_csv(filepath_or_buffer=r"G:\computational_physics\projects\percolation_model\results.csv", names=["p", "sp"], comment="#")

ps = df["p"].to_numpy()
sp = df["sp"].to_numpy()

plt.plot(ps, sp, 'rx-')
plt.ylabel("Spanning Probability")
plt.xlabel("Probabilty of Site Inhabitance")
plt.show()
