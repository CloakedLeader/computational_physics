import numpy as np
from matplotlib import pyplot as plt
from scipy.special import ellipk as elp
import pandas as pd

def onsager_energy(T, J = 1.0):
    K = J / T
    k = 2 * np.sinh(2*K) / np.cosh(2*K)**2

    m = k**2

    return (
        -J / np.tanh(2*K)
        * (
            1
            + (2 / np.pi)
            * (2 * np.tanh(2*K)**2 - 1)
            * elp(m)
        )
    )


def plot_onsager_errors():
    data = pd.read_csv("/home/adamp/Documents/coding/computational_physics/projects/2d_ising_model/onsager_results_w_errors.csv",
                       skipinitialspace=True)

    T = np.linspace(1.2, 4.0, 500)

    plt.plot(T, onsager_energy(T), label="Onsager exact solution")

    plt.errorbar(
        data["T"],
        data["EnergyPerSpin"],
        yerr=data["StdError"],
        fmt="o",
        markersize=3,
        capsize=4,
        capthick=1,
        elinewidth=1,
        label="Metropolis (20 independent runs)",
        zorder=3
    )

    plt.axvline(
        2 / np.log(1 + np.sqrt(2)),
        linestyle="--",
        color="gray",
        label=r"$T_c$"
    )

    plt.xlabel(r"Temperature $T$")
    plt.ylabel(r"Average energy per spin $\langle E\rangle_N$")
    plt.title("2D Ising Model: Metropolis vs Onsager")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

plot_onsager_errors()

