import os
import numpy as np
from matplotlib import pyplot as plt

folder = os.path.join(os.path.split(__file__)[0], 'output')

data = np.loadtxt(folder + "/" + "Sterile_neutrino_(Dirac)Densities.txt", unpack=True)
T = data[0]
scalef = data[1]
density_HNL = data[2]
density_muon = data[3]
density_neutrino = data[4]

ll = len(T)
plt.plot(T[0:ll], density_HNL[0:ll], color="blue", linestyle="-", alpha=0.5, label="HNL")
plt.plot(T[0:ll], density_muon[0:ll], color="red", linestyle="-", alpha=0.5, label="Muon")
plt.plot(T[0:ll], 2 * density_neutrino[0:ll], color="green", linestyle="--", alpha=0.5, label="Muon neutrino")

plt.xlabel("Temperature [MeV]", fontsize=16)
plt.ylabel(r"n$_\mathrm{comoving}$", fontsize=16)
plt.legend()
plt.semilogx()
plt.gca().invert_xaxis()

plt.show(block=True)