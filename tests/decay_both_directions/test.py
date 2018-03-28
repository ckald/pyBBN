import os
import numpy as np
from matplotlib import pyplot as plt

folder = os.path.join(os.path.split(__file__)[0], 'output')

names = ["Sterile_neutrino_(Dirac)Densities.txt", "Sterile_neutrino_(Dirac)Densities_original.txt"]

# datao = np.loadtxt(folder + "/" + names[1], unpack=True)
# To = datao[0]
# scalefo = datao[1]
# density_HNLo = np.diff(datao[2])
# density_muono = np.diff(datao[3])

# llo = len(To)
# plt.plot(To[1:llo], density_HNLo[0:llo], color="blue", linestyle="--")
# plt.plot(To[1:llo], density_muono[0:llo], color="red", linestyle="--")



data = np.loadtxt(folder + "/" + names[0], unpack=True)
T = data[0]
scalef = data[1]
density_HNL = data[2]
density_muon = data[3]

ll = len(T)
plt.plot(T[0:ll], density_HNL[0:ll], color="blue", label="HNL")
plt.plot(T[0:ll], density_muon[0:ll], color="red", label="Muon")

plt.xlabel("Temperature [MeV]", fontsize=16)
plt.ylabel(r"n$_{i}$ - n$_{i-1}$", fontsize=16)
#plt.semilogy()
#plt.axis(ymin=0, ymax=2, xmin=0.1, xmax=50)
plt.title("Difference in comoving number density between step i and i-1")
plt.legend()
plt.semilogx()
plt.gca().invert_xaxis()

plt.show(block=True)
