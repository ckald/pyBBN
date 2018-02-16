import os
import numpy as np
import sys

from matplotlib import pyplot as plt
from scipy.integrate import simps, cumtrapz

from matplotlib import pyplot as plt
from scipy.integrate import simps, cumtrapz
import matplotlib

folder = os.path.join(os.path.split(__file__)[0], 'output')
Theo_value = 3.92472e-21 # MeV

names = ["Sterile_neutrino_(Dirac).decay_rate.txt"]


#parameters = np.arange(0.4, 20.4, 0.4)
#norm = matplotlib.colors.Normalize(
#    vmin=np.min(parameters),
#    vmax=np.max(parameters))
#c_m = matplotlib.cm.copper
#s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
#s_m.set_array([])

data = np.loadtxt(folder + "/" + names[0], unpack=True)
decay_rates = data[3:,:]
T = data[0]
scalef = data[1]

for j in range(len(decay_rates[:,0])):

    decay_rate = decay_rates[45,:]
    plt.plot(T[0:2000], decay_rate[0:2000] / Theo_value, color="red")#s_m.to_rgba(parameters[j]))

#cbar = plt.colorbar(s_m)
#cbar.set_label(r"$p_{HNL}$ [MeV]", labelpad=+25, rotation=270)
plt.xlabel("Temperature [MeV]", fontsize=16)
plt.ylabel(r"$\Gamma_{CM, code} / \Gamma_{CM, Theory}$", fontsize=16)
plt.title("Ratio between decay rate calculated by code and theoretical one")
#plt.semilogy()
#plt.axis(ymin=1.1e-21, ymax=1.5e-21)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.semilogx()
plt.gca().invert_xaxis()
#plt.legend(loc="upper right")

plt.show(block=True)
