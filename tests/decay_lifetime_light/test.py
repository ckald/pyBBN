import os
import numpy as np
import sys

#sys.path.insert(0, '/home/nashwan/Master Project/pyBBN')
#from common import UNITS, CONST

from matplotlib import pyplot as plt
from scipy.integrate import simps, cumtrapz

folder = os.path.join(os.path.split(__file__)[0], 'output')
#print("\n\n\n", folder, "\n\n\n")
length = 3500
"""
names = ["Sterile_neutrino_(Dirac).decay_rate.txt", "Sterile_neutrino_(Dirac).decay_rate_nu.txt"]

T = (np.loadtxt(folder + "/" + "evolution.txt", unpack=True)[1])[0:length]
T_nu = (np.loadtxt(folder + "/" + "evolution_nu.txt", unpack=True)[1])[0:length]

decay_rate = np.loadtxt(folder + "/" + names[0], unpack=False)[0:length]
decay_rate_nu = np.loadtxt(folder + "/" + names[1], unpack=False)[0:length]

plt.plot(T , decay_rate, color="red")
plt.plot(T_nu , decay_rate_nu, color="blue")

plt.xlabel("Temperature [MeV]", fontsize=16)
plt.ylabel("Decay rate [MeV]", fontsize=16)
#plt.title("Elastic processes only")
plt.semilogy()
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.semilogx()
plt.gca().invert_xaxis()
#plt.legend(loc="upper right")
plt.show(block=True)

"""

from matplotlib import pyplot as plt
from scipy.integrate import simps, cumtrapz
import matplotlib

folder = os.path.join(os.path.split(__file__)[0], 'output')
#print("\n\n\n", folder, "\n\n\n")
#row = 1900
Dolgov_value = 1.27031e-21 # MeV^-1

names = ["Sterile_neutrino_(Dirac).decay_rate2.txt"]

#T_dec = (np.loadtxt(folder + "/" + "evolution.txt", unpack=True)[1])[0:length]
parameters = np.arange(0.4, 20.4, 0.4)
norm = matplotlib.colors.Normalize(
    vmin=np.min(parameters),
    vmax=np.max(parameters))
c_m = matplotlib.cm.copper
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])





data = np.loadtxt(folder + "/" + names[0], unpack=True)
decay_rates = data[2:,:]
T = data[0]
decay_rate_mean = np.mean(decay_rates, axis=0)

for j in range(len(decay_rates[:,0])):

    decay_rate = decay_rates[j,:]
#    plt.plot(T , decay_rate, color="red")
    plt.plot(T[10:1300], decay_rate[10:1300] / Dolgov_value, color=s_m.to_rgba(parameters[j]))
#    print(decay_rate)
cbar = plt.colorbar(s_m)
cbar.set_label(r"$p_{HNL}$ [MeV]", labelpad=+25, rotation=270)
plt.xlabel("Temperature [MeV]", fontsize=16)
plt.ylabel(r"$\Gamma_{CM, code} / \Gamma_{CM, Dolgov}$", fontsize=16)
#plt.ylabel(r"$\Gamma_{CM}$ [MeV]", fontsize=16)
plt.title("Ratio between decay rate calculated by code and theoretical one")
#plt.semilogy()
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.semilogx()
plt.gca().invert_xaxis()
#plt.legend(loc="upper right")
#plt.show(block=True)
"""


#plt.figure()
names = [folder + '/' + "/Sterile_neutrino_(Dirac).collision_integral.txt"]
data = np.loadtxt(names[0], unpack=True)
T = (data[0])
col = (data[1])
n_c = (data[2])
plt.figure()


plt.subplot(211)
plt.plot(T, -col)
plt.ylabel("col")
plt.semilogx()
plt.semilogy()
plt.gca().invert_xaxis()

plt.subplot(212)
plt.plot(T, n_c)
plt.ylabel("n_c")
plt.semilogx()
plt.semilogy()
plt.xlabel("T")
plt.gca().invert_xaxis()
"""

"""

names = [folder + '/' + "/Sterile_neutrino_(Dirac).collision_integral.txt"]
data = np.loadtxt(names[0], unpack=False)
mom = data[0,2:] 
T = data[1:,0]
time = data[1:,1]
dist = data[1:, 2:]
dfdt = []
for i in np.arange(1,len(time),1):
    dfdt.append((dist[i,:] - dist[i-1,:]) / (-time[i] + time[i-1]))
Gamma = dfdt / dist[:-1,:]
GammaMean = np.mean(Gamma, axis=1)



names2 = [folder + '/' + "/Sterile_neutrino_(Dirac).collision_integral2.txt"]
data2 = np.loadtxt(names[0], unpack=False)
mom2 = data[0,2:] 
T2 = data[1:,0]
time2 = data[1:,1]
dist2 = data[1:, 2:]
dfdt2 = []
for i in np.arange(1,len(time2),1):
    dfdt2.append((dist2[i,:] - dist2[i-1,:]) / (-time2[i] + time2[i-1]))
Gamma2 = dfdt2 / dist2[:-1,:]
GammaMean2 = np.mean(Gamma2, axis=1)




plt.plot(T[:-1], Gamma[:,0], color="red")

plt.plot(T2[:-1], GammaMean2, color="blue")
"""

#plt.xlabel("Temperature [MeV]", fontsize=16)
#plt.ylabel("Decay rate [MeV]", fontsize=16)
#plt.title("Elastic processes only")
#plt.semilogy()
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.semilogx()

#plt.legend(loc="upper right")

plt.show(block=True)
