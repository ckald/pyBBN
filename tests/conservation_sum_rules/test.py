import os
import numpy as np

from matplotlib import pyplot as plt
from scipy.integrate import simps


folder = os.path.join(os.path.split(__file__)[0], 'output')

names = [
    ("Electron_neutrino.collision_integrals.txt", "Electron_neutrino.number_and_energy_density.txt", []),
    ("Muon_neutrino.collision_integrals.txt", "Muon_neutrino.number_and_energy_density.txt", []),
    ("Tau_neutrino.collision_integrals.txt", "Tau_neutrino.number_and_energy_density.txt", [])
]

g = 2
scale_factor = (np.loadtxt(folder + "/" + "evolution.txt", unpack=True)[2])[:5000]

for namecol, nameden, integral in names:

    data = np.loadtxt(folder + "/" + namecol, unpack=False)
    momenta = data[0]
    densities = np.loadtxt(folder + "/" + nameden, unpack=True)
    number_den = (densities[0])[:5000]

    # \frac{1}{2\pi^2 a^3}\int{dy y^2 \tilde{I}_{code}} / (na^3)
    for N, integrand in enumerate(data[1:5001]):  # row 1 is momenta
        temp = (g / (2 * np.pi**2) * simps(momenta**2 * integrand, momenta)
                / (number_den[N] * scale_factor[N]**3))
        integral.append(temp)

a, b, integrals_e = names[0]
c, d, integrals_mu = names[1]
e, f, integrals_tau = names[2]


plt.plot(scale_factor[1:5000], (integrals_e), color="red", label="e")
plt.plot(scale_factor[1:5000], (integrals_mu), color="blue", label="mu")
plt.plot(scale_factor[1:5000], (integrals_tau), color="green", label="tau")

plt.xlabel("Scale factor", fontsize=16)
plt.ylabel(r"$\frac{d (n a^3)}{dy (n a^3)} =\ \frac{g}{2\pi^2 (na^3)}\int{dy\ y^2 I_{code}}$", fontsize=16)
plt.title("Elastic processes only")
# plt.semilogy()
# plt.yscale('symlog')
plt.semilogx()
plt.legend(loc="upper right")
plt.show()
input("...")
