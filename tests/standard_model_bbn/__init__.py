"""
## Standard Model BBN test

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

This test checks that in the universe filled with photons, electrons and neutrinos:

  * $a * T$ is not conserved by a factor around `1.401` and precise details of this process
  * neutrino non-equilibrium corrections reproduce the results of the Dolgov-Hansen-Semikoz papers

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os
import numpy
import matplotlib

from plotting import plt, RadiationParticleMonitor
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params, GRID


folder = os.path.split(__file__)[0]

params = Params(T_initial=10. * UNITS.MeV,
                T_final=0.0008 * UNITS.MeV)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

neutron = Particle(**SMP.hadrons.neutron)
proton = Particle(**SMP.hadrons.proton)

universe.add_particles([
    photon,
    electron,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    neutron,
    proton
])

universe.interactions += (
    # [SMI.baryons_interaction(neutron=neutron, proton=proton,
    #                          electron=electron, neutrino=neutrino_e)]
    # +
    SMI.neutrino_interactions(leptons=[electron],
                              neutrinos=[neutrino_e, neutrino_mu, neutrino_tau])
)

universe.init_kawano(electron=electron, neutrino=neutrino_e)

universe.graphics.monitor([
    (neutrino_e, RadiationParticleMonitor),
    (neutrino_mu, RadiationParticleMonitor),
    (neutrino_tau, RadiationParticleMonitor)
])


universe.evolve()

""" ## Plots for comparison with articles """

plt.ion()

"""
### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% /> """

plt.figure(9)
plt.title('Figure 9')
plt.xlabel('MeV/T')
plt.ylabel(u'aT')
plt.xscale('log')
plt.xlim(0.5, UNITS.MeV/universe.params.T_final)
plt.xticks([1, 2, 3, 5, 10, 20])
plt.axes().get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.plot(UNITS.MeV / numpy.array(universe.data['T']), numpy.array(universe.data['aT']) / UNITS.MeV)
plt.show()
plt.savefig(os.path.join(folder, 'figure_9.svg'))

"""
### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% /> """

plt.figure(10)
plt.title('Figure 10')
plt.xlabel('Conformal momentum y = pa')
plt.ylabel('f/f_eq')
plt.xlim(0, 20)

# Distribution functions arrays
distributions_file = open(os.path.join(folder, 'distributions.txt'), "w")

for neutrino in (neutrino_e, neutrino_mu, neutrino_tau):
    f = neutrino._distribution
    feq = neutrino.equilibrium_distribution()
    plt.plot(GRID.TEMPLATE/UNITS.MeV, f/feq, label=neutrino.flavour)

    numpy.savetxt(distributions_file, (f, feq, f/feq), header=str(neutrino),
                  footer='-'*80, fmt="%1.5e")

plt.legend()
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10_full.svg'))

plt.xlim(0, 10)
plt.ylim(0.99, 1.06)
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10.svg'))


distributions_file.close()

# Just to be sure everything is okay
# import ipdb
# ipdb.set_trace()

# raw_input("...")
