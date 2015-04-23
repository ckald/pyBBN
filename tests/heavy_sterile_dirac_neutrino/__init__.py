# -*- coding: utf-8 -*-
"""
## Heavy sterile dirac neutrino

$$ M = 33.9 MeV $$

$$ \theta_\tau \approx 7.6 10^{-4} \sim \tau_N \approx 0.3 sec $$

http://arxiv.org/pdf/hep-ph/0002223v2.pdf

<img src="plots.png" width=100% />
<img src="particles.png" width=100% />

"""

import os
import numpy
import matplotlib
from collections import defaultdict

from plotting import plt, RadiationParticleMonitor, MassiveParticleMonitor
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, GRID


folder = os.path.split(__file__)[0]

params = Params(T_initial=100. * UNITS.MeV,
                T_final=1. * UNITS.MeV,
                dy=0.05)

universe = Universe(params=params, logfile=os.path.join(folder, 'log.txt'))

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)
sterile = Particle(**NuP.sterile_neutrino(33.9 * UNITS.MeV))

sterile.decoupling_temperature = params.T_initial

universe.add_particles([
    photon,
    electron,
    muon,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
    sterile,
])

thetas = defaultdict(float, {
    'tau': 7.6 * 1e-4,
})

universe.interactions += (
    SMI.neutrino_interactions(
        leptons=[electron],
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau]
    ) + NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron, muon]
    )
)

universe.graphics.monitor([
    (neutrino_e, RadiationParticleMonitor),
    (neutrino_mu, RadiationParticleMonitor),
    (neutrino_tau, RadiationParticleMonitor),
    (sterile, MassiveParticleMonitor)
])

universe.evolve()

universe.graphics.save(__file__)


""" ## Plots for comparison with articles"""

plt.ion()

"""
### JCAP10(2012)014, Figure 9
<img src="figure_9.png" width=100% />
"""

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
plt.savefig(os.path.join(folder, 'figure_9.png'))

"""
### JCAP10(2012)014, Figure 10
<img src="figure_10.png" width=100% />
<img src="figure_10_full.png" width=100% />
"""

plt.figure(10)
plt.title('Figure 10')
plt.xlabel('Conformal momentum y = pa')
plt.ylabel('f/f_eq')
plt.xlim(0, 20)

# Distribution functions arrays
distributions_file = open(os.path.join(folder, 'distributions.txt'), "w")

for neutrino in [neutrino_e, neutrino_mu, neutrino_tau, sterile]:
    f = neutrino_e._distribution
    feq = neutrino_e.equilibrium_distribution()
    plt.plot(GRID.TEMPLATE/UNITS.MeV, f/feq, label=str(neutrino))

    numpy.savetxt(distributions_file, (f, feq, f/feq), header=str(neutrino),
                  footer='-'*80, fmt="%1.5e")

plt.legend()
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10_full.png'))

plt.xlim(0, 10)
plt.ylim(0.99, 1.06)
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10.png'))

distributions_file.close()

# Just to be sure everything is okay
import ipdb
ipdb.set_trace()

raw_input("...")
