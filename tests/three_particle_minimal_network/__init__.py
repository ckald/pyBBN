# -*- coding: utf-8 -*-
"""
## Minimal network with two-body decays

<img src="plots.png" width=100% />
<img src="particles.png" width=100% />

Study a network of Boltzmann equation to get control over instabilities related to the addition of
decays. For simplicity investigate a minimal network which exhibits the and has the following
features:

  * ms > mπ so the decay channel into pions opens
  * neglect back-reaction → standard cosmological evolution
  * minimal set of interactions for equilibration of sterile neutrinos and all other species
  * no mixing

Take into account sterile neutrinos (single generation) νs, pions π0 (may be treated as photons)
and a single generation of active neutrinos να (possibly also treated as photons)

[Log file](log.txt)
[Distribution functions](distributions.txt)


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

params = Params(T_initial=200 * UNITS.MeV,
                T_final=50 * UNITS.MeV,
                dy=0.025)

universe = Universe(params=params, logfile=os.path.join(folder, 'log.txt'))

photon = Particle(params=params, **SMP.photon)

electron = Particle(params=params, **SMP.leptons.electron)

neutrino_e = Particle(params=params, **SMP.leptons.neutrino_e)

neutral_pion = Particle(params=params, **SMP.hadrons.neutral_pion)

sterile = Particle(params=params, **NuP.sterile_neutrino(300 * UNITS.MeV))
sterile.decoupling_temperature = params.T_initial

universe.particles += [
    photon,

    electron,

    neutrino_e,

    neutral_pion,

    sterile,
]

thetas = defaultdict(float, {
    'electron': 1e-3,
})

universe.interactions += (
    SMI.neutrino_interactions(
        leptons=[electron],
        neutrinos=[neutrino_e]
    ) +
    NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[electron]
    )
    + NuI.sterile_hadrons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[electron],
        hadrons=[neutral_pion]
    )
)

universe.graphics.monitor([
    (neutrino_e, RadiationParticleMonitor),
    (sterile, RadiationParticleMonitor),
    (sterile, MassiveParticleMonitor),
])


universe.evolve()

universe.graphics.save(__file__)


""" ## Plots for comparison with articles """

plt.ion()

""" ### JCAP10(2012)014, Figure 9
    <img src="figure_9.png" width=100% /> """

plt.figure(9)
plt.title('Figure 9')
plt.xlabel('MeV/T')
plt.ylabel(u'aT')
plt.xscale('log')
plt.xlim(UNITS.MeV/universe.params.T_initial, UNITS.MeV/universe.params.T_final)
plt.plot(UNITS.MeV / numpy.array(universe.data['T']), numpy.array(universe.data['aT']) / UNITS.MeV)
plt.show()
plt.savefig(os.path.join(folder, 'figure_9.png'))

""" ### JCAP10(2012)014, Figure 10
    <img src="figure_10.png" width=100% />
    <img src="figure_10_full.png" width=100% /> """

plt.figure(10)
plt.title('Figure 10')
plt.xlabel('Conformal momentum y = pa')
plt.ylabel('f/f_eq')
plt.xlim(0, 20)

f_sterile = sterile._distribution
feq_sterile = sterile.equilibrium_distribution()
plt.plot(GRID.TEMPLATE/UNITS.MeV, f_sterile/feq_sterile, label="sterile")

plt.legend()
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10_full.png'))

plt.xlim(0, 10)
# plt.ylim(0.99, 1.06)
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10.png'))

# Distribution functions arrays
distributions_file = open(os.path.join(folder, 'distributions.txt'), "w")
numpy.savetxt(distributions_file, (f_sterile, feq_sterile, f_sterile/feq_sterile),
              header=str(sterile), footer='-'*80, fmt="%1.5e")

distributions_file.close()

# Just to be sure everything is okay
import ipdb
ipdb.set_trace()

raw_input("...")
