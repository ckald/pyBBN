"""
## CnuB temperature test

<img src="plots.svg" width=100% />

This test checks that in the universe filled with photons, electrons and neutrinos \
(which decouple at temperatures around $\sim 3 MeV$:

  * $a * T$ is not conserved by a factor around `1.401` in presence of additional species

At temperatures $\sim 1 MeV$ electron-positron pairs annihilate, effectively disappearing from the\
plasma. This leads to increasing of the entropy and plasma temperature approximately `1.401` times.

[Log file](log.txt)


"""

import os
from particles import Particle
from evolution import Universe
from common import Params, UNITS
from library.SM import particles as SMP


folder = os.path.join(os.path.split(__file__)[0], 'output')

params = Params(T=10 * UNITS.MeV,
                dy=0.003125)
T_final = 0.008 * UNITS.MeV


Particles = []
photon = Particle(**SMP.photon)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)
electron = Particle(**SMP.leptons.electron)

Particles += [
    photon,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
    electron,
]

universe = Universe(params=params, folder=folder)
universe.add_particles(Particles)
universe.init_kawano(electron=electron, neutrino=neutrino_e)

universe.evolve(T_final)


print ("""
    Cosmic photon background temperature is {:e} times bigger than cosmic neutrinos temperature.
    Relative error is {:e} %
    """.format(universe.params.aT / UNITS.MeV,
               (universe.params.aT / UNITS.MeV - 1.40102) / 1.40102 * 100))
