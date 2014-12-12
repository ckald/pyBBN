"""
== C$\nu$B temperature test ==

<img src="plots.png" width=100% />

This test checks that in the universe filled with photons, neutrons, protons, electrons\
and neutrinos (which decouple at temperatures around $\sim 3 MeV$:

  * $a * T$ is not conserved by a factor around `1.401` in presence of additional species

At temperatures $\sim 1 MeV$ electron-positron pairs annihilate, effectively disappearing from the\
plasma. This leads to increasing of the entropy and plasma temperature approximately `1.401` times.

[Log file](log.txt)


"""

from particles import Particle
from evolution import Universe
from common import PARAMS, UNITS
from library import StandardModelParticles as SMP


PARAMS.T_initial = 10 * UNITS.MeV
PARAMS.T_final = 0.05 * UNITS.MeV
PARAMS.dx = 1e-2 * UNITS.MeV
PARAMS.infer()


Particles = []
photon = Particle(**SMP.photon)
neutron = Particle(**SMP.neutron)
proton = Particle(**SMP.proton)
neutrino_e = Particle(**SMP.neutrino_e)
neutrino_mu = Particle(**SMP.neutrino_mu)
neutrino_tau = Particle(**SMP.neutrino_tau)
electron = Particle(**SMP.electron)
muon = Particle(**SMP.muon)
tau = Particle(**SMP.tau)

Particles += [
    photon,
    neutron,
    proton,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
    electron,
    muon,
    tau
]

universe = Universe(Particles, logfile='tests/cosmic_neutrino_temperature/log.txt')
universe.evolve()

print """
    Cosmic photon background temperature is {:.3f} times bigger than cosmic neutrinos temperature.
    Relative error is {:.3f} %
    """.format(PARAMS.aT / UNITS.MeV, (PARAMS.aT / UNITS.MeV - 1.401) / 1.401 * 100)

if universe.graphics:
    universe.graphics.save(__file__)

raw_input("...")
