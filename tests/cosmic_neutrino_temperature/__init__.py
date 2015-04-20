"""
## C$\nu$B temperature test

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
from common import Params, UNITS
from library.SM import particles as SMP

params = Params(T_initial=10 * UNITS.MeV,
                T_final=0.05 * UNITS.MeV,
                dx=1e-2 * UNITS.MeV)

Particles = []
photon = Particle(**SMP.photon)
neutron = Particle(**SMP.hadrons.neutron)
proton = Particle(**SMP.hadrons.proton)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)
electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
tau = Particle(**SMP.leptons.tau)

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

universe = Universe(params=params, logfile='tests/cosmic_neutrino_temperature/log.txt')
universe.add_particles(Particles)
universe.evolve()

print """
    Cosmic photon background temperature is {:.3f} times bigger than cosmic neutrinos temperature.
    Relative error is {:.3f} %
    """.format(universe.params.aT / UNITS.MeV,
               (universe.params.aT / UNITS.MeV - 1.401) / 1.401 * 100)

if universe.graphics:
    universe.graphics.save(__file__)

raw_input("...")
