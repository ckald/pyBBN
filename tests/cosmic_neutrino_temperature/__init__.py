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

from particles import Particle, STATISTICS
from evolution import Universe
from common import PARAMS, UNITS


PARAMS.T_initial = 100 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-2 * UNITS.MeV
PARAMS.infer()


Particles = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON,
                  dof=2)
Particles.append(photon)

neutron = Particle(name='Neutron',
                   statistics=STATISTICS.FERMION,
                   mass=0.939 * UNITS.GeV,
                   dof=4)
Particles.append(neutron)

proton = Particle(name='Proton',
                  statistics=STATISTICS.FERMION,
                  mass=0.938 * UNITS.GeV,
                  dof=4)
Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=2,
                    decoupling_temperature=3 * UNITS.MeV)
Particles.append(neutrino)

electron = Particle(name='Electron',
                    mass=0.511 * UNITS.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)

universe = Universe(Particles, logfile='tests/cosmic_neutrino_temperature/log.txt')
universe.evolve(dx=PARAMS.dx, T_final=PARAMS.T_final)

print """
    Cosmic photon background temperature is {:.3f} times bigger than cosmic neutrinos temperature.
    Relative error is {:.3f} %
    """.format(PARAMS.aT / UNITS.MeV, (PARAMS.aT / UNITS.MeV - 1.401) / 1.401 * 100)

universe.graphics.save(__file__)

raw_input("...")
