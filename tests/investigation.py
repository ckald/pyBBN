from particles import Particle, STATISTICS
from evolution import Universe
from common import UNITS, PARAMS


PARAMS.T_initial = 100 * UNITS.MeV
PARAMS.T_final = 1e-1 * UNITS.MeV
PARAMS.dx = 1e-2 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=2 * UNITS.MeV)
Particles.append(neutrino)

neutrino_noneq = Particle(name='Neutrino Noneq',
                          statistics=STATISTICS.FERMION,
                          dof=4,
                          decoupling_temperature=75 * UNITS.MeV)
Particles.append(neutrino_noneq)


universe = Universe(Particles)
universe.graphics.monitor([neutrino, neutrino_noneq])
universe.evolve(dx=PARAMS.dx, T_final=PARAMS.T_final)

for particle in Particles:
    print particle

raw_input("...")
