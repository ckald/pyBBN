from particles import Particle, STATISTICS
from evolution import Universe
from common import UNITS, Params
from library.SM import particles as SMP


params = Params(T_initial=100 * UNITS.MeV,
                T_final=0.1 * UNITS.MeV,
                dx=1e-2 * UNITS.MeV)


photon = Particle(params=params, **SMP.photon)
neutrino = Particle(params=params,
                    name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=2 * UNITS.MeV)

neutrino_noneq = Particle(params=params,
                          name='Neutrino Noneq',
                          statistics=STATISTICS.FERMION,
                          dof=4,
                          decoupling_temperature=75 * UNITS.MeV)

Particles = [photon, neutrino, neutrino_noneq]

universe = Universe(params=params)
universe.particles += Particles
universe.graphics.monitor([neutrino, neutrino_noneq])
universe.evolve()

for particle in Particles:
    print particle

raw_input("...")
