from interaction import Interaction
from particles import Particle, STATISTICS
from evolution import Universe
from common import CONST, UNITS, PARAMS, GRID


PARAMS.T_initial = 3 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-1 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=3 * UNITS.MeV
                    )
Particles.append(neutrino)


neutrino_scattering = Interaction(
    in_particles=[neutrino, neutrino],
    out_particles=[neutrino, neutrino],
    decoupling_temperature=0 * UNITS.MeV,
    K1=128 * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3),
    symmetry_factor=0.5
)
Interactions.append(neutrino_scattering)

import numpy
neutrino._distribution += numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-5)**2),
                                          otypes=[numpy.float_])(GRID.TEMPLATE)

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[neutrino])
universe.evolve()

raw_input("...")
