from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, PARAMS, GRID


PARAMS.T_initial = 5 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 0.5 * 1e-3 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []

photon = Particle(**SMP.photon)
neutrino_e = Particle(**SMP.neutrino_e)
neutrino_mu = Particle(**SMP.neutrino_mu)

Particles += [photon, neutrino_e, neutrino_mu]

Interactions += [SMI.neutrino_scattering(neutrino_e, neutrino_e)]

import numpy
neutrino_e._distribution += numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-5)**2),
                                            otypes=[numpy.float_])(GRID.TEMPLATE)

universe = Universe(Particles, Interactions)
if universe.graphics:
    universe.graphics.monitor(particles=[neutrino_e, neutrino_mu])

# from plotting import plot_points
# for p0 in GRID.TEMPLATE:
#     Interactions[0].initialize()
#     i = neutrino_e.collision_integrals[0]
#     plot_points(i.bounds(p0), p0)

universe.evolve()

raw_input("...")
