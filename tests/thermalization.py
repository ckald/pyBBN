from interaction import Interaction
from particles import Particle, STATISTICS
from evolution import Universe
from common import CONST, UNITS, PARAMS, GRID


PARAMS.T_initial = 3 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-6 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

# neutron = Particle(name='Neutron',
#                    statistics=STATISTICS.FERMION,
#                    mass=0.939 * UNITS.GeV,
#                    decoupling_temperature=1.3 * UNITS.MeV)
# Particles.append(neutron)

# proton = Particle(name='Proton',
#                   statistics=STATISTICS.FERMION,
#                   mass=0.938 * UNITS.GeV,
#                   )
# Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=3 * UNITS.MeV
                    )
Particles.append(neutrino)

# electron = Particle(name='Electron',
#                     mass=0.511 * UNITS.MeV,
#                     statistics=STATISTICS.FERMION,
#                     dof=4)
# Particles.append(electron)

neutrino_scattering = Interaction(
    in_particles=[neutrino, neutrino],
    out_particles=[neutrino, neutrino],
    decoupling_temperature=0 * UNITS.MeV,
    K1=lambda (i, j, k, l):
        128 * CONST.G_F**2
        if sorted([i, j]) == [0, 1] and sorted([k, l]) == [2, 3]
        else 0.,
    K2=lambda (i, j, k, l): 0.,
    symmetry_factor=0.5
)
Interactions.append(neutrino_scattering)

# import numpy
# neutrino._distribution += numpy.vectorize(lambda x: 0.3 * numpy.exp(-(x-10)**2),
#                                           otypes=[numpy.float_])(GRID.TEMPLATE)

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[neutrino])
universe.evolve(dx=PARAMS.dx, T_final=PARAMS.T_final)

for particle in Particles:
    print particle

raw_input("...")
