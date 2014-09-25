from interaction import Interaction
from particles import Particle, STATISTICS
from evolution import Universe
from common import CONST, UNITS, PARAMS, GRID


PARAMS.T_initial = 3 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-3 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

neutrino_e = Particle(name='Electron neutrino',
                      statistics=STATISTICS.FERMION,
                      dof=4,
                      decoupling_temperature=3 * UNITS.MeV)
Particles.append(neutrino_e)

neutrino_mu = Particle(name='Muon neutrino',
                       statistics=STATISTICS.FERMION,
                       dof=4,
                       decoupling_temperature=3 * UNITS.MeV)
Particles.append(neutrino_mu)


neutrino_e_scattering = Interaction(
    in_particles=[neutrino_e, neutrino_e],
    out_particles=[neutrino_e, neutrino_e],
    decoupling_temperature=0 * UNITS.MeV,
    K1=128 * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3),
    symmetry_factor=0.5
)
Interactions.append(neutrino_e_scattering)

neutrino_mu_scattering = Interaction(
    in_particles=[neutrino_mu, neutrino_mu],
    out_particles=[neutrino_mu, neutrino_mu],
    decoupling_temperature=0 * UNITS.MeV,
    K1=128 * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3),
    symmetry_factor=0.5
)
Interactions.append(neutrino_mu_scattering)

neutrino_e_mu_scattering = Interaction(
    in_particles=[neutrino_e, neutrino_mu],
    out_particles=[neutrino_e, neutrino_mu],
    decoupling_temperature=0 * UNITS.MeV,
    K1=128 * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3),
    symmetry_factor=1.
)
Interactions.append(neutrino_e_mu_scattering)

neutrino_mu_e_scattering = Interaction(
    in_particles=[neutrino_mu, neutrino_e],
    out_particles=[neutrino_e, neutrino_mu],
    decoupling_temperature=0 * UNITS.MeV,
    K1=128 * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3),
    symmetry_factor=1.
)
Interactions.append(neutrino_mu_e_scattering)

import numpy
neutrino_e._distribution += \
    numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-5)**2))(GRID.TEMPLATE)
neutrino_mu._distribution += \
    numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-2.5)**2))(GRID.TEMPLATE)

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[neutrino_e, neutrino_mu])
universe.evolve()

raw_input("...")
