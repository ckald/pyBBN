import numpy
from nose import with_setup

from common import GRID, PARAMS
from evolution import Universe
from particles import Particle
from library import StandardModelParticles as SMP, StandardModelInteractions as SMI

from . import eps, setup


universe = None


def non_equilibium_setup():
    global universe
    setup()

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.neutrino_e)
    neutrino_mu = Particle(**SMP.neutrino_mu)

    neutrino_self_scattering = SMI.neutrino_self_scattering(neutrino=neutrino_e)

    universe = Universe(particles=[photon, neutrino_e, neutrino_mu],
                        interactions=[neutrino_self_scattering],
                        plotting=False)


@with_setup(non_equilibium_setup)
def free_non_equilibrium_test():

    PARAMS.update(universe.total_energy_density())

    photon, neutrino_e, neutrino_mu = universe.particles

    photon_distribution = photon._distribution
    neutrino_e_distribution = neutrino_e._distribution
    neutrino_mu_distribution = neutrino_mu._distribution

    universe.update_particles()
    universe.init_interactions()
    universe.calculate_collisions()

    assert all(photon.collision_integral == 0), "Equilibrium particle integral is non-zero"
    assert all(neutrino_e.collision_integral * PARAMS.dx < eps), "Integral do not cancel"
    assert all(neutrino_mu.collision_integral == 0), "Free particle integral is non-zero"

    universe.update_distributions()

    print photon.collision_integral, neutrino_e.collision_integral, neutrino_mu.collision_integral
    assert all(photon.collision_integral
               + neutrino_e.collision_integral
               + neutrino_mu.collision_integral == 0), "Collision integrals were not cleared"

    assert all(photon._distribution == photon_distribution),\
        "Equilibrium particle distribution changed"
    assert all(neutrino_e._distribution - neutrino_e_distribution < eps),\
        "Interacting particle distribution changed"
    assert all(neutrino_mu._distribution == neutrino_mu_distribution),\
        "Free particle distribution changed"
