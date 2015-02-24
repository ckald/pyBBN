from evolution import Universe
from particles import Particle
from library.SM import particles as SMP, interactions as SMI

from . import eps, setup, with_setup_args


universe = None


def non_equilibium_setup():
    global universe
    args, kwargs = setup()
    params = args[0]

    photon = Particle(params=params, **SMP.photon)
    neutrino_e = Particle(params=params, **SMP.neutrino_e)
    neutrino_mu = Particle(params=params, **SMP.neutrino_mu)

    neutrino_self_scattering = SMI.neutrino_scattering(neutrino_e, neutrino_e)

    universe = Universe(params=params, plotting=False)
    universe.particles += [photon, neutrino_e, neutrino_mu]
    universe.interactions += [neutrino_self_scattering]

    return args, kwargs


@with_setup_args(non_equilibium_setup)
def free_non_equilibrium_test(params):

    params.update(universe.total_energy_density())

    photon, neutrino_e, neutrino_mu = universe.particles

    photon_distribution = photon._distribution
    neutrino_e_distribution = neutrino_e._distribution
    neutrino_mu_distribution = neutrino_mu._distribution

    universe.update_particles()
    universe.init_interactions()
    universe.calculate_collisions()

    assert all(photon.collision_integral == 0), "Equilibrium particle integral is non-zero"
    assert all(neutrino_e.collision_integral * params.dx < eps), "Integral do not cancel"
    assert all(neutrino_mu.collision_integral == 0), "Free particle integral is non-zero"

    universe.update_distributions()

    assert all(photon._distribution == photon_distribution),\
        "Equilibrium particle distribution changed"
    assert all(neutrino_e._distribution - neutrino_e_distribution < eps),\
        "Interacting particle distribution changed"
    assert all(neutrino_mu._distribution == neutrino_mu_distribution),\
        "Free particle distribution changed"
