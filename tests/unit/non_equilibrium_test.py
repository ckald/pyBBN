import numpy
from collections import defaultdict
import environment
import os
from . import non_equilibium_setup, with_setup_args, setup
from common import CONST, UNITS
from evolution import Universe
from particles import Particle
from library.SM import particles as SMP
from library.NuMSM import particles as NuP, interactions as NuI
from interactions.four_particle.cpp.integral import CollisionIntegralKind

@with_setup_args(non_equilibium_setup)
def four_particle_free_non_equilibrium_test(params, universe):
    params.update(universe.total_energy_density(), universe.total_entropy())
    eps = 1e-14
    photon, neutrino_e, neutrino_mu = universe.particles

    photon_distribution = photon._distribution
    neutrino_e_distribution = neutrino_e._distribution
    neutrino_mu_distribution = neutrino_mu._distribution

    universe.update_particles()
    universe.init_interactions()
    universe.calculate_collisions()

    assert all(photon.collision_integral == 0), "Equilibrium particle integral is non-zero"
    assert all(numpy.abs(neutrino_e.collision_integral * params.h) < eps), "Integrals do not cancel"
    assert all(neutrino_mu.collision_integral == 0), "Free particle integral is non-zero"

    universe.update_distributions()

    assert all(photon._distribution == photon_distribution),\
        "Equilibrium particle distribution changed"
    assert all(neutrino_e._distribution - neutrino_e_distribution < eps),\
        "Interacting particle distribution changed"
    assert all(neutrino_mu._distribution == neutrino_mu_distribution),\
        "Free particle distribution changed"


@with_setup_args(non_equilibium_setup)
def unit_non_equilibrium_test(params, universe):
    params.update(universe.total_energy_density(), universe.total_entropy())
    photon, neutrino_e, neutrino_mu = universe.particles

    universe.update_particles()
    universe.init_interactions()

    integral = neutrino_e.collision_integrals[0]

    if environment.get('SPLIT_COLLISION_INTEGRAL'):
        A, B = integral.integrate(neutrino_e.grid.TEMPLATE)
        collision_integral = (A + neutrino_e._distribution * B)
    else:
        collision_integral = integral.integrate(neutrino_e.grid.TEMPLATE)

    universe.calculate_collisions()

    assert numpy.allclose(collision_integral - neutrino_e.collision_integral, 0)


@with_setup_args(setup)
def four_particle_decay_test(params):
    os.environ['SPLIT_COLLISION_INTEGRAL'] = ''

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=200 * UNITS.MeV))

    theta = 1e-4
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[],
        kind = CollisionIntegralKind.F_f_vacuum_decay
    )

    for inter in interaction:
        inter.integrals = [integral for integral in inter.integrals
                            if sum(reactant.side for reactant in integral.reaction) in [2]]

    universe = Universe(params=params)
    universe.add_particles([photon, neutrino_e, sterile])
    universe.interactions += interaction

    params.update(universe.total_energy_density(), universe.total_entropy())

    universe.update_particles()
    universe.init_interactions()

    collision_integral = sterile.collision_integrals[0].integrate(sterile.grid.TEMPLATE)

    theo_value = (CONST.G_F * theta)**2 * sterile.mass**5 / (192 * numpy.pi**3) / UNITS.MeV

    decay_rate = -(collision_integral * sterile.params.H \
                * sterile.conformal_energy(sterile.grid.TEMPLATE) / sterile.conformal_mass \
                / sterile._distribution) / UNITS.MeV

    ratio = decay_rate / theo_value

    assert any(numpy.abs(val) - 1 < 1e-2 for val in ratio), "Four-particle decay test failed"


@with_setup_args(setup)
def three_particle_free_non_equilibrium_test(params):
    eps = 1e-14

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=140 * UNITS.MeV))
    neutral_pion = Particle(**SMP.hadrons.neutral_pion)

    theta = 1e-3
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],
        leptons = [],
        mesons = [neutral_pion]
    )

    universe = Universe(params=params)
    universe.add_particles([photon, neutrino_e, neutrino_mu, sterile, neutral_pion])
    universe.interactions += interaction

    params.update(universe.total_energy_density(), universe.total_entropy())

    universe.update_particles()
    universe.init_interactions()

    photon_distribution = photon._distribution
    neutrino_e_distribution = neutrino_e._distribution
    neutrino_mu_distribution = neutrino_mu._distribution
    sterile_distribution = sterile._distribution
    neutral_pion_distribution = neutral_pion._distribution

    universe.calculate_collisions()

    assert all(photon.collision_integral == 0), "Equilibrium particle integral is non-zero"
    assert all(numpy.abs(neutrino_e.collision_integral * params.h) < eps), "Integrals do not cancel"
    assert all(numpy.abs(sterile.collision_integral * params.h) < eps), "Integrals do not cancel"
    assert all(numpy.abs(neutral_pion.collision_integral * params.h) < eps), "Integrals do not cancel"
    assert all(neutrino_mu.collision_integral == 0), "Free particle integral is non-zero"

    universe.update_distributions()

    assert all(photon._distribution == photon_distribution),\
        "Equilibrium particle distribution changed"
    assert all(neutrino_e._distribution - neutrino_e_distribution < eps),\
        "Interacting particle distribution changed"
    assert all(sterile._distribution - sterile_distribution < eps),\
        "Interacting particle distribution changed"
    assert all(neutral_pion._distribution - neutral_pion_distribution < eps),\
        "Interacting particle distribution changed"
    assert all(neutrino_mu._distribution == neutrino_mu_distribution),\
        "Free particle distribution changed"


@with_setup_args(setup)
def three_particle_decay_test(params):
    os.environ['SPLIT_COLLISION_INTEGRAL'] = ''

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=200 * UNITS.MeV))
    neutral_pion = Particle(**SMP.hadrons.neutral_pion)

    theta = 1e-2
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],
        leptons = [],
        mesons = [neutral_pion],
        kind = CollisionIntegralKind.F_f_vacuum_decay
    )

    universe = Universe(params=params)
    universe.add_particles([photon, neutrino_e, sterile, neutral_pion])
    universe.interactions += interaction

    params.update(universe.total_energy_density(), universe.total_entropy())

    universe.update_particles()
    universe.init_interactions()

    collision_integral = sterile.collision_integrals[0].integrate(sterile.grid.TEMPLATE)

    theo_value = (CONST.G_F * theta * neutral_pion.decay_constant)**2 \
                * sterile.mass**3 * (1 - (neutral_pion.mass / sterile.mass)**2)**2 \
                / (16 * numpy.pi) / UNITS.MeV

    decay_rate = -(collision_integral * sterile.params.H \
                * sterile.conformal_energy(sterile.grid.TEMPLATE) / sterile.conformal_mass \
                / sterile._distribution) / UNITS.MeV

    ratio = decay_rate / theo_value

    assert any(numpy.abs(val) - 1 < 1e-2 for val in ratio), "Three-particle decay test failed"