from . import non_equilibium_setup, with_setup_args, setup
from common import CONST, UNITS
from evolution import Universe
from particles import Particle
from library.SM import particles as SMP
from library.NuMSM import particles as NuP, interactions as NuI


@with_setup_args(non_equilibium_setup)
def neutrino_scattering_amplitude_test(params, universe):

    params.update(universe.total_energy_density(), universe.total_entropy())

    photon, neutrino_e, neutrino_mu = universe.particles

    assert len(universe.interactions) == 1
    assert len(universe.interactions[0].integrals) == 1
    integral = universe.interactions[0].integrals[0]
    assert len(integral.Ms) == 2

    assert any(M.order == (0, 3, 1, 2) and M.K2 == 0 and M.K1 == 4 * 32 * CONST.G_F**2
               for M in integral.Ms)
    assert any(M.order == (0, 1, 2, 3) and M.K2 == 0 and M.K1 == 2 * 32 * CONST.G_F**2
               for M in integral.Ms)


@with_setup_args(setup)
def three_particle_integral_test(params):

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=200 * UNITS.MeV))
    neutral_pion = Particle(**SMP.hadrons.neutral_pion)

    pion_interactions = NuI.sterile_pion_neutral(theta=1e-3, sterile=sterile,
                                                 active=neutrino_e, pion=neutral_pion)

    universe = Universe(params=params)
    universe.add_particles([photon, neutrino_e, sterile, neutral_pion])
    universe.interactions += pion_interactions

    params.update(universe.total_energy_density())

    print universe.interactions
    assert len(universe.interactions) == 2
    print universe.interactions[0].integrals
    print universe.interactions[1].integrals
    assert len(universe.interactions[0].integrals) == 2
    assert len(universe.interactions[1].integrals) == 3
    integral = universe.interactions[0].integrals[0]
    assert len(integral.Ms) == 1
