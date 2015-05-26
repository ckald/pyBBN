from common import CONST
from . import non_equilibium_setup, with_setup_args


@with_setup_args(non_equilibium_setup)
def neutrino_scattering_amplitude(params, universe):

    params.update(universe.total_energy_density())

    photon, neutrino_e, neutrino_mu = universe.particles

    assert len(universe.interactions) == 1
    assert len(universe.interactions[0].integrals) == 1
    integral = universe.interactions[0].integrals[0]
    assert len(integral.Ms) == 2

    assert any(M.order == (0, 3, 1, 2) and M.K2 == 0 and M.K1 == 4 * 32 * CONST.G_F**2
               for M in integral.Ms)
    assert any(M.order == (0, 1, 2, 3) and M.K2 == 0 and M.K1 == 2 * 32 * CONST.G_F**2
               for M in integral.Ms)
