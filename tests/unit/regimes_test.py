from nose import with_setup

from common import PARAMS, UNITS
from particles import Particle, REGIMES
from library import StandardModelParticles as SMP

from . import eps, setup


@with_setup(setup)
def params_inferrence_test():
    assert PARAMS.T_initial == 5. * UNITS.MeV
    assert PARAMS.T == PARAMS.T_initial
    assert PARAMS.aT == PARAMS.m, "Initial `aT` must be equal to 1"
    assert PARAMS.x - PARAMS.m * PARAMS.aT / PARAMS.T < eps


# @with_setup(setup)
# def params_copy_test():
#     other_params = PARAMS.copy()
#     assert id(other_params) != id(PARAMS), "PARAMS were not copied"
#     other_params.dx *= 1000
#     assert other_params.dx != PARAMS.dx


@with_setup(setup)
def radiation_regime_test():

    photon = Particle(**SMP.photon)
    assert photon.in_equilibrium, "Photon must always stay in equilibrium"
    assert not photon.decoupling_temperature, "Photon can't decouple"
    assert photon.regime == REGIMES.RADIATION, "Photon is a relativistic particle"
    assert photon.mass == 0, "Photon is massless"
    assert photon.eta == -1, "Photon is a boson, it's eta must be equal to -1"
    assert photon.numerator() == 0, "Photon does not contribute to the numerator"
    assert photon.denominator() != 0, "Photon does contribute to the denominator"


@with_setup(setup)
def intermediate_regime_test():

    electron = Particle(**SMP.electron)
    assert electron.in_equilibrium, "Electron must always stay in equilibrium"
    assert not electron.decoupling_temperature, "Electron can't decouple"
    assert electron.regime == REGIMES.INTERMEDIATE, \
        "Electron is not strictly relativistic at {} MeV".format(PARAMS.T/UNITS.MeV)
    assert electron.mass != 0, "Electron is massive"
    assert electron.eta == 1, "Electron is a fermion, it's eta must be equal to 1"
    assert electron.numerator() != 0, "Massive particles contribute to the numerator"
    assert electron.denominator() != 0, "Massive particles contribute to the denominator"


@with_setup(setup)
def dust_regime_test():

    proton = Particle(**SMP.proton)
    assert proton.in_equilibrium, "Proton must always stay in equilibrium"
    assert not proton.decoupling_temperature, "Proton can't decouple"
    assert proton.regime == REGIMES.DUST, \
        "Proton non-relativistic at {} MeV".format(PARAMS.T/UNITS.MeV)
    assert proton.mass != 0, "Proton is massive"
    assert proton.eta == 1, "Proton is a fermion, it's eta must be equal to 1"
    assert proton.numerator() != 0, "Massive particles contribute to the numerator"
    assert proton.denominator() != 0, "Massive particles contribute to the denominator"


@with_setup(setup)
def mass_regime_switching_test():

    electron = Particle(**SMP.electron)
    assert electron.regime == REGIMES.INTERMEDIATE

    electron.mass = 0
    assert electron.regime == REGIMES.RADIATION

    electron.mass = 1e5 * UNITS.MeV
    assert electron.regime == REGIMES.DUST


@with_setup(setup)
def temperature_regime_switching_test():

    electron = Particle(**SMP.electron)
    assert electron.regime == REGIMES.INTERMEDIATE

    PARAMS.T = 100 * UNITS.MeV
    electron.update()
    assert electron.regime == REGIMES.RADIATION

    PARAMS.T = 1 * UNITS.keV
    electron.update()
    assert electron.regime == REGIMES.DUST


@with_setup(setup)
def statistics_consistency_test():

    photon = Particle(**SMP.photon)
    neutrino = Particle(**SMP.neutrino_e)

    assert 7./8. - neutrino.energy_density() / photon.energy_density() < eps


@with_setup(setup)
def decoupling_test():

    PARAMS.T_initial *= 2
    PARAMS.infer()

    neutrino = Particle(**SMP.neutrino_e)

    assert neutrino.in_equilibrium
    eq_distribution = neutrino._distribution

    PARAMS.T_initial /= 2
    PARAMS.infer()

    assert neutrino.in_equilibrium, "Neutrino should not depend on global temperature change"
    neutrino.update()

    assert not neutrino.in_equilibrium, "Neutrino should have decoupled"
    noneq_distribution = neutrino._distribution

    assert all(eq_distribution == noneq_distribution), \
        "Free massless particle distribution should be conserved in conformal coordinates"


@with_setup(setup)
def homeostasis_test():

    PARAMS.T_initial *= 2
    PARAMS.infer()

    neutrino = Particle(**SMP.neutrino_e)

    energy_density = neutrino.energy_density()
    density = neutrino.density()
    pressure = neutrino.pressure()
    numerator = neutrino.numerator()
    denominator = neutrino.denominator()

    PARAMS.T /= 2
    PARAMS.aT *= 7
    PARAMS.a += 4

    assert neutrino.energy_density() == energy_density \
        and neutrino.density() == density \
        and neutrino.pressure() == pressure \
        and neutrino.numerator() == numerator \
        and neutrino.denominator() == denominator, "Particle state should be persistent"


@with_setup(setup)
def smooth_decoupling_test():

    neutrino = Particle(**SMP.neutrino_e)

    energy_density = REGIMES.RADIATION.energy_density(neutrino)
    density = REGIMES.RADIATION.density(neutrino)
    pressure = REGIMES.RADIATION.pressure(neutrino)
    numerator = REGIMES.RADIATION.numerator(neutrino)
    denominator = REGIMES.RADIATION.denominator(neutrino)
    assert neutrino.energy_density() - energy_density < eps \
        and neutrino.density() - density < eps \
        and neutrino.pressure() - pressure < eps \
        and neutrino.numerator() - numerator < eps \
        and neutrino.denominator() - denominator < eps
