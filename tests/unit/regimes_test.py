from common import Params, UNITS
from particles import Particle, REGIMES
from library.SM import particles as SMP

from . import eps, setup, with_setup_args


@with_setup_args(setup)
def params_inferrence_test(params):
    assert params.T == 5. * UNITS.MeV
    assert params.aT == params.m, "Initial `aT` must be equal to 1"
    assert params.x - params.m * params.aT / params.T < eps


@with_setup_args(setup)
def radiation_regime_test(params):

    photon = Particle(params=params, **SMP.photon)
    assert photon.in_equilibrium, "Photon must always stay in equilibrium"
    assert not photon.decoupling_temperature, "Photon can't decouple"
    assert photon.regime == REGIMES.RADIATION, "Photon is a relativistic particle"
    assert photon.mass == 0, "Photon is massless"
    assert photon.eta == -1, "Photon is a boson, it's eta must be equal to -1"
    assert photon.numerator() == 0, "Photon does not contribute to the numerator"
    assert photon.denominator() != 0, "Photon does contribute to the denominator"


@with_setup_args(setup)
def intermediate_regime_test(params):

    electron = Particle(params=params, **SMP.leptons.electron)
    assert electron.in_equilibrium, "Electron must always stay in equilibrium"
    assert not electron.decoupling_temperature, "Electron can't decouple"
    assert electron.regime == REGIMES.INTERMEDIATE, \
        "Electron is not strictly relativistic at {} MeV".format(params.T/UNITS.MeV)
    assert electron.mass != 0, "Electron is massive"
    assert electron.eta == 1, "Electron is a fermion, it's eta must be equal to 1"
    assert electron.numerator() != 0, "Massive particles contribute to the numerator"
    assert electron.denominator() != 0, "Massive particles contribute to the denominator"


@with_setup_args(setup)
def dust_regime_test(params):

    proton = Particle(params=params, **SMP.hadrons.proton)
    assert proton.in_equilibrium, "Proton must always stay in equilibrium"
    assert not proton.decoupling_temperature, "Proton can't decouple"
    assert proton.regime == REGIMES.DUST, \
        "Proton non-relativistic at {} MeV".format(params.T/UNITS.MeV)
    assert proton.mass != 0, "Proton is massive"
    assert proton.eta == 1, "Proton is a fermion, it's eta must be equal to 1"
    assert proton.numerator() != 0, "Massive particles contribute to the numerator"
    assert proton.denominator() != 0, "Massive particles contribute to the denominator"


@with_setup_args(setup)
def mass_regime_switching_test(params):

    electron = Particle(params=params, **SMP.leptons.electron)
    assert electron.regime == REGIMES.INTERMEDIATE

    electron.mass = 0
    assert electron.regime == REGIMES.RADIATION

    electron.mass = 1e5 * UNITS.MeV
    assert electron.regime == REGIMES.DUST


@with_setup_args(setup)
def temperature_regime_switching_test(params):

    electron = Particle(params=params, **SMP.leptons.electron)
    assert electron.regime == REGIMES.INTERMEDIATE

    params.T = 100 * UNITS.MeV
    electron.update()
    assert electron.regime == REGIMES.RADIATION

    params.T = 1 * UNITS.keV
    electron.update()
    assert electron.regime == REGIMES.DUST


@with_setup_args(setup)
def statistics_consistency_test(params):

    photon = Particle(params=params, **SMP.photon)
    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    assert 7./8. - neutrino.energy_density() / photon.energy_density() < eps


def decoupling_test():

    params = Params(T=SMP.leptons.neutrino_e['decoupling_temperature'] * 2,
                    dy=0.025)
    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    assert neutrino.in_equilibrium
    eq_distribution = neutrino._distribution

    params.T /= 2
    params.a = params.m / params.T
    params.infer()

    assert neutrino.in_equilibrium, "Neutrino should not depend on global temperature change"
    neutrino.update()

    assert not neutrino.in_equilibrium, "Neutrino should have decoupled"
    noneq_distribution = neutrino._distribution

    assert all(eq_distribution == noneq_distribution), \
        "Free massless particle distribution should be conserved in conformal coordinates"


@with_setup_args(setup)
def homeostasis_test(params):

    params.T *= 2
    params.infer()

    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    energy_density = neutrino.energy_density()
    density = neutrino.density()
    pressure = neutrino.pressure()
    numerator = neutrino.numerator()
    denominator = neutrino.denominator()

    params.T /= 2
    params.aT *= 7
    params.a += 4

    assert neutrino.energy_density() == energy_density \
        and neutrino.density() == density \
        and neutrino.pressure() == pressure \
        and neutrino.numerator() == numerator \
        and neutrino.denominator() == denominator, "Particle state should be persistent"


@with_setup_args(setup)
def smooth_decoupling_test(params):

    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

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
