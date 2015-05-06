import math
from collections import namedtuple
from common import GRID, CONST
from common.integrators import integrate_1D


from library.SM.particles.hadrons import neutron, proton
from library.SM.particles.leptons import electron
q = (neutron['mass'] - proton['mass'])
m_e = electron['mass']

a = None

Particles = namedtuple("Particles", "electron neutrino neutron proton")
particles = None


default_bounds = (GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,)


def init_kawano(electron=None, neutrino=None, neutron=None, proton=None):
    global particles
    particles = Particles(electron=electron, neutrino=neutrino, neutron=neutron, proton=proton)


def _rate1(y):
    """ n + ν_e ⟶  e + p """
    E_e = q*a + y
    y_e = math.sqrt(E_e**2 - (m_e*a)**2)
    return (
        (y**2 * y_e * E_e
         * (1. - particles.electron.distribution(y_e)) * particles.neutrino.distribution(y)),
        default_bounds
    )


def _rate1b(y):
    """ e + p ⟶  n + ν_e """
    E_e = q*a + y
    y_e = math.sqrt(E_e**2 - (m_e*a)**2)
    return (
        (y**2 * y_e * E_e
         * particles.electron.distribution(y_e) * (1. - particles.neutrino.distribution(y))),
        default_bounds
    )


def _rate2(y):
    """ n ⟶  e + ν_e' + p """
    E_e = q*a - y
    y_e = math.sqrt(E_e**2 - (m_e*a)**2)
    return (
        (y**2 * y_e * E_e
         * (1. - particles.electron.distribution(y_e))
         * (1. - particles.neutrino.distribution(y))),
        (GRID.MIN_MOMENTUM, q*a - (m_e*a))
    )


def _rate2b(y):
    """ e + ν_e' + p ⟶  n """
    E_e = q*a - y
    y_e = math.sqrt(E_e**2 - (m_e*a)**2)
    return (
        (y**2 * y_e * E_e
         * particles.electron.distribution(y_e) * particles.neutrino.distribution(y)),
        (GRID.MIN_MOMENTUM, q*a - (m_e*a))
    )


def _rate3(y):
    """ n + e' ⟶  ν_e' + p """
    E_e = -q*a + y
    y_e = math.sqrt(E_e**2 - (m_e*a)**2)
    return (
        (y_e**2 * y_e * E_e
         * particles.electron.distribution(y_e) * (1. - particles.neutrino.distribution(y))),
        default_bounds
    )


def _rate3b(y):
    """ ν_e' + p ⟶  n + e' """
    E_e = math.sqrt(y**2 + (m_e*a)**2)
    y_n = q*a + E_e

    return (
        (y**2 * y_n**2
         * (1. - particles.electron.distribution(y)) * particles.neutrino.distribution(y_n)),
        default_bounds
    )


def baryonic_rates(_a):
    global a
    a = _a

    return [
        CONST.rate_normalization / a**5 * integrate_1D(rate[0], bounds=rate[1])
        for rate in [_rate1, _rate1b, _rate2, _rate2b, _rate3, _rate3b]
    ]
