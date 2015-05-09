# -*- coding: utf-8 -*-

import numpy
from collections import namedtuple
from common import GRID, CONST
from common.integrators import integrate_1D
from library import SM


q = (SM.particles.hadrons.neutron['mass'] - SM.particles.hadrons.proton['mass'])
m_e = SM.particles.leptons.electron['mass']

a = None

Particles = namedtuple("Particles", "electron neutrino")
particles = None


default_bounds = (GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,)


def init_kawano(electron=None, neutrino=None):
    global particles
    particles = Particles(electron=electron, neutrino=neutrino)


def _rate1(y):
    """ n + ν_e ⟶  e + p """
    E_e = q*a + y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * (1. - particles.electron.distribution(y_e)) * particles.neutrino.distribution(y))


def _rate1b(y):
    """ e + p ⟶  n + ν_e """
    E_e = q*a + y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * particles.electron.distribution(y_e) * (1. - particles.neutrino.distribution(y)))


def _rate2(y):
    """ n ⟶  e + ν_e' + p """
    E_e = q*a - y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * (1. - particles.electron.distribution(y_e))
            * (1. - particles.neutrino.distribution(y)))


def _rate2b(y):
    """ e + ν_e' + p ⟶  n """
    E_e = q*a - y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * particles.electron.distribution(y_e) * particles.neutrino.distribution(y))


def _rate3(y):
    """ n + e' ⟶  ν_e' + p """
    E_e = -q*a + y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y_e**2 * y_e * E_e
            * particles.electron.distribution(y_e) * (1. - particles.neutrino.distribution(y)))


def _rate3b(y):
    """ ν_e' + p ⟶  n + e' """
    E_e = numpy.sqrt(y**2 + (m_e*a)**2)
    y_n = q*a + E_e

    return (y**2 * y_n**2
            * (1. - particles.electron.distribution(y)) * particles.neutrino.distribution(y_n))


def baryonic_rates(_a):
    global a
    a = _a

    return [
        CONST.rate_normalization / a**5 * integrate_1D(numpy.vectorize(integrand), bounds=bounds)[0]
        if bounds[0] < bounds[1] else 0.
        for integrand, bounds in [
            (_rate1, default_bounds),
            (_rate1b, default_bounds),
            (_rate2, (GRID.MIN_MOMENTUM, q*a - (m_e*a))),
            (_rate2b, (GRID.MIN_MOMENTUM, q*a - (m_e*a))),
            (_rate3, (q*a + (m_e*a), GRID.MAX_MOMENTUM)),
            (_rate3b, default_bounds)]
    ]
