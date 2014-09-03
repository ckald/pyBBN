"""
Ultra-relativistic simplifications of density, energy density and pressure calculations
"""
import numpy
from common import STATISTICS, PARAMS


def density(particle):
    density = PARAMS.aT**3 * particle.dof / 2. / numpy.pi**2 * 2. * 1.2
    if particle.statistics == STATISTICS.FERMION:
        density *= 3./4.
    return density


def energy_density(particle):
    density = PARAMS.T**4 * particle.dof * numpy.pi**2 / 30.
    if particle.statistics == STATISTICS.FERMION:
        density *= 7./8.
    return density


def pressure(particle):
    return 1. * energy_density(particle) / 3.


def numerator(particle):
    return 0.


def denominator(particle):
    return 2. * numpy.pi**2 / 15. * particle.dof * PARAMS.aT**3
