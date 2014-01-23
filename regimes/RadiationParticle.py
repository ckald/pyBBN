"""
Ultra-relativistic simplifications of density, energy density and pressure calculations
"""
import numpy
from common import STATISTICS


def density(particle):
    density = particle.temperature**3 * particle.dof / 2. / numpy.pi**2 * 2. * 1.2
    if particle.statistics == STATISTICS.FERMION:
        density *= 3./4.
    return density


def energy_density(particle):
    density = particle.temperature**4 * particle.dof * numpy.pi**2 / 30.
    if particle.statistics == STATISTICS.FERMION:
        density *= 7./8.
    return density


def pressure(particle):
    return 1. * energy_density(particle) / 3.


def energy_density_rate(particle):
    return 4. * energy_density(particle)
