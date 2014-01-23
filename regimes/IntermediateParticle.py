"""
For intermediate regime equilibrium particles, density, energy density and pressure
are obtained through integration of distribution function
"""
from scipy import integrate
import numpy
from common import GRID


def density(particle):
    return integrate.quad(
        lambda p: (
            particle.distribution_function(
                particle.energy(p) / particle.temperature
            )
            * p**2 * particle.dof / 2. / numpy.pi**2
        ), GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
    )[0]


def energy_density(particle):
    def integrand(p, particle):
        E = particle.energy(p)
        return (
            particle.distribution_function(
                E / particle.temperature
            )
            * p**2 * E
            * particle.dof / 2. / numpy.pi**2
        )

    return integrate.quad(
        lambda p: integrand(p, particle),
        GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
    )[0]


def pressure(particle):
    def integrand(p, particle):
        E = particle.energy(p)
        return (
            particle.distribution_function(
                E / particle.temperature
            )
            * p**4 / E
            * particle.dof / 6. / numpy.pi**2
        )
    return integrate.quad(
        lambda p: integrand(p, particle),
        GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
    )[0]


def energy_density_rate(particle):
    def integrand(p, particle):
        E = particle.energy(p)
        return (
            particle.dof / 2. / numpy.pi**2
            / (
                numpy.exp(E / particle.temperature)
                + particle.eta
            )**2
            * numpy.exp(E / particle.temperature)
            * p**2 * E**2 / particle.temperature
        )
    return integrate.quad(
        lambda p: integrand(p, particle),
        GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
    )[0]
