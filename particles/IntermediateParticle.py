"""
For intermediate regime equilibrium particles, density, energy density and pressure
are obtained through integration of distribution function
"""
from scipy import integrate
import numpy
from common import GRID, benchmark, PARAMS


def density(particle):
    if not particle._density:
        particle._density = integrate.quad(
            lambda p: (
                particle.distribution_function(
                    particle.energy(p) / particle.temperature
                ) * p**2 * particle.dof / 2. / numpy.pi**2
            ), GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
        )[0]
    return particle._density


def energy_density_integrand(p, particle):
    E = particle.energy(p)
    return (
        particle.distribution_function(E / particle.temperature)
        * p**2 * E
        * particle.dof / 2. / numpy.pi**2
    )


def energy_density(particle):
    if not particle._energy_density:
        particle._energy_density = integrate.quad(
            lambda p: energy_density_integrand(p, particle),
            GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
        )[0]
    return particle._energy_density


def pressure_integrand(p, particle):
    E = particle.energy(p)
    return (
        particle.distribution_function(E / particle.temperature)
        * p**4 / E
        * particle.dof / 6. / numpy.pi**2
    )


def pressure(particle):
    if not particle._pressure:
        particle._pressure = integrate.quad(
            lambda p: pressure_integrand(p, particle),
            GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
        )[0]
    return particle._pressure


def numerator(particle):
    return particle.mass**2 * PARAMS.x / PARAMS.m**2 / PARAMS.aT * I(particle)


def denominator(particle):
    return (J(particle) + particle.mass_normalized**2 * I(particle)) / PARAMS.aT**2


def I(particle):
    return particle.dof / 2 / numpy.pi**2 * integrate.quad(
        lambda y: (
            y**2 * numpy.exp(-particle.energy_normalized(y) / PARAMS.aT)
            / (numpy.exp(-particle.energy_normalized(y) / PARAMS.aT) + particle.eta) ** 2
        ), GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
    )[0]


def J(particle):
    return particle.dof / 2 / numpy.pi**2 * integrate.quad(
        lambda y: (
            y**4 * numpy.exp(-particle.energy_normalized(y) / PARAMS.aT)
            / (numpy.exp(-particle.energy_normalized(y) / PARAMS.aT) + particle.eta) ** 2
        ), GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM
    )[0]
