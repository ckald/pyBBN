"""
General non-equilibrium particles
"""
import numpy
from common import lambda_integrate, PARAMS, momentum_to_index


name = 'non-equilibrium'


@lambda_integrate
def density(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 2
        * particle.dof / 2. / numpy.pi**2 / PARAMS.a**3
    ), otypes=[numpy.float_])


@lambda_integrate
def energy_density(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 2
        * particle.energy(p)
        * particle.dof / 2. / numpy.pi**2 / PARAMS.a**4
    ), otypes=[numpy.float_])


@lambda_integrate
def pressure(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 4
        / particle.energy(p)
        * particle.dof / 6. / numpy.pi**2 / PARAMS.a**4
    ), otypes=[numpy.float_])


@lambda_integrate
def numerator(particle):
    return numpy.vectorize(lambda y: (
        -1. * particle.dof / 2. / numpy.pi**2
        * y**2 * particle.energy_normalized(y)
        * particle.collision_integral[momentum_to_index(y)]
    ), otypes=[numpy.float_])


def denominator(particle):
    return 0.
