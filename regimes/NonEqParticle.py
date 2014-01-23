"""
General non-equilibrium particles
"""
import numpy
from common import lambda_integrate


@lambda_integrate
def density(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 2
        * particle.dof / 2. / numpy.pi**2
    ))


@lambda_integrate
def energy_density(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 2
        * particle.energy(p)
        * particle.dof / 2. / numpy.pi**2
    ))


@lambda_integrate
def pressure(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 4
        / particle.energy(p)
        * particle.dof / 6. / numpy.pi**2
    ))


@lambda_integrate
def energy_density_rate(particle):
    return numpy.vectorize(lambda p: 0.)
    # return lambda p: (
    #     particle.dof / 2. / numpy.pi**2 / particle.temperature
    #     * p**2 * particle.energy(p)
    #     * numpy.exp(particle.energy(p) / particle.temperature)
    #     / (
    #         numpy.exp(particle.energy(p) / particle.temperature)
    #         + particle.eta
    #     )**2
    # )

# @lambda_integrate
# def density(particle):
#     return lambda i: (
#         particle.distribution(i, by_index=True) * index_to_momentum(i) ** 2
#         * particle.dof / 2. / numpy.pi**2
#     )


# @lambda_integrate
# def energy_density(particle):
#     return lambda i: (
#         particle.distribution(i, by_index=True) * index_to_momentum(i) ** 2
#         * particle.energy(index_to_momentum(i))
#         * particle.dof / 2 / numpy.pi**2
#     )


# @lambda_integrate
# def pressure(particle):
#     return lambda i: (
#         particle.distribution(i, by_index=True) * index_to_momentum(i) ** 4
#         / particle.energy(index_to_momentum(i))
#         * particle.dof / 6 / numpy.pi**2
#     )


# @lambda_integrate
# def energy_density_rate(particle):
#     return lambda i: (
#         particle.dof / 2. / numpy.pi**2 / particle.temperature
#         * index_to_momentum(i)**2 * particle.energy(index_to_momentum(i))
#         * numpy.exp(particle.energy(index_to_momentum(i)) / particle.temperature)
#         / (
#             numpy.exp(particle.energy(index_to_momentum(i)) / particle.temperature)
#             + particle.eta
#         )**2
#     )
