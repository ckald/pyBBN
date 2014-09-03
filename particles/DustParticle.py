"""
Non-relativistic simplifications of density, energy density and pressure calculations
"""
import numpy
import IntermediateParticle


def density(particle):
    return (
        particle.dof
        * numpy.sqrt(particle.mass * particle.temperature / 2. / numpy.pi)**3
        * numpy.exp(- particle.mass / particle.temperature)
    )


def energy_density(particle):
    return (particle.mass + 3./2. * particle.temperature) * density(particle)


def pressure(particle):
    return density(particle) * particle.temperature


numerator = IntermediateParticle.numerator
denominator = IntermediateParticle.denominator
