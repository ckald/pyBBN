"""
General non-equilibrium particles
"""
import numpy
from scipy import interpolate
from common import lambda_integrate, PARAMS, momentum_to_index, GRID


name = 'non-equilibrium'


@lambda_integrate
def density(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 2
        * particle.dof / 2. / numpy.pi**2 / PARAMS.a**3
    ), otypes=[numpy.float_])


@lambda_integrate
def energy_density(particle):
    """ == Energy density ==

        \begin{equation}
            \rho = \frac{g}{2 \pi^2} \frac{m^4}{x^4} \int dy y^2 \sqrt{y^2 +\
            \frac{M_N^2 x^2}{m^2}} f(y)
        \end{equation}
    """
    return numpy.vectorize(lambda y: (
        particle.distribution(y) * y ** 2
        * particle.energy_normalized(y)
        * particle.dof / 2. / numpy.pi**2 / PARAMS.a**4
    ), otypes=[numpy.float_])


@lambda_integrate
def pressure(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 4
        / particle.energy_normalized(p)
        * particle.dof / 6. / numpy.pi**2 / PARAMS.a**4
    ), otypes=[numpy.float_])


@lambda_integrate
def numerator(particle):
    integral = interpolate.interp1d(GRID.TEMPLATE, particle.collision_integral,
                                    kind='quadratic', assume_sorted=True, copy=False)
    func = numpy.vectorize(lambda y: (
        -1. * particle.dof / 2. / numpy.pi**2
        * y**2 * particle.energy_normalized(y)
        * integral(y)
    ), otypes=[numpy.float_])
    return func


def denominator(particle):
    return 0.
