"""
# Non-equilibrium particles
"""
from __future__ import division
import numpy
from common import linear_interpolation
from common.integrators import lambda_integrate


name = 'non-equilibrium'


@lambda_integrate()
def density(particle):
    return numpy.vectorize(lambda y: (
        particle.distribution(y) * y**2
        * particle.dof / 2. / numpy.pi**2 / particle.params.a**3
    ), otypes=[numpy.float_])


@lambda_integrate()
def energy_density(particle):
    """ ### Energy density

        \begin{equation}
            \rho = \frac{g}{2 \pi^2} \frac{m^4}{x^4} \int dy y^2 \sqrt{y^2 +\
            \frac{M_N^2 x^2}{m^2}} f(y)
        \end{equation}
    """
    return numpy.vectorize(lambda y: (
        particle.distribution(y)
        * y**2 * particle.conformal_energy(y)
        * particle.dof / 2. / numpy.pi**2 / particle.params.a**4
    ), otypes=[numpy.float_])


@lambda_integrate()
def pressure(particle):
    """ ### Pressure

        \begin{equation}
            p = \frac{g}{6 \pi^2} \frac{m^4}{x^4} \int \frac{dy \, y^4 f(y)}\
            { \sqrt{y^2 + \frac{M_N^2 x^2}{m^2}} }
        \end{equation}
    """

    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p**4 / particle.conformal_energy(p)
        * particle.dof / 6. / numpy.pi**2 / particle.params.a**4
    ), otypes=[numpy.float_])


@lambda_integrate()
def entropy(particle):
    """ ## Entropy

        \begin{equation}
            s = - \int_0^\inf p^2 dp \left{ f(p) \ln f(p) \mp (1 \pm f(p)) \ln (1 \pm f(p)) \right}
        \end{equation}
    """

    def integrand(p):
        f = particle.distribution(p)
        eta = particle.eta

        if f == 0:
            return 0.

        return (- particle.dof / 2 / numpy.pi**2 / particle.params.a**3
                * p**2 * (f * numpy.log(f) + eta * (1 - eta * f) * numpy.log(1 - eta * f)))

    return numpy.vectorize(integrand, otypes=[numpy.float_])


""" ## Master equation terms """

""" ### Numerator

    \begin{equation}
        -\frac{g}{2 \pi^2} \int dy \, y^2 \sqrt{y^2 + \frac{M_N^2 x^2}{m^2}} I_{coll}(y)
    \end{equation}
"""


@lambda_integrate()
def numerator(particle):
    integral = linear_interpolation(particle.collision_integral / particle.params.x,
                                    particle.grid.TEMPLATE)
    return numpy.vectorize(lambda y: (
        -1. * particle.dof / 2. / numpy.pi**2
        * y**2 * particle.conformal_energy(y)
        * integral(y)
    ), otypes=[numpy.float_])


def denominator(particle):
    """
    ### Denominator

    \begin{equation}
        0
    \end{equation}
    """
    return 0.
