"""
= Non-equilibrium particles =
"""
from __future__ import division
import functools
import numpy
from scipy import interpolate, integrate
from common import PARAMS, GRID


name = 'non-equilibrium'


def lambda_integrate(func):
    """ Scipy integration over the momentum space of the lambda function applied to the \
        grid template """

    method = ['simps', 'quad', 'romberg'][1]

    @functools.wraps(func)
    def wrapper(*args, **kw):
        if method == 'simps':
            fpp = func(*args, **kw)(GRID.TEMPLATE)
            result = integrate.simps(fpp, dx=GRID.MOMENTUM_STEP)

        elif method == 'quad':
            fpp = func(*args, **kw)
            result, err = integrate.quad(fpp, GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM)

        elif method == 'romberg':
            fpp = func(*args, **kw)
            result = integrate.romberg(fpp, GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM)

        return result
    return wrapper


@lambda_integrate
def density(particle):
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 2
        * particle.dof / 2. / numpy.pi**2 / PARAMS.a**3
    ), otypes=[numpy.float_])


@lambda_integrate
def energy_density(particle):
    """ === Energy density ===

        \begin{equation}
            \rho = \frac{g}{2 \pi^2} \frac{m^4}{x^4} \int dy y^2 \sqrt{y^2 +\
            \frac{M_N^2 x^2}{m^2}} f(y)
        \end{equation}
    """
    return numpy.vectorize(lambda y: (
        particle.distribution(y) * y ** 2
        * particle.energy_normalized(y)
        * particle.dof / 2. / numpy.pi**2 / PARAMS.a**4
    ))


@lambda_integrate
def pressure(particle):
    """ === Pressure ===

        \begin{equation}
            p = \frac{g}{6 \pi^2} \frac{m^4}{x^4} \int \frac{dy \, y^4 f(y)}\
            { \sqrt{y^2 + \frac{M_N^2 x^2}{m^2}} }
        \end{equation}
    """
    return numpy.vectorize(lambda p: (
        particle.distribution(p) * p ** 4
        / particle.energy_normalized(p)
        * particle.dof / 6. / numpy.pi**2 / PARAMS.a**4
    ))


""" == Master equation terms == """

""" === Numerator ===

    \begin{equation}
        -\frac{g}{2 \pi^2} \int dy \, y^2 \sqrt{y^2 + \frac{M_N^2 x^2}{m^2}} I_{coll}(y)
    \end{equation}
"""


@numpy.vectorize
def numerator_lambda(y, particle, integral):
    return (
        -1. * particle.dof / 2. / numpy.pi**2
        * y**2 * particle.energy_normalized(y)
        * integral(y)
    )


@lambda_integrate
def numerator(particle):
    integral = interpolate.interp1d(GRID.TEMPLATE, particle.collision_integral,
                                    kind='quadratic', assume_sorted=True, copy=False)
    return functools.partial(numerator_lambda, particle=particle, integral=integral)


def denominator(particle):
    """
    === Denominator ===

    \begin{equation}
        0
    \end{equation}
    """
    return 0.
