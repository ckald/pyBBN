"""
# Non-equilibrium particles
"""
from __future__ import division
import numpy
from scipy.integrate import simps

import environment
from common import linear_interpolation
from common.integrators import lambda_integrate


name = 'non-equilibrium'

if not environment.get('SIMPSONS_NONEQ_PARTICLES'):

    def density(particle):
        @lambda_integrate(particle.grid.BOUNDS)
        def den():
            return numpy.vectorize(lambda y: (
                particle.distribution(y) * y**2
                * particle.dof / 2. / numpy.pi**2 / particle.params.a**3
            ), otypes=[numpy.float_])
        return den()


    def energy_density(particle):
        @lambda_integrate(particle.grid.BOUNDS)
        def energy_den():
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
        return energy_den()


    def pressure(particle):
        @lambda_integrate(particle.grid.BOUNDS)
        def press():
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
        return press()


    def entropy(particle):
        @lambda_integrate(particle.grid.BOUNDS)
        def ent():
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
        return ent()


    # @lambda_integrate()
    # def inverse_gamma_factor(particle):
    #     """ ## Average gamma factor of the distribution """

    #     def integrand(p):
    #         if particle.mass == 0:
    #             return 0.

    #         return (particle.dof / 2 / numpy.pi**2 * p**2 / particle.params.a**3
    #                 * particle.distribution(p) * particle.conformal_mass / particle.conformal_energy(p))

    #     return numpy.vectorize(integrand, otypes=[numpy.float_])


    """ ## Master equation terms """

    """ ### Numerator

        \begin{equation}
            -\frac{g}{2 \pi^2} \int dy \, y^2 \sqrt{y^2 + \frac{M_N^2 x^2}{m^2}} I_{coll}(y)
        \end{equation}
    """


    def numerator(particle):
        integral = linear_interpolation(particle.collision_integral / particle.params.x,
                                        particle.grid.TEMPLATE)
        return lambda_integrate(particle.grid.BOUNDS)(lambda particle: numpy.vectorize(lambda y: (
            -1. * particle.dof / 2. / numpy.pi**2
            * y**2 * particle.conformal_energy(y)
            * integral(y)
        ), otypes=[numpy.float_]))(particle)


    def denominator(particle):
        """
        ### Denominator

        \begin{equation}
            0
        \end{equation}
        """
        return 0.



else:

    def density(particle):
        temp = particle.grid.TEMPLATE
        return simps((
            particle.distribution(temp) * temp**2
            * particle.dof / 2. / numpy.pi**2 / particle.params.a**3), temp
        )


    def energy_density(particle):
        """ ### Energy density

            \begin{equation}
                \rho = \frac{g}{2 \pi^2} \frac{m^4}{x^4} \int dy y^2 \sqrt{y^2 +\
                \frac{M_N^2 x^2}{m^2}} f(y)
            \end{equation}
        """
        temp = particle.grid.TEMPLATE
        return simps((
            particle.distribution(temp) * temp**2 * particle.conformal_energy(temp)
            * particle.dof / 2. / numpy.pi**2 / particle.params.a**4), temp
        )


    def pressure(particle):
        """ ### Pressure

            \begin{equation}
                p = \frac{g}{6 \pi^2} \frac{m^4}{x^4} \int \frac{dy \, y^4 f(y)}\
                { \sqrt{y^2 + \frac{M_N^2 x^2}{m^2}} }
            \end{equation}
        """
        temp = particle.grid.TEMPLATE
        return simps((
            particle.distribution(temp) * temp**4 / particle.conformal_energy(temp)
            * particle.dof / 6. / numpy.pi**2 / particle.params.a**4), temp
        )


    def entropy(particle):
        """ ## Entropy

            \begin{equation}
                s = - \int_0^\inf p^2 dp \left{ f(p) \ln f(p) \mp (1 \pm f(p)) \ln (1 \pm f(p)) \right}
            \end{equation}
        """
        temp = particle.grid.TEMPLATE
        eta = particle.eta
        integrand = numpy.zeros(len(temp))

        for num, mom in enumerate(temp):
            f = particle.distribution(mom)

            if f == 0:
                integrand[num] = 0.
            else:
                integrand[num] = (- particle.dof / 2 / numpy.pi**2 / particle.params.a**3
                                * mom**2 * (f * numpy.log(f) + eta * (1 - eta * f) * numpy.log(1 - eta * f)))

        return simps(integrand, temp)


    def numerator(particle):
        temp = particle.grid.TEMPLATE
        return simps((
            -1. * particle.dof / 2. / numpy.pi**2
            * temp**2 * particle.conformal_energy(temp)
            * particle.collision_integral / particle.params.x), temp
        )


    def denominator(particle):
        """
        ### Denominator

        \begin{equation}
            0
        \end{equation}
        """
        return 0.
