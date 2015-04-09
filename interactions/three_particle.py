# -*- coding: utf-8 -*-
import numpy
from common import GRID, integrators
from interactions.boltzmann import BoltzmannIntegral


class ThreeParticleM(object):

    """ ## Three-particle interaction matrix element
        Matrix elements of the interest for three-particle interactions are constant """

    K = 0.

    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        return "K={: .2e}".format(self.K)

    def __iadd__(self, M):
        self.K += M.K
        return self

    def __idiv__(self, div):
        self.K /= div
        return self

    def __imul__(self, mul):
        self.K *= mul
        return self


class ThreeParticleIntegral(BoltzmannIntegral):

    def initialize(self):
        """
        Initialize collision integral constants and save them to the first involved particle
        """
        params = self.in_particles[0].params
        if params.T > self.decoupling_temperature and not self.in_particles[0].in_equilibrium:
            MM = 0
            for M in self.Ms:
                MM += M.K
            self.constant = MM / 16. / numpy.pi * params.m / params.x / params.H
            self.particles[0].collision_integrals.append(self)

    @staticmethod
    def integrate(p0, integrand, bounds=None, kwargs=None):
        kwargs = kwargs if kwargs else {}

        if bounds is None:
            bounds = (GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM)

        if isinstance(integrand, list):
            def prepared_integrand(p1):
                return sum([i(p0, p1, **kwargs) for i in integrand])
        else:
            def prepared_integrand(p1):
                return integrand(p0, p1, **kwargs)

        integral, error = integrators.integrate_1D(
            prepared_integrand,
            bounds=bounds
        )

        return integral, error

    def integrand(self, p0, p1, fau=None):

        """
        Collision integral interior.
        """

        p = [p0, p1, 0]
        p, E, m = self.calculate_kinematics(p)

        if not self.in_bounds(p, E, m):
            return 0.

        integrand = self.constant / p[0] / E[0]

        # Avoid rounding errors and division by zero
        if m[1] != 0:
            integrand *= p[1] / E[1]

        if integrand == 0:
            return 0

        integrand *= fau(p)

        return integrand

    """ ### Integration region bounds methods """

    def in_bounds(self, p, E=None, m=None):
        """ $D$-functions involved in the interactions imply a cut-off region for the collision\
            integrand. In the general case of arbitrary particle masses, this is a set of \
            irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
            not succeed in this). To avoid excessive computations, it is convenient to do an early\
            `return 0` when the particles kinematics lay out of the cut-off region """
        if not E or not m:
            p, E, m = self.calculate_kinematics(p)

        is_in = E[2] >= m[2] \
            and p[0] + p[1] > p[2] \
            and p[0] + p[2] > p[1] \
            and p[1] + p[2] > p[0]

        return is_in
