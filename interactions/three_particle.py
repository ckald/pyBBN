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
        params = self.particle.params
        if params.T > self.decoupling_temperature and not self.particle.in_equilibrium:
            MM = 0
            for M in self.Ms:
                MM += M.K
            self.constant = MM / 16. / numpy.pi * params.m / params.x / params.H
            self.particle.collision_integrals.append(self)

    @staticmethod
    def integrate(E0, integrand, bounds=None, kwargs=None):
        kwargs = kwargs if kwargs else {}

        if bounds is None:
            bounds = (GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM)

        if isinstance(integrand, list):
            def prepared_integrand(E1):
                return sum([i(E0, E1, **kwargs) for i in integrand])
        else:
            def prepared_integrand(E1):
                return integrand(E0, E1, **kwargs)

        integral, error = integrators.integrate_1D(
            prepared_integrand,
            bounds=bounds
        )

        return integral, error

    def integrand(self, E0, E1, fau=None):

        """
        Collision integral interior.
        """

        E = [E0, E1, 0]
        E, p, m = self.calculate_kinematics(E)

        integrand = self.in_bounds(p, E, m) * self.constant / p[0] / E[0]

        integrand *= fau(E)

        return integrand

    """ ### Integration region bounds methods """

    def in_bounds(self, E, p=None, m=None):
        """ The kinematically allowed region in momentum space """
        if not p or not m:
            E, p, m = self.calculate_kinematics(E)

        is_in = (
            (E[2] >= m[2])
            * (p[0] + p[1] > p[2])
            * (p[0] + p[2] > p[1])
            * (p[1] + p[2] > p[0])
        )

        return is_in
