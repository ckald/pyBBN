# -*- coding: utf-8 -*-
import numpy
from common import integrators
from interactions.boltzmann import BoltzmannIntegral
from interactions.ds import D, Db


class FourParticleM(object):

    """
    ## Matrix element
    All four-particle interactions of the interest can be rewritten in a form

    \begin{equation}
        |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
             + K_2 m_i m_j (p_k \cdot p_l)
    \end{equation}
    """
    K1 = 0.
    K2 = 0.
    # Order defines the values of the $(i, j, k, l)$ indices
    order = (0, 1, 2, 3)

    def __init__(self, **kwargs):
        """ Configure matrix element, check that it makes sense """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        if len(set(self.order)) != len(self.order):
            raise Exception("Meaningless order of momenta of the matrix element: {}"
                            .format(self.order))

    def __iadd__(self, M):
        self.K1 += M.K1
        self.K2 += M.K2
        return self

    def __idiv__(self, div):
        self.K1 /= div
        self.K2 /= div
        return self

    def __imul__(self, mul):
        self.K1 *= mul
        self.K2 *= mul
        return self

    def __setattr__(self, name, value):
        if name != 'order':
            return super(FourParticleM, self).__setattr__(name, value)

        self.__dict__['order'] = tuple(sorted(value[:2]) + sorted(value[2:]))

    def __str__(self):
        """ String-like representation of the matrix element """
        ret = ""
        if self.K1:
            ret += "K1={: .2e} ".format(self.K1)
        if self.K2:
            ret += "K2={: .2e} ".format(self.K2)
        return ret + "{}".format(self.order)

    def apply_order(self, order, reaction):
        self.order = order

        # Change the sign of the mass term if related particle has crossed (but not both)
        if bool(reaction[order[0]].crossed) ^ bool(reaction[order[1]].crossed):
            self.K2 *= -1.

    def stackable(self, other):
        return self.order == other.order or \
            (tuple(self.order[2:] + self.order[:2]) == other.order and self.K2 == other.K2 == 0)


class FourParticleIntegral(BoltzmannIntegral):

    def __init__(self, **kwargs):
        super(FourParticleIntegral, self).__init__(**kwargs)

        if self.grids is None:
            self.grids = tuple([self.reaction[1].specie.grid, self.reaction[2].specie.grid])

    def initialize(self):
        """
        Initialize collision integral constants and save them to the first involved particle
        """
        params = self.particle.params
        if params.T > self.decoupling_temperature and not self.particle.in_equilibrium:
            self.particle.collision_integrals.append(self)

    def integrate(self, p0, fau=None, bounds=None):
        if bounds is None:
            bounds = (
                self.grids[0].BOUNDS,
                (lambda p1: self.grids[1].MIN_MOMENTUM,
                 lambda p1: min(p0 + p1, self.grids[1].MAX_MOMENTUM)),
            )

        def prepared_integrand(p1, p2):
            return self.integrand(p0, p1, p2, fau)

        params = self.particle.params
        constant = (params.m / params.x)**5 / 64. / numpy.pi**3
        integral, error = integrators.integrate_2D(
            prepared_integrand,
            bounds=bounds
        )

        return constant * integral

    def integrand(self, p0, p1, p2, fau=None):

        """
        Collision integral interior.
        """

        p = [p0, p1, p2, 0]
        p, E, m = self.calculate_kinematics(p)

        if not self.in_bounds(p, E, m):
            return 0.

        integrand = 1.

        ds = 0.
        if p[0] != 0:
            for M in self.Ms:
                ds += D(p=p, E=E, m=m, M=M, reaction=self.reaction)
            ds = ds / p[0] / E[0]
        else:
            for M in self.Ms:
                ds += Db(p=p, E=E, m=m, M=M, reaction=self.reaction)
        integrand *= ds

        # Avoid rounding errors and division by zero
        for i in [1, 2, 3]:
            if m[i] != 0:
                integrand *= p[i] / E[i]

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

        q1, q2 = (p[0], p[1]) if p[0] > p[1] else (p[1], p[0])
        q3, q4 = (p[2], p[3]) if p[2] > p[3] else (p[3], p[2])

        is_in = (E[3] >= m[3]
                 and q1 <= q2 + q3 + q4
                 and q3 <= q1 + q2 + q4)

        return is_in
