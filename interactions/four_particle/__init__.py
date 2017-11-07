# -*- coding: utf-8 -*-
import numpy
from array import array

import environment
from common.integrators import paired as paired_integrators
from interactions.boltzmann import BoltzmannIntegral
# from interactions.four_particle.integral import integrand
from interactions.four_particle.cpp.integral import integrand, M_t, grid_t, particle_t, reaction_t


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

    creaction = None
    cMs = None

    def __init__(self, **kwargs):
        super(FourParticleIntegral, self).__init__(**kwargs)

    def initialize(self):
        """
        Initialize collision integral constants and save them to the first involved particle
        """
        params = self.particle.params
        if params.T > self.decoupling_temperature and not self.particle.in_equilibrium:
            self.particle.collision_integrals.append(self)

        if self.grids is None:
            self.grids = tuple([self.reaction[1].specie.grid, self.reaction[2].specie.grid])

    def integrate(self, p0, bounds=None):
        if bounds is None:
            bounds = (
                self.grids[0].BOUNDS,
                (lambda p1: self.grids[1].MIN_MOMENTUM, lambda p1: self.grids[1].MAX_MOMENTUM)
            )

        if not self.creaction:
            self.creaction = [
                reaction_t(
                    specie=particle_t(
                        m=particle.specie.conformal_mass,
                        grid=grid_t(
                            grid=particle.specie.grid.TEMPLATE,
                            distribution=particle.specie._distribution,
                            size=particle.specie.grid.MOMENTUM_SAMPLES
                        ),
                        eta=int(particle.specie.eta),
                        in_equilibrium=int(particle.specie.in_equilibrium),
                        aT=particle.specie.aT
                    ),
                    side=particle.side
                )
                for particle in self.reaction
            ]
        if not self.cMs:
            self.cMs = [M_t(list(M.order), M.K1, M.K2) for M in self.Ms]

        def prepared_integrand(p1, p2):
            integrand_1, integrand_f = integrand(p0, p1, p2, self.creaction, self.cMs)
            return numpy.reshape(integrand_1, p1.shape + p2.shape), numpy.reshape(integrand_f, p1.shape + p2.shape)

        params = self.particle.params

        constant = (params.m / params.x)**5 / 64. / numpy.pi**3

        integral_1, integral_f = paired_integrators.integrate_2D(
            prepared_integrand,
            bounds=bounds
        )

        return constant * integral_1, constant * integral_f
