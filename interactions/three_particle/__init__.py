# -*- coding: utf-8 -*-
import numpy

import environment
from interactions.boltzmann import BoltzmannIntegral
from interactions.three_particle.cpp.integral import (
    integration_3, grid_t3, particle_t3, reaction_t3
)
from interactions.four_particle.cpp.integral import CollisionIntegralKind

class ThreeParticleM(object):

    """ ## Three-particle interaction matrix element
        Matrix elements of the interest for three-particle interactions are constant """

    K = 0.

    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        return "|M|² ={: .2e}".format(self.K)

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

    def __init__(self, **kwargs):
        super(ThreeParticleIntegral, self).__init__(**kwargs)

    def initialize(self):
        """
        Initialize collision integral constants and save them to the first involved particle
        """
        self.par = self.particle.params
        if self.par.T > self.washout_temperature and not self.particle.in_equilibrium:
            self.particle.collision_integrals.append(self)

        if self.grids is None:
            self.grids = self.reaction[1].specie.grid

    def integrate(self, ps, stepsize=None, bounds=None):
        params = self.particle.params

        if self.reaction[0].specie.mass == 0 and self.reaction[1].side == 1:
            return numpy.zeros(len(ps))

        bounds = tuple(b1 / params.aT for b1 in self.grids.BOUNDS) + (self.reaction[2].specie.grid.MAX_MOMENTUM / params.aT, )

        if stepsize is None:
            stepsize = params.h

        if not environment.get('LOGARITHMIC_TIMESTEP'):
            stepsize /= params.aT

        self.creaction = [
            reaction_t3(
                specie=particle_t3(
                    m=particle.specie.conformal_mass / params.aT,
                    grid=grid_t3(
                        grid=particle.specie.grid.TEMPLATE / params.aT,
                        distribution=particle.specie._distribution
                    ),
                    eta=int(particle.specie.eta),
                    in_equilibrium=int(particle.specie.in_equilibrium),
                    T=particle.specie.aT / params.aT
                ),
                side=particle.side
            )
            for particle in self.reaction
        ]

        ps = ps / params.aT
        self.cMs = sum(M.K for M in self.Ms)

        if self.particle.mass == 0:
            constant_0 = 0
        else:
            constant_0 = self.cMs * params.aT / params.a / 8. / numpy.pi / params.H / self.particle.mass**2

        constant_else = self.cMs * params.a / params.aT / 32. / numpy.pi / params.H

        constant = numpy.append(constant_0, numpy.repeat(constant_else, len(ps) - 1)) / self.particle.dof

        if not environment.get('LOGARITHMIC_TIMESTEP'):
            constant /= params.x

        stepsize *= constant_else

        if environment.get('SPLIT_COLLISION_INTEGRAL'):
            A = integration_3(ps, *bounds, self.creaction, stepsize, CollisionIntegralKind.F_1)
            B = integration_3(ps, *bounds, self.creaction, stepsize, CollisionIntegralKind.F_f)
            return numpy.array(A) * constant, numpy.array(B) * constant

        fullstack = integration_3(ps, *bounds, self.creaction, stepsize, self.kind)
        if self.kind == CollisionIntegralKind.F_f or self.kind == CollisionIntegralKind.F_f_vacuum_decay:
            constant *= self.particle._distribution

        return numpy.array(fullstack) * constant