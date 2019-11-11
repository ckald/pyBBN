# -*- coding: utf-8 -*-
import numpy
from scipy.integrate import simps
from collections import Counter

import environment
from common import kinematics, UNITS
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
        if self.particle.params.T < self.particle.decoupling_temperature and not self.particle.in_equilibrium:
            self.particle.collision_integrals.append(self)

        if self.grids is None:
            self.grids = self.reaction[1].specie.grid

    def integrate(self, ps, stepsize=None, bounds=None):
        params = self.particle.params

        if kinematics.Neglect3pInteraction(self, ps):
            return kinematics.return_function(self, ps)

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

        constant = numpy.append(constant_0, numpy.repeat(constant_else, len(ps) - 1))

        if self.particle.majorana:
            constant /= self.particle.dof
        else:
            constant /= self.particle.dof / 2

        constant *= kinematics.CollisionMultiplier3p(self)

        if not environment.get('LOGARITHMIC_TIMESTEP'):
            constant /= params.x

        stepsize *= constant_else

        try:
            ps, slice_1, slice_3 = kinematics.grid_cutoff_3p(self, ps)
        except:
            return kinematics.return_function(self, ps)

        if self.kind in [CollisionIntegralKind.Full, CollisionIntegralKind.Full_vacuum_decay] and not hasattr(self.particle, 'fast_decay'):
            C = integration_3(ps, *bounds, self.creaction, stepsize, CollisionIntegralKind.Full)
            B = integration_3(ps, *bounds, self.creaction, stepsize, CollisionIntegralKind.F_f)
            return numpy.array(slice_1 + C + slice_3) * constant, numpy.array(slice_1 + B + slice_3) * constant

        fullstack = integration_3(ps, *bounds, self.creaction, stepsize, self.kind)
        fullstack = numpy.array(slice_1 + fullstack + slice_3)

        scaled_output = kinematics.scaling(self, fullstack, constant)
        try:
            if not scaled_output:
                return kinematics.return_function(self, fullstack)
        except:
            fullstack = scaled_output

        if hasattr(self.particle, 'fast_decay'):
            if self.kind in [CollisionIntegralKind.F_decay, CollisionIntegralKind.F_f_vacuum_decay]:
                return numpy.zeros(len(fullstack)), fullstack * constant
            if self.kind in [CollisionIntegralKind.F_creation, CollisionIntegralKind.F_1_vacuum_decay]:
                return fullstack * constant, numpy.zeros(len(fullstack))

        if self.kind in [CollisionIntegralKind.F_f, CollisionIntegralKind.F_decay, CollisionIntegralKind.F_f_vacuum_decay]:
            constant *= self.particle._distribution

        return fullstack * constant
