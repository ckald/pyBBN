# -*- coding: utf-8 -*-
import numpy

import environment
from common.integrators import paired as paired_integrators
from interactions.boltzmann import BoltzmannIntegral
from interactions.three_particle.cpp.integral import integrand_3, grid_t3, particle_t3, reaction_t3


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

        self.creaction = None
        self.cMs = None

        if not self.cMs:
            self.cMs = sum(M.K for M in self.Ms)

    def integrate(self, ps, bounds=None):
        if self.reaction[0].specie.mass == 0:
            return 0, 0
        
        if bounds is None:
            bounds = tuple(b1 * self.par.a for b1 in self.grids.BOUNDS)

        if not self.creaction:
            self.creaction = [
                reaction_t3(
                    specie=particle_t3(
                        m=particle.specie.conformal_mass,
                        grid=grid_t3(
                            grid=particle.specie.grid.TEMPLATE,
                            distribution=particle.specie._distribution
                        ),
                        eta=int(particle.specie.eta),
                        in_equilibrium=int(particle.specie.in_equilibrium),
                        T=particle.specie.aT
                    ),
                    side=particle.side
                )
                for particle in self.reaction
            ]

        def prepared_integrand(p1):
            integrand_1, integrand_f = integrand_3(ps, p1, self.cMs, self.creaction)
            return integrand_1, integrand_f


        if ps == 0:
            constant = self.par.m / self.par.x / 8 / numpy.pi / self.par.H / self.particle.mass**2

            integral_1, integral_f = prepared_integrand(numpy.zeros(len(self.grids.TEMPLATE)))
            integral_1 = integral_1[0]
            integral_f = integral_f[0]

        else:
            constant = self.par.x / self.par.m / 32. / numpy.pi / self.par.H

            if environment.get('SIMPSONS_INTEGRATION'):
                integral_1, integral_f = paired_integrators.integrate_1D_simpsons(
                    prepared_integrand(self.grids.TEMPLATE),
                    grid=self.grids.TEMPLATE
                )
            else:
                integral_1, integral_f = paired_integrators.integrate_1D(
                    prepared_integrand,
                    bounds=bounds
                )

        if not environment.get('LOGARITHMIC_TIMESTEP'):
            constant /= self.par.x

        return constant * integral_1, constant * integral_f
