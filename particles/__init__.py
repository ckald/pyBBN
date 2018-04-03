# -*- coding: utf-8 -*-
"""
# Particles

This file contains `Particle` class definition and code governing the switching of the dynamical\
regimes
"""
from __future__ import division

import numpy

import os
import environment
from common import GRID, UNITS, statistics as STATISTICS
from common.integrators import (
    adams_bashforth_correction, adams_moulton_solver, implicit_euler,
    MAX_ADAMS_BASHFORTH_ORDER, MAX_ADAMS_MOULTON_ORDER
)
from common.utils import Dynamic2DArray, DynamicRecArray

from particles import DustParticle, RadiationParticle, IntermediateParticle, NonEqParticle
# from interactions.four_particle.integral import distribution_interpolation
from interactions.four_particle.cpp.integral import distribution_interpolation


class REGIMES(dict):
    """ ## Particle dynamical regimes
        Radiation (ultra-relativistic) `RADIATION`
        :   particle mass is neglected, all values obtained analytically
        Dust (non-relativistic) `DUST`
        :   Boltzman law is used as species distribution function, simplifying the computations
        Intermediate regime `INTERMEDIATE`
        :   values are computed explicitly using the precise form of the Bose-Einstein and\
            Fermi-Dirac distributions including particle mass term
        Non-equilibrium `NONEQ`
        :   particle quantum interactions have to be computed explicitly by solving Boltzman \
            equation; individual distribution function of the species becomes utilized to obtain\
            any species-related values
        """
    RADIATION = RadiationParticle
    DUST = DustParticle
    INTERMEDIATE = IntermediateParticle
    NONEQ = NonEqParticle


class AbstractParticle:

    _defaults = {
        'mass': 0 * UNITS.eV,
        'decoupling_temperature': 0 * UNITS.eV,
        'name': 'Particle',
        'symbol': 'p',
        'dof': 2,
        'statistics': STATISTICS.FERMION,
        'params': None,
        'grid': GRID,
        'thermal': True
    }

    def __init__(self, **kwargs):

        settings = dict(self._defaults)
        settings.update(kwargs)

        for key, value in settings.items():
            setattr(self, key, value)

        self.eta = 1. if self.statistics == STATISTICS.FERMION else -1.
        self.equilibrium_distribution_function = self.statistics
        self.set_grid(self.grid)

        self.collision_integrals = []
        self.data = {
            'distribution': Dynamic2DArray(self.grid.TEMPLATE),
            'collision_integral': Dynamic2DArray(self.grid.TEMPLATE),
            'params': DynamicRecArray([
                ['a', '', 1],
                ['t', 's', UNITS.s],
                ['T', 'MeV', UNITS.MeV],
                ['density', 'MeV^3', UNITS.MeV**3],
                ['energy_density', 'MeV^4', UNITS.MeV**4]
            ])
        }
        # self.data['distribution'].append(self._distribution)

        if not self.thermal:
            self.decoupling_temperature = 1e20 * UNITS.MeV

        if self.params:
            self.set_params(self.params)

    def __str__(self):
        """ String-like representation of particle species it's regime and parameters """
        return "{} ({}, {})\nn = {:e} MeV^3, rho = {:e} MeV^4\n".format(
            self.name,
            "eq" if self.in_equilibrium else "non-eq",
            self.regime.name,
            self.density / UNITS.MeV**3,
            self.energy_density / UNITS.MeV**4
        ) + ("-" * 80)

    def __repr__(self):
        return self.symbol

    def __gt__(self, other):
        return self.name > other.name

    def populate_methods(self):
        regime = self.regime
        self.density = regime.density(self)
        self.energy_density = regime.energy_density(self)
        self.pressure = regime.pressure(self)
        self.entropy = regime.entropy(self)
        self.numerator = lambda: regime.numerator(self)
        self.denominator = lambda: regime.denominator(self)

    @property
    def regime(self):
        """
        ### Regime-switching ratio
        For ultra-relativistic particles the mass is effectively `0`. This implies that all\
        computed numerically values can be as well obtained analytically: energy density, pressure,\
        etc.

        Let particle mass be equal $M$ and regime factor equal $\gamma$. As soon as the \
        temperature of the system $T$ drops to the value about $ M \gamma $, particle should be \
        switched to the computation regime where its mass is also considered: \
        `REGIMES.INTERMEDIATE`. When $T$ drops down even further to the value $ M / \gamma $,\
        particle species can be treated as `REGIMES.DUST` with a Boltzmann distribution function.
        """
        regime_factor = environment.get('REGIME_SWITCHING_FACTOR')

        if not self.in_equilibrium:
            return REGIMES.NONEQ

        if self.T > self.mass * regime_factor:
            regime = REGIMES.RADIATION
        elif self.T * regime_factor < self.mass:
            regime = REGIMES.DUST
        else:
            regime = REGIMES.INTERMEDIATE

        return regime

    @property
    def in_equilibrium(self):
        """ Simple check for equilibrium """
        return self.T > self.decoupling_temperature

    def energy(self, p):
        """ Physical energy of the particle

            \begin{equation}
                E = \sqrt{p^2 + M^2}
            \end{equation}
        """

        if self.mass > 0:
            return numpy.sqrt(p**2 + self.mass**2)
        else:
            return numpy.absolute(p)

    def conformal_energy(self, y):
        """ Conformal energy of the particle in comoving coordinates with evolving mass term

            \begin{equation}
                E_N = \sqrt{y^2 + (M a)^2}
            \end{equation}
        """
        if self.mass > 0:
            return numpy.sqrt(y**2 + self.conformal_mass**2)
        else:
            return numpy.absolute(y)

    @property
    def conformal_mass(self):
        """ In the expanding Universe, distribution function of the massive particle in the\
            conformal coordinates changes because of the evolving mass term $M a$ """
        return self.mass * self.params.a


class DistributionParticle(AbstractParticle):

    """ ## Particle
        Master-class for particle species. Realized as finite state machine that switches to\
        different regime when temperature becomes comparable to the particle mass or drops below\
        particle `decoupling_temperature`
    """

    def set_params(self, params):
        """ Set internal parameters using arguments or default values """
        self.params = params
        self.T = params.T
        self.aT = params.aT

        self.update()
        self.init_distribution()

        self.populate_methods()

    def set_grid(self, grid):
        self.grid = grid
        """ For equilibrium particles distribution function is by definition given by its\
            statistics and will not be used until species comes into non-equilibrium regime """
        self._distribution = numpy.zeros(self.grid.MOMENTUM_SAMPLES, dtype=numpy.float_)
        """ Particle collision integral is not effective in the equilibrium as well """
        self.collision_integral = numpy.zeros(self.grid.MOMENTUM_SAMPLES, dtype=numpy.float_)

    def update(self, force_print=False):
        """ Update the particle parameters according to the new state of the system """
        oldregime = self.regime
        oldeq = self.in_equilibrium

        # Update particle internal params only while it is in equilibrium
        if self.in_equilibrium:
            # Particle species has temperature only when it is in equilibrium
            self.aT = self.params.aT

        self.T = self.aT / self.params.a

        if self.in_equilibrium != oldeq and not self.in_equilibrium:
            # Particle decouples, have to init the distribution function array for kinetics
            self.init_distribution()

        self.collision_integral = numpy.zeros(self.grid.MOMENTUM_SAMPLES, dtype=numpy.float_)

        self.populate_methods()

        self.data['params'].append({
            'a': self.params.a,
            't': self.params.t,
            'T': self.params.T,
            'density': self.density,
            'energy_density': self.energy_density
        })

        if force_print or self.regime != oldregime or self.in_equilibrium != oldeq:
            print("\n" + "\t"*2 + "{} decoupled at T_dec = {:.2f} MeV \n"
                  .format(self.name, self.decoupling_temperature / UNITS.MeV)
                  + "\t"*2 + "-"*50)

    def update_distribution(self):
        """ Apply collision integral to modify the distribution function """
        if self.in_equilibrium:
            return

        assert numpy.all(numpy.isfinite(self.collision_integral))

        self.old_distribution = self._distribution.copy()
        self._distribution += self.collision_integral * self.params.h

        # assert all(self._distribution >= 0), self._distribution
        self._distribution = numpy.maximum(self._distribution, 0)

        # Clear collision integrands for the next computation step
        self.collision_integrals = []
        self.data['collision_integral'].append(self.collision_integral)
        self.data['distribution'].append(self._distribution)

    def integrate_collisions(self):
        return self.calculate_collision_integral(self.grid.TEMPLATE)

    def calculate_collision_integral(self, ps):
        """ ### Particle collisions integration """

        if not self.collision_integrals:
            return numpy.zeros(len(ps))

        # # Adams-Bashforth integrator is unstable for unknown reason at the moment
        # if not environment.get('SPLIT_COLLISION_INTEGRAL'):
        #     fs = [sum([integral.integrate(ps, stepsize=self.params.h)
        #                for integral in self.collision_integrals])]

        #     fs = list(self.data['collision_integral'][-MAX_ADAMS_BASHFORTH_ORDER:]) + fs

        #     return adams_bashforth_correction(fs=fs, h=self.params.h) / self.params.h

        if not environment.get('SPLIT_COLLISION_INTEGRAL'):
            fs = sum([integral.integrate(ps, stepsize=self.params.h)
                       for integral in self.collision_integrals])
            return fs

        As = []
        Bs = []

        for integral in self.collision_integrals:
            A, B = integral.integrate(ps, stepsize=self.params.h)
            As.append(A)
            Bs.append(B)

        A = sum(As)
        B = sum(Bs)

        if environment.get('ADAMS_MOULTON_DISTRIBUTION_CORRECTION'):
            fs = list(self.data['collision_integral'][-MAX_ADAMS_MOULTON_ORDER:])

            prediction = adams_moulton_solver(y=self.distribution(ps), fs=fs,
                                              A=A, B=B, h=self.params.h)
        else:
            prediction = implicit_euler(y=self.distribution(ps), t=None,
                                        A=A, B=B, h=self.params.h)

        return (prediction - self.distribution(ps)) / self.params.h

    def distribution(self, p):
        """
        ## Distribution function interpolation

        Two possible strategies for the distribution function interpolation are implemented:

          * linear interpolation
          * exponential interpolation

        While linear interpolation is exceptionally simple, the exponential interpolation gives\
        exact results for the equilibrium functions - thus collision integral for them almost\
        exactly cancels out unlike the case of linear interpolation.

        Distribution functions are continuously evaluated each during the simulation, so to avoid\
        excessive evaluation of the costly logarithms and exponentials, this function first checks\
        if the current momenta value coincides with any of the grid points.
        """

        return distribution_interpolation(
            self.grid.TEMPLATE,
            self._distribution,
            p, self.conformal_mass,
            int(self.eta),
            self.aT,
            self.in_equilibrium
        )

    def equilibrium_distribution(self, y=None, aT=None):

        """ Equilibrium distribution that corresponds to the particle internal temperature """
        if aT is None:
            aT = self.aT
        if y is None:
            return self.equilibrium_distribution_function(
                self.conformal_energy(self.grid.TEMPLATE) / aT
            )
        else:
            return self.equilibrium_distribution_function(self.conformal_energy(y) / aT)

    def init_distribution(self):
        if not self.thermal:
            self._distribution = numpy.zeros(self.grid.MOMENTUM_SAMPLES)
        else:
            self._distribution = self.equilibrium_distribution()
        return self._distribution


class AdaptiveDistributionParticle(DistributionParticle):
    energy_limit = numpy.inf

    def integrate_collisions(self):
        collision_integral = self.grid.TEMPLATE.copy()
        for ix, p0 in enumerate(self.grid.TEMPLATE):
            if p0 < self.energy_limit * self.params.a:
                collision_integral[ix] = self.calculate_collision_integral(p0)
            else:
                collision_integral[ix] = 0
        return collision_integral


Particle = DistributionParticle
