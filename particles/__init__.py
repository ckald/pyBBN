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
from common import GRID, UNITS, kinematics, statistics as STATISTICS
from common.integrators import (
    adams_bashforth_correction, adams_moulton_solver, implicit_euler, backward_differentiation, heun_method,
    MAX_ADAMS_BASHFORTH_ORDER, MAX_ADAMS_MOULTON_ORDER, MAX_BACKWARD_DIFF_ORDER
)
from common.utils import Dynamic2DArray, DynamicRecArray
from scipy.integrate import simps
from scipy.interpolate import interp1d
from collections import Counter

from particles import DustParticle, RadiationParticle, IntermediateParticle, NonEqParticle
from interactions.four_particle.cpp.integral import distribution_interpolation, CollisionIntegralKind


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
        'thermal_dyn': True
    }

    def __init__(self, **kwargs):

        settings = dict(self._defaults)
        settings.update(kwargs)

        for key, value in settings.items():
            setattr(self, key, value)

        self.eta = 1. if self.statistics == STATISTICS.FERMION else -1.
        self.equilibrium_distribution_function = self.statistics
        self.set_grid(self.grid)
        self.decayed= False

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
        return (self.T > self.decoupling_temperature) * self.thermal_dyn

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

    def conformal_energy(self, y, conformal_mass=None):
        """ Conformal energy of the particle in comoving coordinates with evolving mass term

            \begin{equation}
                E_N = \sqrt{y^2 + (M a)^2}
            \end{equation}
        """
        if conformal_mass is None:
            conformal_mass = self.conformal_mass
        if self.mass > 0:
            return numpy.sqrt(y**2 + conformal_mass**2)
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
        self.t_decoupling = 0

        if not self.thermal_dyn and self.decoupling_temperature == 0:
            self.decoupling_temperature = params.T
        if hasattr(self, 'fast_decay'):
            self.num_creation = 0

        self.init_distribution()
        self.data['distribution'].append(self._distribution)
        self.populate_methods()
        self.update()

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
        self.oldeq = self.in_equilibrium

        # Update particle internal params only while it is in equilibrium
        if self.in_equilibrium or (not self.thermal_dyn and self.Q != 0):
            # Particle species has temperature only when it is in equilibrium
            self.init_distribution()
            self.aT = self.params.aT

        self.T = self.aT / self.params.a

        if self.in_equilibrium != self.oldeq and not self.in_equilibrium:
            # Particle decouples, have to init the distribution function array for kinetics
            self.t_decoupling = self.params.t
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

        if self.regime != oldregime:
            print("\n" + "\t"*2 + "{} changed regime at T = {:.2f} MeV from {} to {}\n"
                  .format(self.name, self.T / UNITS.MeV, oldregime.name, self.regime.name)
                  + "\t"*2 + "-"*72)

        if self.in_equilibrium != self.oldeq:
            print("\n" + "\t"*2 + "{} decoupled at T_dec = {:.2f} MeV \n"
                  .format(self.name, self.decoupling_temperature / UNITS.MeV)
                  + "\t"*2 + "-"*50)

    def update_distribution(self):
        """ Apply collision integral to modify the distribution function """
        if self.in_equilibrium:
            return

        if not numpy.all(numpy.isfinite(self.collision_integral)):
            self.collision_integral[numpy.isnan(self.collision_integral)] = 0.

        assert numpy.all(numpy.isfinite(self.collision_integral))

        if not hasattr(self, 'fast_decay'):
            self.old_distribution = self._distribution.copy()
            self._distribution += self.collision_integral * self.params.h
        else:
            if self.t_decoupling == 0:
                self.t_decoupling = self.data['params']['t'][1]
            if not self.thermal_dyn:
                self._distribution = numpy.zeros(self.grid.MOMENTUM_SAMPLES)
            else:
                decoupled_distribution = self.equilibrium_distribution(conf_mass=self.mass * self.aT / self.decoupling_temperature)
                self._distribution = decoupled_distribution * numpy.exp(-(self.params.t - self.t_decoupling) / self.lifetime)

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

        if not self.collision_integrals or self.decayed:
            return numpy.zeros(len(ps))

        if kinematics.has_decayed(self, ps):
            return numpy.zeros(len(ps))

        if hasattr(self, 'fast_decay'):
            return kinematics.Icoll_fast_decay(self, ps)

        else:
            Fs = []
            ABs = []
            Bs = []

            for integral in self.collision_integrals:
                if integral.kind in [CollisionIntegralKind.Full, CollisionIntegralKind.Full_vacuum_decay]:
                    C,B = integral.integrate(ps, stepsize=self.params.h)
                    ABs.append(C)
                    Bs.append(B)
                elif integral.kind in [CollisionIntegralKind.F_f_vacuum_decay, CollisionIntegralKind.F_decay]:
                    G = integral.integrate(ps, stepsize=self.params.h)
                    Bs.append(G)
                    ABs.append(G)
                else:
                    ABs.append(integral.integrate(ps, stepsize=self.params.h))

            AB = sum(ABs)
            B = sum(Bs)

            # Adams-Moulton method
            fs = list(self.data['collision_integral'][-MAX_ADAMS_MOULTON_ORDER:])

            I_coll = adams_moulton_solver(y=self.distribution(ps), fs=fs,
                                            A=AB, B=B, h=self.params.h)

            # # Backward differentiation method
            # ys = list(self.data['distribution'][-MAX_BACKWARD_DIFF_ORDER:])
            # I_coll = backward_differentiation(ys=ys, AB=AB, B=B, h=self.params.h)

            # # Heun method
            # I_coll =  heun_method(Is=self.data['collision_integral'], AB=AB)

            # # Implicit Euler method
            # I_coll = implicit_euler(AB=AB, B=B, h=self.params.h)

            # # Euler method
            # I_coll = AB

            return I_coll

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

        if self.in_equilibrium:
            return 1. / (
            numpy.exp(self.conformal_energy(p, self.conformal_mass) / self.aT)
            + self.eta
            )

        conformal_mass = self.mass * self.aT / self.decoupling_temperature # Approximate

        return distribution_interpolation(
            self.grid.TEMPLATE,
            self._distribution,
            p, conformal_mass,
            int(self.eta),
            self.aT,
            self.in_equilibrium
        )

    def equilibrium_distribution(self, y=None, aT=None, conf_mass=None):

        """ Equilibrium distribution that corresponds to the particle internal temperature """
        if aT is None:
            aT = self.aT
        if y is None:
            return self.equilibrium_distribution_function(
                self.conformal_energy(self.grid.TEMPLATE, conf_mass) / aT
            )
        else:
            return self.equilibrium_distribution_function(self.conformal_energy(y, conf_mass) / aT)

    def init_distribution(self, conf_mass=None):
        if not self.thermal_dyn:
            self._distribution = numpy.zeros(self.grid.MOMENTUM_SAMPLES)
        else:
            self._distribution = self.equilibrium_distribution(conf_mass=conf_mass)
            # if self.mass == 0. and self.grid.MAX_MOMENTUM > environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV:
            #     self._distribution[self.grid.TEMPLATE > environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV] = 0.
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
