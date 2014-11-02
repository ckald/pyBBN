# -*- coding: utf-8 -*-
"""
= Particles =

This file contains `Particle` class definition and code governing the switching of the dynamical\
regimes
"""
from __future__ import division
import functools
import math

import numpy
from skmonaco import mcquad

from scipy import integrate

from plotting import plot_integrand
from common import GRID, PARAMS, UNITS, benchmark, smoothe_function

from particles import DustParticle, RadiationParticle, IntermediateParticle, NonEqParticle
# from particles.interpolation import distribution_interpolation


class STATISTICS:
    """ === Particle species statistics === """

    @staticmethod
    @numpy.vectorize
    def Fermi(e):
        """ Fermi-Dirac:
            \begin{equation}
                \frac{1}{e^E + 1}
            \end{equation}
        """
        return 1. / (numpy.exp(e) + 1.)

    @staticmethod
    @numpy.vectorize
    def Bose(e):
        """ Bose-Einstein:
            \begin{equation}
                \frac{1}{e^E - 1}
            \end{equation}
        """
        return 1. / (numpy.exp(e) - 1.)

    BOSON = Bose
    FERMION = Fermi


class REGIMES(dict):
    """ === Particle dynamical regimes ===
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


class Particle():

    """ == Particle ==
        Master-class for particle species. Realized as finite state machine that switches to\
        different regime when temperature becomes comparable to the particle mass or drops below\
        particle `decoupling_temperature`
    """

    COLLISION_INTEGRATION_METHOD = ['mcquad', 'dblquad'][0]

    def __init__(self, *args, **kwargs):

        """ Set internal parameters using arguments or default values """
        self.T = PARAMS.T
        self.aT = PARAMS.aT
        self.mass = kwargs.get('mass', 0 * UNITS.eV)
        self.decoupling_temperature = kwargs.get('decoupling_temperature', 0 * UNITS.eV)
        self.name = kwargs.get('name', 'Particle')

        self.dof = kwargs.get('dof', 2)  # particle species degeneracy (e.g., spin-degeneracy)

        self.statistics = kwargs.get('statistics', STATISTICS.FERMION)
        if self.statistics == STATISTICS.FERMION:
            self.eta = 1.
            self.distribution_function = STATISTICS.Fermi
        else:
            self.eta = -1.
            self.distribution_function = STATISTICS.Bose

        """ For equilibrium particles distribution function is by definition given by its\
            statistics and will not be used until species comes into non-equilibrium regime """
        self._distribution = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float_)
        """ Particle collision integral is not effective in the equilibrium as well """
        self.collision_integral = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float_)

        self.collision_integrands = []

        self.update()
        self.init_distribution()

    def __str__(self):
        """ String-like representation of particle species it's regime and parameters """
        return "%s (%s, %s)\nn = %s, rho = %s\n" % (
            self.name,
            "eq" if self.in_equilibrium else "non-eq",
            self.regime.name,
            self.density() / UNITS.eV**3,
            self.energy_density() / UNITS.eV**4
        ) + ("-" * 80)

    def __repr__(self):
        return self.name

    def update(self, force_print=False):
        """ Update the particle parameters according to the new state of the system """
        oldregime = self.regime
        oldeq = self.in_equilibrium

        # Clear saved values of density, energy_density and pressure
        self._density = None
        self._energy_density = None
        self._pressure = None

        # Update particle internal params only while it is in equilibrium
        if self.in_equilibrium:
            # Particle species has temperature only when it is in equilibrium
            self.T = PARAMS.T
            self.aT = PARAMS.aT

        if self.in_equilibrium != oldeq and not self.in_equilibrium:
            # Particle decouples, have to init the distribution function array for kinetics
            self.init_distribution()

        if force_print or self.regime != oldregime or self.in_equilibrium != oldeq:
            print self

    def check_step_size(self, delta):
        """
        Scale factor step size controller.

        TODO: not usable now, need to orchestrate between many collision integrations - \
              i.e., select smallest suggested
        """
        relative_delta = numpy.absolute(delta / self._distribution).max()
        if relative_delta < 0.1:
            PARAMS.dx *= 2
            delta *= 2
            print "//// Step size increased to", PARAMS.dx / UNITS.MeV
        elif relative_delta > 0.2:
            PARAMS.dx /= 2
            delta /= 2
            print "//// Step size decreased to", PARAMS.dx / UNITS.MeV

    def update_distribution(self):
        if self.in_equilibrium:
            return

        # Smooth the collision integral a bit to eliminate the Monte Carlo errors
        # self.collision_integral = smoothe_function(self.collision_integral, GRID.MOMENTUM_SAMPLES)

        delta = self.collision_integral * PARAMS.dx

        prediction = self._distribution + delta
        if numpy.any(prediction < 0):
            print("Holy cow, negative distribution!")

        # self.check_step_size(delta)

        # Force minimal distribution function value to 0
        self._distribution = numpy.fmax(prediction, 0)

        # Clear collision integrands for the next computation step
        self.collision_integrands = []

    def integrate_collisions(self, p0):
        """ == Particle collisions integration == """

        if not self.collision_integrands:
            return 0

        tmp = 1./64. / numpy.pi**3 * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H

        integrand = lambda (p1, p2): sum([func(p0, p1, p2) for func in self.collision_integrands])

        # plot_integrand(integrand, self, p0)

        with benchmark("p0 = {:.2f}\t".format(p0 / UNITS.MeV)):
            if self.COLLISION_INTEGRATION_METHOD == 'dblquad':
                integral, error = integrate.dblquad(
                    lambda p1, p2: integrand((p1, p2)),
                    GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
                    lambda p1: GRID.MIN_MOMENTUM, lambda p1: min(p0 + p1, GRID.MAX_MOMENTUM)
                )
            elif self.COLLISION_INTEGRATION_METHOD == 'mcquad':
                integral, error = mcquad(
                    integrand,
                    xl=[GRID.MIN_MOMENTUM, GRID.MIN_MOMENTUM],
                    xu=[GRID.MAX_MOMENTUM, GRID.MAX_MOMENTUM],
                    npoints=1e4
                )
            print '{name:}\t{integral: .5e}\t±{error: .5e}\t({relerror: 8.2f}%)\t'\
                .format(name=self.name,
                        integral=integral / UNITS.MeV,
                        error=error / UNITS.MeV,
                        relerror=(error / integral * 100 if integral else 0)
                        ),

        tmp *= integral
        return tmp

    @property
    def regime(self):
        """
        === Regime-switching ratio ===
        For ultra-relativistic particles the mass is effectively `0`. This implies that all\
        computed numerically values can be as well obtained analytically: energy density, pressure,\
        etc.

        Let particle mass be equal $M$ and regime factor equal $\gamma$. As soon as the \
        temperature of the system $T$ drops to the value about $ M \gamma $, particle should be \
        switched to the computation regime where its mass is also considered: \
        `REGIMES.INTERMEDIATE`. When $T$ drops down even further to the value $ M / \gamma $,\
        particle species can be treated as `REGIMES.DUST` with a Boltzmann distribution function.
        """
        regime_factor = 1e2

        if not self.in_equilibrium:
            return REGIMES.NONEQ

        if self.T > self.mass * regime_factor:
            regime = REGIMES.RADIATION
        elif self.T * regime_factor < self.mass:
            regime = REGIMES.DUST
        else:
            regime = REGIMES.INTERMEDIATE

        return regime

    # TODO: probably just remove these and always use `particle.regime.something()`?
    def density(self):
        return self.regime.density(self)

    def energy_density(self):
        return self.regime.energy_density(self)

    def pressure(self):
        return self.regime.pressure(self)

    def numerator(self):
        return self.regime.numerator(self)

    def denominator(self):
        return self.regime.denominator(self)

    def distribution(self, p):
        """
        Returns interpolated value of distribution function by momentum.
        """
        exponential_interpolation = False  # True
        p = abs(p)
        if self.in_equilibrium or p > GRID.MAX_MOMENTUM:
            return self.distribution_function(self.energy_normalized(p) / PARAMS.aT)

        index = numpy.where(GRID.TEMPLATE == p)
        if len(index[0]):
            return self._distribution[index[0][0]]

        # return distribution_interpolation(GRID.TEMPLATE, self._distribution, p,
        #                                   # energy_normalized=self.energy_normalized,
        #                                   eta=self.eta)

        p_low = None  # GRID.MIN_MOMENTUM
        p_high = None  # GRID.MAX_MOMENTUM

        i = 0
        for point in GRID.TEMPLATE:
            if point == p:
                return self._distribution[i]
            elif point < p:
                p_low = point
                i_low = i
            elif point > p:
                p_high = point
                i_high = i
                break
            i += 1

        if p_low is None:
            raise Exception("Outside of interpolated range: {}".format(p))

        if exponential_interpolation:
            E_p = self.energy_normalized(p)
            E_low = self.energy_normalized(p_low)
            E_high = self.energy_normalized(p_high)

            """
            \begin{equation}
                g = \frac{ (E_p - E_low) g_high + (E_high - E_p) g_low }{ (E_high - E_low) }
            \end{equation}
            """

            g_high = self._distribution[i_high]
            g_low = self._distribution[i_low]

            if g_high > 0:
                g_high = (1. / g_high - 1.)
                if g_high > 0:
                    g_high = math.log(g_high)

            if g_low > 0:
                g_low = (1. / g_low - 1.)
                if g_low > 0:
                    g_low = math.log(g_low)

            g = ((E_p - E_low) * g_high + (E_high - E_p) * g_low) / (E_high - E_low)

            return 1. / (numpy.exp(g) + self.eta)

        else:
            return (
                self._distribution[i_low] * (p_high - p) + self._distribution[i_high] * (p - p_low)
            ) / (p_high - p_low)

    def init_distribution(self):
        self._distribution = self.distribution_function(
            numpy.vectorize(self.energy_normalized)(GRID.TEMPLATE) / self.aT
        )
        return self._distribution

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
            return numpy.sqrt(p**2 + self.mass**2, dtype=numpy.float_)
        else:
            return abs(p)

    def energy_normalized(self, y):
        """ Normalized energy of the particle in comoving coordinates with evolving mass term

            \begin{equation}
                E_n = \sqrt{y^2 + (M a)^2}
            \end{equation}
        """
        if self.mass > 0:
            return numpy.sqrt(y**2 + self.mass_normalized**2, dtype=numpy.float_)
        else:
            return abs(y)

    @property
    def mass_normalized(self):
        return self.mass * PARAMS.a
