"""
= Particles =

This file contains `Particle` class definition and code governing the switching of the dynamical\
regimes
"""

import numpy
from scipy import interpolate, integrate

from common import index_to_momentum, GRID, PARAMS, UNITS

from particles import DustParticle, RadiationParticle, IntermediateParticle, NonEqParticle


class STATISTICS:
    """ === Particle species statistics === """

    @staticmethod
    def Fermi(e):
        """ Fermi-Dirac:
            \begin{equation}
                \frac{1}{e^E + 1}
            \end{equation}
        """
        return 1. / (numpy.exp(e) + 1.)

    @staticmethod
    def Bose(e):
        """ Bose-Einstein:
            \begin{equation}
                \frac{1}{e^E - 1}
            \end{equation}
        """
        return 1. / (numpy.exp(e) - 1.)

    BOSON = Bose
    FERMION = Fermi

# Vectorized distribution functions. Just for convenience - no performance gain.
STATISTICS.FermiV = numpy.vectorize(STATISTICS.Fermi, otypes=[numpy.float_])
STATISTICS.BoseV = numpy.vectorize(STATISTICS.Bose, otypes=[numpy.float_])


class REGIMES(dict):
    """ === Particle dynamical regimes ===
        Radiation (ultra-relativistic) `RADIATION`
        :   particle mass is neglected, all values obtained analytically
        Dust (non-relativistic) `DUST`
        :   Boltzman law is used as species distribution function, simplifying the computations
        Tntermediate regime `INTERMEDIATE`
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

    def __init__(self, *args, **kwargs):

        # Functions are vectorized for convenience
        self.energy_vectorized = numpy.vectorize(self.energy, otypes=[numpy.float_])
        self.energy_normalized_vectorized = numpy.vectorize(self.energy_normalized,
                                                            otypes=[numpy.float_])
        self.integrate_collisions_vectorized = numpy.vectorize(self.integrate_collisions,
                                                               otypes=[numpy.float_])

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
            self.distribution_function_vectorized = STATISTICS.FermiV
        else:
            self.eta = -1.
            self.distribution_function = STATISTICS.Bose
            self.distribution_function_vectorized = STATISTICS.BoseV

        """ For equilibrium particles distribution function is by definition given by its\
            statistics and will not be used until species comes into non-equilibrium regime """
        self._distribution = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float_)
        self._distribution_interpolation = self.interpolate_distribution(self._distribution)
        """ Particle collision integral as well is not effective in equilibrium """
        self.collision_integral = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float_)

        """
        == Collision integral constants ==
        Collision integral for a particle starting to go out of the equilibrium is a difference\
        of 2 big numbers that cancel each other almost precisely.

        Collision integrand functional

        \begin{equation}
            \mathcal{F}(\left\\{ f_{\alpha} \right\\} ) = \
                \prod_j \prod_k f_{j_{in}} (1 \pm f_{k_{out}}) \
                - \prod_j \prod_k (1 \pm f_{j_{in}}) f_{k_{out}}
        \end{equation}

        can be basically reorganized in 2 parts: one that contains $f_{\alpha}$ and one that\
        doesn't:

        \begin{equation}
            \mathcal{F}(\left\\{ f_{\alpha} \right\\} ) = F_1 + F_f f_{\alpha}
        \end{equation}

        It is convenient because typical values of distribution functions can be very small at\
        the high momentum sector and operation of addition $(1 \pm \epsilon)$ can lead to\
        losing the $\epsilon$ and, hence, to numerical error.
        """
        self.F_f = []
        self.F_1 = []

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

        if not self.in_equilibrium:
            # For non-equilibrium particle calculate the collision integrals
            self.collision_integral = self.integrate_collisions_vectorized(GRID.TEMPLATE)
            print(self.collision_integral * UNITS.MeV)
            self._distribution += self.collision_integral * PARAMS.dx
            self._distribution_interpolation = self.interpolate_distribution(self._distribution)
            # Clear [collision integral constants](#collision-integral-constants) \
            # until the next computation step
            self.F_f = []
            self.F_1 = []

        if force_print or self.regime != oldregime or self.in_equilibrium != oldeq:
            print self

    def integrate_collisions(self, p0):
        """ == Particle collisions integration == """

        integral = None
        if self.F_1 and self.F_f:

            tmp = 1./64. / numpy.pi**3 / p0 / self.energy_normalized(p0) \
                * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H

            integrand = lambda (p1, p2): (
                self.distribution(p0) * numpy.sum([F(p0, p1, p2) for F in self.F_f])
                + numpy.sum([F(p0, p1, p2) for F in self.F_1])
            ),
            GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
            lambda x: GRID.MIN_MOMENTUM, lambda x: GRID.MAX_MOMENTUM
        )
        tmp *= integral[0] * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H
        # print p0 / UNITS.MeV, '\t', \
        #     tmp * PARAMS.dx, '\t',\
        #     integral[1] / integral[0] if integral[0] else 0

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
        particle species can be treated as `REGIMES.DUST` with a Boltzman distribution function.
        """
        regime_factor = 1e4

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
        """
        Number density of particle species. See particles/*.py for implementations
        """
        return self.regime.density(self)

    def energy_density(self):
        """
        Energy density of particle species. See particles/*.py for implementations
        """
        return self.regime.energy_density(self)

    def pressure(self):
        """
        Pressure of particle species. See particles/*.py for implementations
        """
        return self.regime.pressure(self)

    def numerator(self):
        """
        Numerator term of particle species in equation for the temperature.
        See particles/*.py for implementations
        """
        return self.regime.numerator(self)

    def denominator(self):
        """
        Denominator term of particle species in equation for the temperature.
        See particles/*.py for implementations
        """
        return self.regime.denominator(self)

    def distribution(self, p, by_index=False):
        """
        Returns interpolated value of distribution function by momentum or grid element index.

        Checks boundaries, neglects `p > GRID.MAX_MOMENTUM` and sets `p < GRID.MIN_MOMENTUM` to\
        `GRID.MIN_MOMENTUM`

        :param by_index: Return the real value of sample function by given grid point index `p`
        """
        p = abs(p)
        if self.in_equilibrium:
            if by_index:
                p = index_to_momentum(p)
            return self.distribution_function(self.energy_normalized(p) / PARAMS.aT)

        if not by_index:
            if p < GRID.MIN_MOMENTUM:
                p = GRID.MIN_MOMENTUM
            if p > GRID.MAX_MOMENTUM:
                p = GRID.MAX_MOMENTUM
            return self._distribution_interpolation(p)

        if by_index:
            if p < GRID.MOMENTUM_SAMPLES:
                if p < 0:
                    p = 0
                return self._distribution[p]
            else:
                return 0.

    def interpolate_distribution(self, distribution):
        return interpolate.interp1d(GRID.TEMPLATE, distribution,
                                    kind='quadratic', assume_sorted=True, copy=False)

    def init_distribution(self):
        self._distribution = self.distribution_function(
            self.energy_normalized_vectorized(GRID.TEMPLATE) / PARAMS.aT
        )
        self._distribution_interpolation = self.interpolate_distribution(self._distribution)

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
