import numpy
import numericalunits as nu
from scipy import interpolate, integrate

from common import Distributions, momentum_to_index, index_to_momentum,\
    GRID, STATISTICS, lambda_integrate, PARAMS, REGIMES, benchmark

from regimes import DustParticle, RadiationParticle, IntermediateParticle, NonEqParticle


"""
Dict for easy access to specific particle regimes functions
"""
ROUTINES = {
    REGIMES.DUST: DustParticle,
    REGIMES.RADIATION: RadiationParticle,
    REGIMES.INTERMEDIATE: IntermediateParticle,
    REGIMES.NONEQ: NonEqParticle
}


class Particle():
    """
    Master-class for particle species. Realized as finite state machine that switches to different
    regime when temperature becomes comparable to particle mass or drop below particle
    decoupling_temperature
    """

    def __init__(self, *args, **kwargs):
        # energy function is vectorized for faster work on numpy.arrays
        self.energy_vectorized = numpy.vectorize(self.energy, otypes=[numpy.float_])
        self.energy_normalized_vectorized = numpy.vectorize(self.energy_normalized,
                                                            otypes=[numpy.float_])

        self.integrate_collisions_vectorized = numpy.vectorize(self.integrate_collisions,
                                                               otypes=[numpy.float_])

        self.temperature = PARAMS.T
        self.mass = kwargs.get('mass', 0 * nu.eV)
        self.decoupling_temperature = kwargs.get('decoupling_temperature', 0 * nu.eV)
        self.name = kwargs.get('name', 'Particle')

        self.dof = kwargs.get('dof', 2)  # particle species degeneracy

        self.statistics = kwargs.get('statistics', STATISTICS.FERMION)
        if self.statistics == STATISTICS.FERMION:
            self.eta = 1.
            self.distribution_function = Distributions.Fermi
            self.distribution_function_vectorized = Distributions.FermiV
        else:
            self.eta = -1.
            self.distribution_function = Distributions.Bose
            self.distribution_function_vectorized = Distributions.BoseV

        self._distribution = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float_)
        self._distribution_interpolation = interpolate.interp1d(GRID.TEMPLATE,
                                                                self._distribution,
                                                                kind='cubic')
        self.collision_integral = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float_)

        self.update(force_print=True)
        self.init_distribution()

        self.F_f = []
        self.F_1 = []

    def __str__(self):
        """ String-like representation of particle species it's regime and parameters """
        return "%s (%s, %s)\nn = %s, rho = %s\n" % (
            self.name,
            "eq" if self.in_equilibrium else "non-eq",
            self.regime,
            self.density() / nu.eV**3,
            self.energy_density() / nu.eV**4
        ) + ("-" * 80)

    def __repr__(self):
        return self.name

    def update(self, force_print=False):
        """ Evolution """
        oldregime = self.regime
        oldeq = self.in_equilibrium

        self.temperature = PARAMS.T

        # clear saved values of density, energy_density and pressure
        self._density = None
        self._energy_density = None
        self._pressure = None

        if self.in_equilibrium != oldeq and not self.in_equilibrium:
            # particle decouples, init distribution function array for kinetics
            self.init_distribution()

        if not self.in_equilibrium:
            self.collision_integral = self.integrate_collisions_vectorized(GRID.TEMPLATE)
            self._distribution += PARAMS.dx * self.collision_integral
            self._distribution_interpolation = interpolate.interp1d(GRID.TEMPLATE,
                                                                    self._distribution,
                                                                    kind='cubic')

            self.F_f = []
            self.F_1 = []

        if force_print or self.regime != oldregime or self.in_equilibrium != oldeq:
            print self

    def integrate_collisions(self, p0):
        tmp = 1./64. / numpy.pi**3 / p0 / self.energy_normalized(p0)

        integral = integrate.dblquad(
            lambda p1, p2: (
                self.distribution(p0) * numpy.sum([F(p0, p1, p2) for F in self.F_f])
                + numpy.sum([F(p0, p1, p2) for F in self.F_1])
            ),
            GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
            lambda x: GRID.MIN_MOMENTUM, lambda x: GRID.MAX_MOMENTUM
        )
        tmp *= integral[0] * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H
        print p0 / nu.MeV, '\t', \
            tmp * PARAMS.dx, '\t',\
            integral[1] / integral[0] if integral[0] else 0

        return tmp

    @property
    def regime(self):
        """
        Uses particle species temperature to check correspondence of particle species to one of the
        dynamic regimes: radiation, dust, intermediate or non-equilibrium
        """
        if not self.in_equilibrium:
            return REGIMES.NONEQ

        if self.temperature > self.mass * PARAMS.regime_factor:
            regime = REGIMES.RADIATION
        elif self.temperature * PARAMS.regime_factor < self.mass:
            regime = REGIMES.DUST
        else:
            regime = REGIMES.INTERMEDIATE

        return regime

    def density(self):
        """
        Number density of particle species. See regimes/*.py for implementations
        """
        return ROUTINES[self.regime].density(self)

    def energy_density(self):
        """
        Energy density of particle species. See regimes/*.py for implementations
        """
        return ROUTINES[self.regime].energy_density(self)

    def pressure(self):
        """
        Pressure of particle species.
        See regimes/*.py for implementations
        """
        return ROUTINES[self.regime].pressure(self)

    def numerator(self):
        """
        Numerator term of particle species in equation for the temperature.
        See regimes/*.py for implementations
        """
        return ROUTINES[self.regime].numerator(self)

    def denominator(self):
        """
        Denominator term of particle species in equation for the temperature.
        See regimes/*.py for implementations
        """
        return ROUTINES[self.regime].denominator(self)

    def distribution(self, p, by_index=False):
        """
        Returns value of distribution function by momentum or grid element index.
        Checks boundaries, neglects p > GRID.MAX_MOMENTUM and sets p < GRID.MIN_MOMENTUM to
        GRID.MIN_MOMENTUM
        """
        p = abs(p)
        if self.in_equilibrium:
            if by_index:
                p = index_to_momentum(p)
            return self.distribution_function(self.energy_normalized(p) / PARAMS.aT)

        if not by_index:
            # i = momentum_to_index(p)
            if p < GRID.MIN_MOMENTUM:
                p = GRID.MIN_MOMENTUM
            if p > GRID.MAX_MOMENTUM:
                p = GRID.MAX_MOMENTUM
            return self._distribution_interpolation(numpy.abs(p))

        if by_index:
            if p < GRID.MOMENTUM_SAMPLES:
                if p < 0:
                    p = 0
                return self._distribution[p]
            else:
                return 0.

    def init_distribution(self):
        if self.statistics == STATISTICS.FERMION:
            distribution_function = Distributions.FermiV
        else:
            distribution_function = Distributions.BoseV

        # with benchmark(self.name + " init " + str(self.statistics)):
        self._distribution = distribution_function(
            self.energy_normalized_vectorized(GRID.TEMPLATE) / self.temperature
        )

    @property
    def in_equilibrium(self):
        """ Simple check for equilibrium """
        return self.temperature > self.decoupling_temperature

    def energy(self, p):
        """ Physical energy of the particle """
        if self.mass > 0:
            return numpy.sqrt(p**2 + self.mass**2, dtype=numpy.float_)
        else:
            return abs(p)

    def energy_normalized(self, y):
        """ Normalized energy of the particle in comoving coordinates with evolving mass term """
        if self.mass > 0:
            return numpy.sqrt(y**2 + self.mass_normalized**2, dtype=numpy.float_)
        else:
            return abs(y)

    @property
    def mass_normalized(self):
        return self.mass * PARAMS.a
