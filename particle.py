import numpy
import numericalunits as nu

from common import Distributions, momentum_to_index, index_to_momentum,\
    GRID, STATISTICS, lambda_integrate, PARAMS, REGIMES, benchmark

from regimes import DustParticle, RadiationParticle, IntermediateParticle, NonEqParticle


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

        self.energy_vectorized = numpy.vectorize(self.energy)

        self.temperature = PARAMS.T
        self.mass = kwargs.get('mass', 0 * nu.eV)
        self.decoupling_temperature = kwargs.get('decoupling_temperature', 0 * nu.eV)
        self.name = kwargs.get('name', 'Particle')
        self.dof = kwargs.get('dof', 2)
        self.statistics = kwargs.get('statistics', STATISTICS.FERMION)
        if self.statistics == STATISTICS.FERMION:
            self.eta = 1.
            self.distribution_function = Distributions.Fermi
        else:
            self.eta = -1.
            self.distribution_function = Distributions.Bose

        self._distribution = numpy.zeros(GRID.MOMENTUM_SAMPLES, dtype=numpy.float64)

        self.update(force_print=True)
        self.init_distribution()

    def __str__(self):
        return "%s (%s, %s)\nn = %s, rho = %s, p = %s\n" % (
            self.name,
            "eq" if self.in_equilibrium else "non-eq",
            self.regime,
            self.density() / nu.eV**3,
            self.energy_density() / nu.eV**4,
            self.pressure() / nu.eV**4
        ) + ("-" * 80)

    def update(self, force_print=False):
        oldregime = self.regime
        oldeq = self.in_equilibrium

        self.temperature = PARAMS.T

        # particle decouples
        if self.in_equilibrium != oldeq and not self.in_equilibrium:
            self.init_distribution()

        if force_print or self.regime != oldregime or self.in_equilibrium != oldeq:
            print
            print self

    @property
    def regime(self):
        if not self.in_equilibrium:
            return REGIMES.NONEQ

        if self.temperature > self.mass * 10:
            regime = REGIMES.RADIATION
        elif self.temperature * 10 < self.mass:
            regime = REGIMES.DUST
        else:
            regime = REGIMES.INTERMEDIATE

        return regime

    def density(self):
        return ROUTINES[self.regime].density(self)

    def energy_density(self):
        return ROUTINES[self.regime].energy_density(self)

    def pressure(self):
        return ROUTINES[self.regime].pressure(self)

    def energy_density_rate(self):
        return ROUTINES[self.regime].energy_density_rate(self)

    def distribution(self, p, by_momentum=True, by_index=False):
        if self.in_equilibrium:
            if by_index:
                p = index_to_momentum(p)
            return self.distribution_function(p)

        if by_momentum:
            if p <= GRID.MAX_MOMENTUM:
                # Treat p < MIN_MOMENTUM as MIN_MOMENTUM
                if p < GRID.MIN_MOMENTUM:
                    p = GRID.MIN_MOMENTUM
                return self._distribution[momentum_to_index(p)]
            # neglect p > MAX_MOMENTUM
            else:
                return 0.

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
            self.energy_vectorized(GRID.TEMPLATE) / self.temperature
        )

    @property
    def in_equilibrium(self):
        return self.temperature > self.decoupling_temperature

    def energy(self, p):
        if self.mass > 0:
            return numpy.sqrt(p**2 + self.mass**2, dtype=numpy.float_)
        else:
            return p
