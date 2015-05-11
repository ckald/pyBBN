# -*- coding: utf-8 -*-

import sys
import numpy
import pandas
from datetime import datetime

from common import UNITS, Params, Grid, CONST, integrators, parallelization, utils

import kawano


class Universe(object):

    """ ## Universe
        The master object that governs the calculation. """

    # System state is rendered to the log file each `log_freq` steps
    log_freq = 1

    # Controls parallelization of the collision integrals calculations
    PARALLELIZE = True

    particles = []
    interactions = []

    kawano = None
    kawano_log = None

    data = pandas.DataFrame(columns=('aT', 'T', 'a', 'x', 't', 'rho', 'fraction'))

    def __init__(self, logfile='logs/' + str(datetime.now()) + '.txt',
                 plotting=True, params=None, grid=None):
        """
        :param particles: Set of `particle.Particle` to model
        :param interactions: Set of `interaction.Interaction` - quantum interactions \
                             between particle species
        :param logfile: Log file path (current `datetime` by default)
        """
        self.params = Params() if not params else params
        self.grid = Grid() if not grid else grid

        self.graphics = None
        if plotting:
            from plotting import Plotting
            self.graphics = Plotting()

        self.logfile = logfile
        self.init_log()

        self.fraction = 0

    def init_kawano(self, datafile='s4.dat', **kwargs):
        kawano.init_kawano(**kwargs)
        self.kawano_log = open(datafile, 'w')
        self.kawano_log.write("\t".join(kawano.heading))
        self.kawano = kawano
        self.kawano_data = pandas.DataFrame(columns=self.kawano.heading)

    def evolve(self):
        """
        ## Main computing routine

        Modeling is carried in steps of scale factor, starting from the initial conditions defined\
        by a single parameter: the initial temperature.

        Initial temperature (e.g. 10 MeV) corresponds to a point in the history of the Universe \
        long before then BBN. Then most particle species are in the thermodynamical equilibrium.

        """

        for particle in self.particles:
            print particle

        for interaction in self.interactions:
            print interaction

        self.step = 0

        self.params.update(self.total_energy_density())
        self.save()

        print '#step\tTime, s\taT, MeV\tT, MeV\tscale factor\tdx, MeV'
        self.log()
        self.step += 1

        while self.params.T > self.params.T_final:
            try:
                self.make_step()
            except KeyboardInterrupt:
                print "Keyboard interrupt!"
                break

        self.log()
        for particle in self.particles:
            print particle

        print "Data saved to file {}".format(self.logfile)

        if self.PARALLELIZE:
            parallelization.pool.close()

        return self.data

    def make_step(self):
        self.params.dx = self.params.x * (numpy.exp(self.params.dy) - 1.)
        self.integrand(self.params.x, self.params.aT)

        self.log()

        order = min(self.step + 1, 5)
        fs = self.data['fraction'].tail(order-1).values.tolist()
        fs.append(self.fraction)

        self.params.aT +=\
            integrators.adams_bashforth_correction(fs=fs, h=self.params.dy, order=order)
        self.params.x += self.params.dx

        self.params.update(self.total_energy_density())

        self.save()

        self.step += 1

    def add_particles(self, particles):
        for particle in particles:
            particle.set_params(self.params)

        self.particles += particles

    def update_particles(self):
        """ ### 1. Update particles state
            Update particle species distribution functions, check for regime switching,\
            update precalculated variables like energy density and pressure. """
        for particle in self.particles:
            particle.update()

    def init_interactions(self):
        """ ### 2. Initialize non-equilibrium interactions
            Non-equilibrium interactions of different particle species are treated by a\
            numerical integration of the Boltzmann equation for distribution functions in\
            the expanding space-time.

            Depending on the regime of the particle species involved and cosmological parameters, \
            each `Interaction` object populates `Particle.collision_integrals` array with \
            currently active `Integral` objects.
        """
        for interaction in self.interactions:
            interaction.initialize()

    def calculate_collisions(self):
        """ ### 3. Calculate collision integrals """

        particles = [particle for particle in self.particles if particle.collision_integrals]

        with utils.printoptions(precision=2):
            if self.PARALLELIZE:
                # Send the tasks
                for particle in particles:
                    particle.collision_integral = parallelization.poolmap(
                        particle, 'calculate_collision_integral',
                        self.grid.TEMPLATE
                    )
                # Collect results
                for particle in particles:
                    with utils.benchmark(lambda: "I(" + particle.symbol + ") = "
                                         + repr(particle.collision_integral)):
                        particle.collision_integral = numpy.array(particle.collision_integral
                                                                  .get(1000))
            else:
                for particle in particles:
                    with utils.benchmark(lambda: "I(" + particle.symbol + ") = "
                                         + repr(particle.collision_integral)):
                        particle.collision_integral = particle.integrate_collisions()

    def update_distributions(self):
        """ ### 4. Update particles distributions """

        for particle in self.particles:
            particle.update_distribution()

    def calculate_temperature_terms(self):
        """ ### 5. Calculate temperature equation terms """

        numerator = 0
        denominator = 0

        for particle in self.particles:
            numerator += particle.numerator()
            denominator += particle.denominator()

        return numerator, denominator

    def integrand(self, t, y):
        """ ## Temperature equation integrand

            Master equation for the temperature looks like

            \begin{equation}
                \frac{d (aT)}{dx} = \frac{\sum_i{N_i}}{\sum_i{D_i}}
            \end{equation}

            Where $N_i$ and $D_i$ represent contributions from different particle species.

            See definitions for different regimes:
              * [[Radiation|particles/RadiationParticle.py#master-equation-terms]]
              * [[Intermediate|particles/IntermediateParticle.py#master-equation-terms]]
              * [[Dust|particles/DustParticle.py#master-equation-terms]]
              * [[Non-equilibrium|particles/NonEqParticle.py#master-equation-terms]]
        """

        # 1\. Update particles states
        self.update_particles()
        # 2\. Initialize non-equilibrium interactions
        self.init_interactions()
        # 3\. Calculate collision integrals
        self.calculate_collisions()
        # 4\. Update particles distributions
        self.update_distributions()
        # 5\. Calculate temperature equation terms
        numerator, denominator = self.calculate_temperature_terms()
        self.fraction = self.params.x * numerator / denominator

        return self.fraction

    def save(self):
        """ Save current Universe parameters into the data arrays or output files """
        self.data = self.data.append({
            'aT': self.params.aT,
            'T': self.params.T,
            'a': self.params.a,
            'x': self.params.x,
            'rho': self.params.rho,
            't': self.params.t,
            'fraction': self.fraction
        }, ignore_index=True)

        if self.kawano:

            #     t[s]         x    Tg[10^9K]   dTg/dt[10^9K/s] rho_tot[g cm^-3]     H[s^-1]
            # n nue->p e  p e->n nue  n->p e nue  p e nue->n  n e->p nue  p nue->n e

            rates = self.kawano.baryonic_rates(self.params.a)

            row = {
                self.kawano.heading[0]: self.params.t / UNITS.s,
                self.kawano.heading[1]: self.params.x / UNITS.MeV,
                self.kawano.heading[2]: self.params.T / UNITS.MeV * CONST.MeV_to_10_9K,
                self.kawano.heading[3]: (self.params.T - self.data['T'].iloc[-2]) / UNITS.MeV
                / (self.params.t - self.data['t'].iloc[-2]) * UNITS.s * CONST.MeV_to_10_9K,
                self.kawano.heading[4]: self.params.rho / UNITS.MeV**4 * CONST.MeV4_to_g_cm_3,
                self.kawano.heading[5]: self.params.H * UNITS.s
            }
            row.update({self.kawano.heading[i]: rate / UNITS.MeV**5
                        for i, rate in enumerate(rates, 6)})

            self.kawano_data = self.kawano_data.append(row, ignore_index=True)
            log_entry = "\t".join("{:e}".format(item) for item in self.kawano_data.iloc[-1])

            print "KAWANO", log_entry
            self.kawano_log.write(log_entry + "\n")

    def init_log(self):
        sys.stdout = utils.Logger(self.logfile)

    def log(self):
        """ Runtime log output """

        # Print parameters every now and then
        if self.step % self.log_freq == 0:
            print '#' + str(self.step), \
                '\tt =', self.params.t / UNITS.s, \
                '\taT =', self.params.aT / UNITS.MeV, \
                '\tT =', self.params.T / UNITS.MeV, \
                '\ta =', self.params.a, \
                '\tdx =', self.params.dx / UNITS.MeV

            if self.graphics:
                self.graphics.plot(self.data)

    def total_energy_density(self):
        return sum(particle.energy_density() for particle in self.particles)
