# -*- coding: utf-8 -*-

import sys
import numpy
import array
from datetime import datetime
from collections import defaultdict

from common import UNITS, Params, Grid, CONST
from common import integrators, parallelization, utils
# from common.utils import PicklableObject


class Universe(object):

    """ ## Universe
        The master object that governs the calculation. """

    # System state is rendered to the log file each `log_freq` steps
    log_freq = 1

    # Controls parallelization of the collision integrals calculations
    PARALLELIZE = True

    particles = []
    interactions = []

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

        baryons_interaction = [interaction
                               for interaction in self.interactions
                               if interaction.name == "Baryons interaction"]

        if baryons_interaction:
            self.baryons_interaction = baryons_interaction[0]

        print "dx = {} MeV".format(self.params.dx / UNITS.MeV)

        self.step = 0

        self.params.update(self.total_energy_density())
        self.data = defaultdict(list, {
            'aT': array.array('f', [self.params.aT]),
            'T': array.array('f', [self.params.T]),
            'a': array.array('f', [self.params.a]),
            'x': array.array('f', [self.params.x]),
            't': array.array('f', [self.params.t]),
            'rho': array.array('f', [self.params.rho]),
            'fraction': array.array('f', [0]),
        })

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

        parallelization.pool.close()

        return self.data

    def make_step(self):
        self.params.dx = self.params.x * (numpy.exp(self.params.dy) - 1.)
        self.integrand(self.params.x, self.params.aT)

        self.log()

        order = min(self.step + 1, 5)
        fs = self.data['fraction'][-order+1:]
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
        for particle in self.particles:
            if not particle.collision_integrals:
                continue

            if self.PARALLELIZE:
                particle.collision_integral = \
                    numpy.array(parallelization.poolmap(particle, 'calculate_collision_integral',
                                                        self.grid.TEMPLATE))
            else:
                particle.collision_integral = particle.integrate_collisions()

            print particle.symbol, "I =", particle.collision_integral

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
        self.data['aT'].append(self.params.aT)
        self.data['T'].append(self.params.T)
        self.data['a'].append(self.params.a)
        self.data['x'].append(self.params.x)
        self.data['rho'].append(self.params.rho)
        self.data['t'].append(self.params.t)
        self.data['fraction'].append(self.fraction)

        #     t[s]         x    Tg[10^9K]   dTg/dt[10^9K/s] rho_tot[g cm^-3]     H[s^-1]
        # n nue->p e  p e->n nue  n->p e nue  p e nue->n  n e->p nue  p nue->n e

        if self.baryons_interaction:

            rates_map = (
                "n + ν_e ⟶  e + p",
                "n ⟶  e + ν_e' + p",
                "n + e' ⟶  ν_e' + p",
            )

            rates = []
            for reaction in rates_map:
                for integral in self.baryons_interaction.integrals:
                    if reaction in str(integral):
                        rs = integral.rates()
                        rs = [rs[0] / self.params.x * UNITS.MeV
                              / CONST.MeV_to_s_1 * CONST.rate_normalization,
                              rs[1] / self.params.x * UNITS.MeV
                              / CONST.MeV_to_s_1 * CONST.rate_normalization]
                        rates += rs
                        print integral, rs

            self.data['kawano'].append(tuple([
                self.params.t / UNITS.s,
                self.params.x / UNITS.MeV,
                self.params.T / UNITS.MeV * CONST.MeV_to_10_9K,
                self.fraction / self.params.x / UNITS.MeV * CONST.MeV_to_10_9K,
                self.params.rho / UNITS.MeV**4 * CONST.MeV4_to_g_cm_3,
                self.params.H * UNITS.s
            ] + rates))
            print "KAWANO", self.data['kawano'][-1]

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
