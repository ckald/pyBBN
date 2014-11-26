# -*- coding: utf-8 -*-

import sys
import copy
import numpy
import array
from datetime import datetime

from common import UNITS, PARAMS, GRID
from common import integrators, parallelization, utils
from plotting import Plotting


class Universe:

    """ == Universe ==
        The master object that governs the calculation. """

    # System state is rendered to the log file each `log_freq` steps
    log_freq = 1
    """ Temperature equation integration method: explicit Euler or two-step Heun's method """
    INTEGRATION_METHOD = ['euler', 'heun'][0]

    # Controls parallelization of the collision integrals calculations
    PARALLELIZE = True

    # Set size increasing multiplier
    step_size_multiplier = 1.1

    def __init__(self, particles=[], interactions=[],
                 logfile='logs/' + str(datetime.now()) + '.txt',
                 plotting=True,
                 adaptive_step_size=True,
                 postmortem_debugger=True):
        """
        :param particles: Set of `particle.Particle` to model
        :param interactions: Set of `interaction.Interaction` - quantum interactions \
                             between particle species
        :param logfile: Log file path (current `datetime` by default)
        :param adaptive_step_size: Boolean, whether to control the step size or not
        :param postmortem_debugger: Boolean, whether to invoke the `pdb` debugger at the end of\
                                    computation
        """
        self.particles = particles
        self.interactions = interactions

        self.graphics = Plotting() if plotting else None
        self.adaptive_step_size = adaptive_step_size
        self.postmortem_debugger = postmortem_debugger

        sys.stdout = utils.Logger(logfile)
        self.logger = sys.stdout

    def evolve(self):
        """
        == Main computing routine ==

        Modeling is carried in steps of scale factor, starting from the initial conditions defined\
        by a single parameter: the initial temperature.

        Initial temperature (e.g. 10 MeV) corresponds to a point in the history of the Universe \
        long before then BBN. Then most particle species are in the thermodynamical equilibrium.

        """

        print "dx = {} MeV".format(PARAMS.dx / UNITS.MeV)

        self.step = 0

        PARAMS.update(self.total_energy_density())
        self.data = {
            'aT': array.array('f', [PARAMS.aT]),
            'T': array.array('f', [PARAMS.T]),
            'a': array.array('f', [PARAMS.a]),
            'x': array.array('f', [PARAMS.x]),
            't': array.array('f', [PARAMS.t]),
            'rho': array.array('f', [PARAMS.rho])
        }

        print '#step\tTime, s\taT, MeV\tT, MeV\tscale factor\tρ energy density, eV^4\tH, GeV'
        self.log()
        self.step += 1

        while PARAMS.T > PARAMS.T_final:
            try:
                self.make_step()
            except KeyboardInterrupt:
                print "Keyboard interrupt!"
                break

        self.log()
        for particle in self.particles:
            print particle

        print "Data saved to file {}".format(self.logger.log.name)

        """
        == Postmortem debugger ==
        Invoking `pdb` debugger after the computation is finished allows to do final adjustments \
        and export the missing data without restarting the whole simulation.
        """

        if self.postmortem_debugger:
            import ipdb
            ipdb.set_trace()

        return self.data

    def make_step(self):
        integrator = integrators.heun_correction if self.INTEGRATION_METHOD == 'heun' \
            else integrators.euler_correction

        PARAMS.aT += integrator(y=PARAMS.aT, f=self.integrand, t=PARAMS.x, h=PARAMS.dx)
        PARAMS.x += PARAMS.dx

        PARAMS.update(self.total_energy_density())

        self.save()

        self.log()
        self.step += 1

    def update_particles(self):
        """ === 1. Update particles state ===
            Update particle species distribution functions, check for regime switching,\
            update precalculated variables like energy density and pressure. """
        for particle in self.particles:
            particle.update()

    def init_interactions(self):
        """ === 2. Initialize non-equilibrium interactions ===
            Non-equilibrium interactions of different particle species are treated by a\
            numerical integration of the Boltzmann equation for distribution functions in\
            the expanding space-time.

            Depending on the regime of the particle species involved and cosmological parameters, \
            each `Interaction` object populates `Particle.collision_integrands` array with \
            currently active `Integral` objects.
        """
        for interaction in self.interactions:
            interaction.calculate()

    def calculate_collisions(self):
        """ === 3. Calculate collision integrals === """
        for particle in self.particles:
            if not particle.collision_integrands:
                continue

            if self.PARALLELIZE:
                particle.collision_integral = \
                    numpy.array(parallelization.parmap(particle.integrate_collisions,
                                                       GRID.TEMPLATE))
            else:
                particle.collision_integral = \
                    numpy.vectorize(particle.integrate_collisions,
                                    otypes=[numpy.float_])(GRID.TEMPLATE)

    def control_step_size(self, maximum_change=0.2, fallback_change=0.1):
        """
        === 4. Control integration step size ===

        :param maximum_change: Maximal allowed relative difference in the distribution functions\
                               value after Boltzmann equation integration. Generally, master\
                               equation of the temperature can be solved with arbitrarily large\
                               step size without significant loss in accuracy. The main problem is\
                               the truncation error of the Boltzmann equation integration.

        :param fallback_change: In the case when distribution functions change too fast, step size\
                                is urgently decreased to compensate this effect. `fallback_change`\
                                smaller than `maximum_change` enhances the integration stability.

        :param exponent: If distribution functions change is insignificant, step size can be\
                         increased by a factor of `exponent`.
        """

        if not self.adaptive_step_size:
            return

        exponent = self.step_size_multiplier

        dx = PARAMS.dx
        multipliers = []
        for particle in self.particles:
            if particle.in_equilibrium:
                continue
            relative_delta = numpy.absolute(particle.collision_integral * dx
                                            / particle._distribution).max()
            if relative_delta > maximum_change:
                multipliers.append(fallback_change / relative_delta)
            else:
                multipliers.append(exponent)

        multiplier = min(multipliers) if multipliers else 1.

        if multiplier != 1.:
            PARAMS.dx *= multiplier
            print "//// Step size changed to dx =", PARAMS.dx / UNITS.MeV

    def update_distributions(self):
        """ === 5. Update particles distributions === """

        for particle in self.particles:
            particle.update_distribution()

    def calculate_temperature_terms(self):
        """ === 6. Calculate temperature equation terms === """

        numerator = 0
        denominator = 0

        for particle in self.particles:
            numerator += particle.numerator()
            denominator += particle.denominator()

        return numerator, denominator

    def integrand(self, t, y):
        """ == Temperature equation integrand ==

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

        # from common import PARAMS
        # old_PARAMS = PARAMS.copy()
        # old_particles = copy.copy(self.particles)
        # PARAMS.aT = y
        # PARAMS.x = t
        # PARAMS.update(self.total_energy_density())

        self.update_particles()

        self.init_interactions()

        self.calculate_collisions()

        self.control_step_size()

        self.update_distributions()

        numerator, denominator = self.calculate_temperature_terms()

        # PARAMS = old_PARAMS
        # self.particles = old_particles

        return numerator / denominator

    def save(self):
        """ Save current Universe parameters into the data arrays or output files """
        self.data['aT'].append(PARAMS.aT)
        self.data['T'].append(PARAMS.T)
        self.data['a'].append(PARAMS.a)
        self.data['x'].append(PARAMS.x)
        self.data['rho'].append(PARAMS.rho)
        self.data['t'].append(PARAMS.t)

    def log(self):
        """ Runtime log output """

        # Print parameters every now and then
        if self.step % self.log_freq == 0:
            print '#' + str(self.step), \
                '\tt =', PARAMS.t / UNITS.s, \
                '\taT =', PARAMS.aT / UNITS.MeV, \
                '\tT =', PARAMS.T / UNITS.MeV, \
                '\ta =', PARAMS.a

            if self.graphics:
                self.graphics.plot(self.data)

    def total_energy_density(self):
        return sum(particle.energy_density() for particle in self.particles)
