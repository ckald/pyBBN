# -*- coding: utf-8 -*-

import sys
import copy
import numpy
import array
from datetime import datetime

from common import UNITS, PARAMS, GRID, CONST
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

        self.logfile = logfile
        self.init_log()

    def evolve(self):
        """
        == Main computing routine ==

        Modeling is carried in steps of scale factor, starting from the initial conditions defined\
        by a single parameter: the initial temperature.

        Initial temperature (e.g. 10 MeV) corresponds to a point in the history of the Universe \
        long before then BBN. Then most particle species are in the thermodynamical equilibrium.

        """

        for particle in self.particles:
            print particle

        for interaction in self.interactions:
            print interaction

        print "dx = {} MeV".format(PARAMS.dx / UNITS.MeV)

        self.step = 0

        PARAMS.update(self.total_energy_density())
        self.data = {
            'aT': array.array('f', [PARAMS.aT]),
            'T': array.array('f', [PARAMS.T]),
            'a': array.array('f', [PARAMS.a]),
            'x': array.array('f', [PARAMS.x]),
            't': array.array('f', [PARAMS.t]),
            'rho': array.array('f', [PARAMS.rho]),
            'fraction': array.array('f', [0]),
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

        print "Data saved to file {}".format(self.logfile)

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
        PARAMS.dx = PARAMS.x * (numpy.exp(PARAMS.dy) - 1.)
        self.integrand(PARAMS.x, PARAMS.aT)

        order = min(self.step + 1, 5)
        fs = self.data['fraction'][-order+1:]
        fs.append(self.fraction)

        PARAMS.aT += integrators.adams_bashforth_correction(fs=fs, h=PARAMS.dy, order=order)
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
            each `Interaction` object populates `Particle.collision_integrals` array with \
            currently active `Integral` objects.
        """
        for interaction in self.interactions:
            interaction.initialize()

    def calculate_collisions(self):
        """ === 3. Calculate collision integrals === """
        for particle in self.particles:
            if not particle.collision_integrals:
                continue

            if self.PARALLELIZE:
                particle.collision_integral = \
                    numpy.array(parallelization.parmap(particle.integrate_collisions,
                                                       GRID.TEMPLATE))
            else:
                particle.collision_integral = \
                    numpy.vectorize(particle.integrate_collisions,
                                    otypes=[numpy.float_])(GRID.TEMPLATE)

            print particle.symbol, "I =", particle.collision_integral * UNITS.MeV

    def update_distributions(self):
        """ === 4. Update particles distributions === """

        for particle in self.particles:
            particle.update_distribution()

    def calculate_temperature_terms(self):
        """ === 5. Calculate temperature equation terms === """

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

        """
        Save system state (unfinished, inactive)

        ```python
            from common import PARAMS
            old_PARAMS = PARAMS.copy()
            old_particles = copy.copy(self.particles)
            PARAMS.aT = y
            PARAMS.x = t
            PARAMS.update(self.total_energy_density())
        ```
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
        self.fraction = PARAMS.x * numerator / denominator

        """
        Load system state (unfinished, inactive)

        ```python
        PARAMS = old_PARAMS
        self.particles = old_particles
        ```
        """

        return self.fraction

    def save(self):
        """ Save current Universe parameters into the data arrays or output files """
        self.data['aT'].append(PARAMS.aT)
        self.data['T'].append(PARAMS.T)
        self.data['a'].append(PARAMS.a)
        self.data['x'].append(PARAMS.x)
        self.data['rho'].append(PARAMS.rho)
        self.data['t'].append(PARAMS.t)
        self.data['fraction'].append(self.fraction)

    def init_log(self):
        sys.stdout = utils.Logger(self.logfile)

    def log(self):
        """ Runtime log output """

        # Print parameters every now and then
        if self.step % self.log_freq == 0:
            print '#' + str(self.step), \
                '\tt =', PARAMS.t / UNITS.s, \
                '\taT =', PARAMS.aT / UNITS.MeV, \
                '\tT =', PARAMS.T / UNITS.MeV, \
                '\ta =', PARAMS.a, \
                '\tdx =', PARAMS.dx / UNITS.MeV

            if self.graphics:
                self.graphics.plot(self.data)

    def total_energy_density(self):
        return sum(particle.energy_density() for particle in self.particles)
