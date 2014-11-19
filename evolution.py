# -*- coding: utf-8 -*-

import sys
import copy
import numpy
import numericalunits as nu
import array
from datetime import datetime

from common import UNITS, PARAMS, GRID
from common import integrators, parallelization, utils
from plotting import Plotting


class Universe:

    """ The master object that governs the calculation. """

    log_freq = 1
    INTEGRATION_METHOD = ['euler', 'heun'][0]
    PARALLELIZE = True

    def __init__(self, particles=[], interactions=[],
                 logfile='logs/' + str(datetime.now()) + '.txt',
                 plotting=True, adaptive_step_size=True):
        """
        :param particles: List of `particle.Particle` objects - set of particles to model
        :param interactions: List of `interaction.Interaction` objects - quantum interactions \
                             between particle species
        """
        self.particles = particles
        self.interactions = interactions

        self.graphics = Plotting() if plotting else None
        self.adaptive_step_size = adaptive_step_size

        sys.stdout = utils.Logger(logfile)
        self.logger = sys.stdout

    def evolve(self):
        """
        == Main computing routine ==

        Modeling is carried in steps of scale factor, starting from the initial conditions defined\
        by a single parameter: the initial temperature.

        Initial temperature (e.g. 10 MeV) corresponds to a point in the history of the Universe \
        long before BBN when most particle species are in the thermodynamical equilibrium

        """

        print "dx =", PARAMS.dx

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
                if self.INTEGRATION_METHOD == 'heun':
                    PARAMS.aT += integrators.heun_correction(y=PARAMS.aT, f=self.integrand,
                                                             t=PARAMS.x, h=PARAMS.dx)

                elif self.INTEGRATION_METHOD == 'euler':
                    PARAMS.aT += integrators.euler_correction(y=PARAMS.aT, f=self.integrand,
                                                              t=PARAMS.x, h=PARAMS.dx)

                PARAMS.x += PARAMS.dx

                PARAMS.update(self.total_energy_density())

                self.save()

                self.log()
                self.step += 1
            except KeyboardInterrupt:
                print "Keyboard interrupt!"
                break

        self.log()
        for particle in self.particles:
            print particle

        print "Data saved to file {}".format(self.logger.log.name)
        return self.data

    def integrand(self, t, y):
        """ Master equation for the temperature looks like

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

        from common import PARAMS
        old_PARAMS = copy.copy(PARAMS)
        old_particles = copy.copy(self.particles)
        PARAMS.aT = y
        PARAMS.x = t
        PARAMS.update(self.total_energy_density())

        """ Update particle species distribution functions, check for regime switching,\
            update precalculated variables like energy density and pressure. """
        for particle in self.particles:
            particle.update()

        """ Non equilibrium interactions of different particle species are treated by a\
            numerical integration of the Boltzman equation for distribution functions in\
            the expanding spacetime.
        """
        for interaction in self.interactions:
            interaction.calculate()

        for particle in self.particles:
            # For non-equilibrium particle calculate the collision integrals
            if particle.collision_integrands:
                if self.PARALLELIZE:
                    particle.collision_integral = \
                        numpy.array(parallelization.parmap(particle.integrate_collisions,
                                                           GRID.TEMPLATE, workers=7))
                else:
                    particle.collision_integral = \
                        numpy.vectorize(particle.integrate_collisions,
                                        otypes=[numpy.float_])(GRID.TEMPLATE)

        self.control_step_size()

        for particle in self.particles:
            particle.update_distribution()

        numerator = 0
        denominator = 0

        for particle in self.particles:
            numerator += particle.numerator()
            denominator += particle.denominator()

        PARAMS = old_PARAMS
        self.particles = old_particles

        return numerator / denominator

    def control_step_size(self):
        """
        Scale factor step size controller.

        TODO: not usable now, need to orchestrate between many collision integrations - \
              i.e., select smallest suggested
        """

        if not self.adaptive_step_size:
            return

        dx = PARAMS.dx
        multipliers = []
        for particle in self.particles:
            if particle.in_equilibrium:
                continue
            relative_delta = numpy.absolute(particle.collision_integral * dx
                                            / particle._distribution).max()
            if relative_delta > 0.2:
                multipliers.append(0.1 / relative_delta)
            else:
                multipliers.append(1.25)

        multiplier = min(multipliers) if multipliers else 1.

        if multiplier != 1.:
            PARAMS.dx *= multiplier
            print "//// Step size changed to dx =", PARAMS.dx / UNITS.MeV

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
                # self.graphics.plot_pipe.send([self.data, False])

    def total_energy_density(self):
        # return sum(parallelization.parmap(lambda particle: particle.energy_density(),
        #                                   self.particles))
        return sum(particle.energy_density() for particle in self.particles)
