# -*- coding: utf-8 -*-

import sys
import numpy
import numericalunits as nu
import array
from datetime import datetime

from common import UNITS, PARAMS, CONST, GRID, parmap, Logger, forward_euler_integrator
from plotting import Plotting


class Universe:

    """ The master object that governs the calculation. """

    log_freq = 1

    def __init__(self, particles=[], interactions=[],
                 logfile='logs/' + str(datetime.now()) + '.txt'):
        """
        :param particles: List of `particle.Particle` objects - set of particles to model
        :param interactions: List of `interaction.Interaction` objects - quantum interactions \
                             between particle species
        """
        self.particles = particles

        self.interactions = interactions
        self.graphics = Plotting()

        sys.stdout = Logger(logfile)
        self.logger = sys.stdout

        for particle in self.particles:
            print particle

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

        self.data = {
            'aT': array.array('f', [PARAMS.aT]),
            'T': array.array('f', [PARAMS.T]),
            'a': array.array('f', [PARAMS.a]),
            't': array.array('f', [PARAMS.t]),
            'rho': array.array('f', [0])
        }
        for particle in self.particles:
            self.data['rho'][0] += particle.energy_density()

        print '# Time, s\taT, MeV\tTemperature, MeV\tscale factor\tρ energy density, eV^4\tH, GeV'
        self.log()

        def integrand(t, y):

            numerator = 0
            denominator = 0

            for particle in self.particles:
                numerator += particle.numerator()
                denominator += particle.denominator()

            return numerator / denominator

        while PARAMS.T > PARAMS.T_final:

            rho = 0

            for particle in self.particles:
                rho += particle.energy_density()

            """ Conformal scale factor $x = a m$ step size is fixed: expansion of the Universe \
                is the main factor that controls the thermodynamical state of the system """
            PARAMS.x += PARAMS.dx

            PARAMS.aT = forward_euler_integrator(y=self.data['aT'],
                                                 t=self.data['a'],
                                                 f=integrand,
                                                 h=PARAMS.dx)
            """ Physical scale factor and temperature for convenience """
            PARAMS.a = PARAMS.x / PARAMS.m
            PARAMS.T = PARAMS.aT / PARAMS.a
            """ Hubble expansion parameter defined by a Friedmann equation:

                \begin{equation}
                    H = \sqrt{\frac{8 \pi}{3} G \rho}
                \end{equation}
            """
            PARAMS.H = numpy.sqrt(8./3.*numpy.pi * CONST.G * rho)

            """ Time step size is inferred from the approximation of the scale factor `a` \
                derivative and a definition of the Hubble parameter `H`:

                \begin{equation}
                    H = \frac{\dot{a}}{a} = \frac{1}{a_{i-1}} \frac{a_i - a_{i-1}}{\Delta t} \
                      = \frac{1}{\Delta t}(\frac{a_i}{a_{i-1}} -1)
                \end{equation}

                \begin{equation}
                    \Delta t = \frac{1}{H} (\frac{a_i}{a_{i-1}} - 1)
                \end{equation}
            """
            dt = (PARAMS.a / self.data['a'][-1] - 1) / PARAMS.H
            PARAMS.t += dt
            PARAMS.rho = rho

            self.save()

            """ Non equilibrium interactions of different particle species are treated by a\
                numerical integration of the Boltzman equation for distribution functions in\
                the expanding spacetime.
            """
            for interaction in self.interactions:
                interaction.calculate()

            for particle in self.particles:
                # For non-equilibrium particle calculate the collision integrals
                particle.collision_integral = \
                    particle.integrate_collisions_vectorized(GRID.TEMPLATE)

            """ Update particle species distribution functions, check for regime switching, update\
                precalculated variables like energy density and pressure. """
            for particle in self.particles:
                particle.update()

            for particle in self.particles:
                particle.update_distribution()

            self.log()
            self.step += 1

        self.log()
        for particle in self.particles:
            print particle

        print "Data saved to file {}".format(self.logger.log.name)
        return self.data

    def save(self):
        """ Save current Universe parameters into the data arrays or output files """
        self.data['aT'].append(PARAMS.aT)
        self.data['T'].append(PARAMS.T)
        self.data['a'].append(PARAMS.a)
        self.data['rho'].append(PARAMS.rho)
        self.data['t'].append(PARAMS.t)

    def log(self):
        """ Runtime log output """

        # Print parameters every now and then
        if self.step % self.log_freq == 0:
            print 't =', PARAMS.t / UNITS.s, \
                '\taT =', PARAMS.aT / UNITS.MeV, \
                "\tT =", PARAMS.T / UNITS.MeV, \
                "\ta =", PARAMS.a, \
                "\tρ =", PARAMS.rho / UNITS.eV**4, \
                "\tH =", PARAMS.H / nu.GeV

            self.graphics.plot(self.data)
            # self.graphics.plot_pipe.send([self.data, False])

    def totals(self):
        """
        Parallel realization of the Friedmann equation parameters computation.
        """
        particle_totals = parmap(lambda particle: (
            particle.energy_density(),
            particle.numerator(),
            particle.denominator()
        ), self.particles)

        rho = 0
        numerator = 0
        denominator = 0

        for particle in particle_totals:
            rho += particle[0]
            numerator += particle[1]
            denominator += particle[2]

        return rho, numerator, denominator
