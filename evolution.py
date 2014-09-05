# -*- coding: utf-8 -*-

import numpy
import numericalunits as nu
import array

from common import UNITS, PARAMS, CONST, parmap
from plotting import Plotting


class Universe:

    """ The master object that governs the calculation. """

    def __init__(self, particles=[], interactions=[]):
        """
        :param particles: List of `particle.Particle` objects - set of particles to model
        :param interactions: List of `interaction.Interaction` objects - quantum interactions \
                             between particle species
        """
        self.particles = particles
        self.interactions = interactions
        self.graphics = Plotting()

    def evolve(self, T_final=PARAMS.T_final, dx=PARAMS.dx):
        """
        == Main computing routine ==

        Modeling is carried in steps of scale factor, starting from the initial conditions defined\
        by a single parameter: the initial temperature.

        Initial temperature (e.g. 10 MeV) corresponds to a point in the history of the Universe \
        long before BBN when most particle species are in the thermodynamical equilibrium

        """

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

        while PARAMS.T > T_final:

            rho = 0  # total energy density
            numerator = 0
            denominator = 0

            for particle in self.particles:
                rho += particle.energy_density()
                numerator += particle.numerator()
                denominator += particle.denominator()

            """ Conformal scale factor $x = a m$ step size is fixed: expansion of the Universe \
                is the main factor that controls the thermodynamical state of the system """
            PARAMS.x += dx
            """ TODO: Sanity check: $a T$ remains constant if the entropy of the system is \
                conserved.

                A bump in the $a T$ can be seen if the number of relativistic degrees of\
                freedom changes as a result of effective disappearance of the particle species \
                from the system. For example, electron-positron pairs annihilate at temperatures\
                close to the electron mass $\sim 511 keV$. Basically, this result in a cosmic \
                photon and neutrino backgrounds temperature ratio around \
                $\frac{T_\nu}{T_\gamma} \sim 1.401$"""
            PARAMS.aT += numerator / denominator * dx
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

            self.save(rho)

            """ Non equilibrium interactions of different particle species are treated by a\
                numerical integration of the Boltzman equation for distribution functions in\
                the expanding spacetime.
            """
            for interaction in self.interactions:
                interaction.calculate()

            """ Update particle species distribution functions, check for regime switching, update\
                precalculated variables like energy density and pressure. """
            for particle in self.particles:
                particle.update()

            self.log(rho)
            self.step += 1

        return self.data

    def save(self, rho):
        """ Save current Universe parameters into the data arrays or output files """
        self.data['aT'].append(PARAMS.aT)
        self.data['T'].append(PARAMS.T)
        self.data['a'].append(PARAMS.a)
        self.data['rho'].append(rho)
        self.data['t'].append(PARAMS.t)

    def log(self, rho):
        """ Runtime log output """

        # Print parameters each 100 steps
        if self.step % 100 == 0:
            print 't =', PARAMS.t / UNITS.s, \
                '\taT =', PARAMS.aT / UNITS.MeV, \
                "\tT =", PARAMS.T / UNITS.MeV, \
                "\ta =", PARAMS.a, \
                "\tρ =", PARAMS.rho / UNITS.eV**4, \
                "\tH =", PARAMS.H / nu.GeV
            self.graphics.plot(self.data)

        # Clean up the data array to free the memory
        # if self.step % 1000 == 0:
        #     for k in self.data.keys():
        #         self.data[k] = self.data[k][-1000:]
        #     self.graphics.plot(self.data, full=True)

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
