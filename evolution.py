# -*- coding: utf-8 -*-

import numpy
import numericalunits as nu
import matplotlib.pyplot as plt
import array
from common import UNITS, PARAMS, CONST, GRID, benchmark, parmap


class Plotting:
    particles = None

    def __init__(self):
        self.params_figure, self.plots = plt.subplots(2, 2, num=1)
        self.params_figure.subplots_adjust(hspace=0.5, wspace=0.5)

        self.plots[0][0].set_title("Temperature")
        self.plots[0][0].set_xlabel("time, s")
        self.plots[0][0].set_ylabel("T, MeV")
        self.plots[0][0].set_yscale("log")

        self.plots[0][1].set_title("Scale factor")
        self.plots[0][1].set_xlabel("time, s")
        self.plots[0][1].set_ylabel("a, 1")

        self.plots[1][0].set_title("T * a")
        self.plots[1][0].set_xlabel("time, s")
        self.plots[1][0].set_ylabel("T * a, MeV")

        self.plots[1][1].set_title("Total energy density")
        self.plots[1][1].set_yscale("log")
        self.plots[1][1].set_xlabel("time, s")
        self.plots[1][1].set_ylabel("rho, eV**4")

        self.params_figure.show()

    def monitor(self, particles=[]):
        self.particles = particles
        if self.particles:
            self.particles_figure, self.particles_plots = plt.subplots(len(particles)*2, 1, num=2)
            self.particles_figure.subplots_adjust(hspace=0.5, wspace=0.5)

            # if len(self.particles) == 1:
                # self.particles_plots = [self.particles_plots]
            for i, particle in enumerate(self.particles):
                self.particles_plots[i*2].set_title(particle.name)
                self.particles_plots[i*2].set_xlabel("y, MeV")
                self.particles_plots[i*2].set_ylabel("f")

                self.particles_plots[i*2 + 1].set_title(particle.name)
                self.particles_plots[i*2 + 1].set_xlabel("y, MeV")
                self.particles_plots[i*2 + 1].set_ylabel("f/f_eq")

            self.particles_figure.show()

    def plot(self, data, full=False):
        if not full:
            ts = data['t'][-1]
            Ts = data['T'][-1]
            As = data['a'][-1]
            rhos = data['rho'][-1]
            style = 'b.'
        else:
            ts = numpy.array(data['t'])
            Ts = numpy.array(data['T'])
            As = numpy.array(data['a'])
            rhos = numpy.array(data['rho'])
            style = 'b-'

        self.plots[0][0].plot(ts / UNITS.s, Ts / nu.MeV, style)
        self.plots[0][1].plot(ts / UNITS.s, As, style)
        self.plots[1][0].plot(ts / UNITS.s, Ts * As / nu.MeV, style)
        self.plots[1][1].plot(ts / UNITS.s, rhos / nu.eV**4, style)

        plt.figure(1)
        plt.draw()

        if self.particles:
            for i, particle in enumerate(self.particles):
                if not particle.in_equilibrium:
                    self.particles_plots[i*2].plot(GRID.TEMPLATE / nu.MeV, particle._distribution)
                    feq = particle.distribution_function_vectorized(
                        particle.energy_normalized_vectorized(GRID.TEMPLATE)
                        / particle.temperature
                    )
                    self.particles_plots[i*2 + 1].plot(
                        GRID.TEMPLATE / nu.MeV,
                        particle._distribution / feq
                    )

            plt.figure(2)
            plt.draw()


class Universe:
    def __init__(self, particles=[], interactions=[]):
        self.particles = particles
        self.interactions = interactions
        self.graphics = Plotting()

    def totals(self):
        totals = parmap(lambda particle: (
            particle.energy_density(),
            particle.numerator(),
            particle.denominator()
        ), self.particles)

        rho = 0
        numerator = 0
        denominator = 0
        for particle in totals:
            rho += particle[0]
            numerator += particle[1]
            denominator += particle[2]

        return rho, numerator, denominator

    def evolve(self, T_final=PARAMS.T_final, dx=PARAMS.dx):

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

            # if PARAMS.T < 2. * nu.MeV:
            # with benchmark("Parallel all"):
            # rho, numerator, denominator = self.totals()

            # with benchmark("energy_density"):
            #     rho = 0
            #     for particle in self.particles:
            #         rho += particle.energy_density()

            for particle in self.particles:
                rho += particle.energy_density()
                numerator += particle.numerator()
                denominator += particle.denominator()

            PARAMS.aT += numerator / denominator * dx
            PARAMS.x += dx
            PARAMS.a = PARAMS.x / PARAMS.m
            PARAMS.T = PARAMS.aT / PARAMS.a
            PARAMS.H = numpy.sqrt(8./3.*numpy.pi * CONST.G * rho)

            self.data['aT'].append(PARAMS.aT)
            self.data['T'].append(PARAMS.T)
            self.data['a'].append(PARAMS.a)
            self.data['rho'].append(rho)

            PARAMS.t += (self.data['a'][-1] / self.data['a'][-2] - 1) / PARAMS.H

            self.data['t'].append(PARAMS.t)

            if self.step % 1 == 0:
                self.graphics.plot(self.data)

            for interaction in self.interactions:
                interaction.calculate()

            for particle in self.particles:
                particle.update()

            if self.step % 1 == 0:
                print 't =', PARAMS.t / UNITS.s, \
                    "\tT =", PARAMS.T / nu.MeV, \
                    "\ta =", PARAMS.a, \
                    "\trho =", self.data['rho'][-1] / nu.MeV**4, \
                    "\tH =", PARAMS.H / nu.GeV
                # self.graphics.plot(self.data)

            self.step += 1

        self.graphics.plot(self.data, full=True)

        return self.data
