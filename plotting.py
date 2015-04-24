# -*- coding: utf-8 -*-

import os
import itertools

import numpy
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from common import UNITS, GRID
from common.utils import ring_deque


class Plotting(object):
    particles = None

    def __init__(self):
        """ Initialize plots, setup basic styling """

        plt.rcParams['toolbar'] = 'None'

        self.params_figure, self.plots = plt.subplots(2, 2, num=1)
        self.plots = [p for p in itertools.chain(*self.plots)]
        self.params_figure.subplots_adjust(hspace=0.5, wspace=0.5)

        self.plot_map = ['T', 'a', 'aT', 'rho']
        self.divider_map = [UNITS.MeV, 1, UNITS.MeV, UNITS.eV**4]

        self.plots[0].set_title("Temperature")
        self.plots[0].set_xlabel("time, s")
        self.plots[0].set_ylabel("T, MeV")
        self.plots[0].set_yscale("log")

        self.plots[1].set_title("Scale factor")
        self.plots[1].set_xlabel("time, s")
        self.plots[1].set_ylabel("a, 1")
        self.plots[1].set_ylim(0, 1)

        self.plots[2].set_title("T * a")
        self.plots[2].set_xlabel("time, s")
        self.plots[2].set_xscale("log")
        self.plots[2].set_ylabel("T * a, MeV")
        self.plots[2].set_ylim(1, 1.1)

        self.plots[3].set_title("Total energy density")
        self.plots[3].set_yscale("log")
        self.plots[3].set_xlabel("time, s")
        self.plots[3].set_ylabel(u"ρ, eV**4")

        self.lines = []
        self.plots_data = []
        self.times = ring_deque([], 1000)
        for plot in self.plots:
            self.lines.append(plot.plot([], [], 'b-')[0])
            self.plots_data.append(ring_deque([], 1000))

        self.params_figure.show()

    def save(self, filename):
        """ Save cosmological and monitored particles plots to the file in the same folder as \
            `filename` """

        folder = os.path.split(filename)[0]
        plt.figure(1)
        plt.savefig(os.path.join(folder, 'plots.png'))
        if self.particles:
            plt.figure(2)
            plt.savefig(os.path.join(folder, 'particles.png'))

    def monitor(self, map):
        """ Setup the detailed distribution function and energy density plots for specific \
            particle species """

        self.particles_figure, self.particles_plots = plt.subplots(len(map), 2, num=2)
        self.particles_figure.subplots_adjust(hspace=0.5, wspace=0.5)

        self.particles = self.particles if self.particles else []

        for i, (particle, monitor) in enumerate(map):
            self.particles.append((particle, monitor(particle, self.particles_plots[i])))

        self.particles_figure.show()

    def plot(self, data):
        """ Plot cosmological parameters and monitored particles distribution functions """

        last_t = data['t'][-1] / UNITS.s
        self.times.append(last_t)

        for i, plot in enumerate(self.plots):
            _, xmax = plot.get_xlim()
            ymin, ymax = plot.get_ylim()

            if last_t >= xmax:
                plot.set_xlim(self.times[0], last_t * 1.1)

            last_data = data[self.plot_map[i]][-1] / self.divider_map[i]
            self.plots_data[i].append(last_data)

            if last_data >= ymax:
                plot.set_ylim(self.plots_data[i].min, 1.1 * last_data)
            if last_data <= ymin:
                plot.set_ylim(last_data / 1.1, self.plots_data[i].max)

            self.lines[i].set_data(self.times, self.plots_data[i])

        plt.figure(1)
        plt.draw()

        if self.particles:
            for i, (particle, monitor) in enumerate(self.particles):
                monitor.plot(data)

            plt.figure(2)
            plt.draw()


class ParticleMonitor(object):
    def __init__(self, particle, plots):
        raise NotImplementedError()

    def plot(self, data):
        raise NotImplementedError()


class RadiationParticleMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].set_xlabel("a")
        self.plots[0].set_xscale("log")
        self.plots[0].set_ylabel("rho/rho_eq")

        self.plots[1].set_xlabel("y, MeV")
        self.plots[1].set_xlim(GRID.MIN_MOMENTUM / UNITS.MeV, GRID.MAX_MOMENTUM / UNITS.MeV)
        self.plots[1].set_ylabel("f/f_eq")

    def comparison_distributions(self, data):
        a = self.particle.params.a
        aT = self.particle.aT

        rhoeq = self.particle.energy_density() / (
            7. * self.particle.dof * numpy.pi**2
            * (self.particle.params.m / a)**4 / 240.
        )
        feq = self.particle.equilibrium_distribution(aT=aT)

        return (a, rhoeq), feq

    def plot(self, data):
        (a, rhoeq), feq = self.comparison_distributions(data)

        self.plots[0].scatter(a, rhoeq, s=1)

        age_lines(self.plots[1].get_axes().lines)

        self.plots[1].plot(
            GRID.TEMPLATE / UNITS.MeV,
            numpy.vectorize(self.particle.distribution)(GRID.TEMPLATE) / feq
        )


class EquilibriumRadiationParticleMonitor(RadiationParticleMonitor):
    def comparison_distributions(self, data):
        a = data['a'][-1]
        aT = data['aT'][-1]

        rhoeq = self.particle.energy_density() / (
            7. * self.particle.dof * numpy.pi**2
            * (self.particle.params.m / a)**4 / 240.
        )
        feq = self.particle.equilibrium_distribution(aT=aT)

        return (a, rhoeq), feq


class MassiveParticleMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].set_xlabel("a")
        self.plots[0].set_xscale("log")
        self.plots[0].set_ylabel("rho/(n M)")

        self.plots[1].set_xlabel("y, MeV")
        self.plots[1].set_xlim(GRID.MIN_MOMENTUM / UNITS.MeV, GRID.MAX_MOMENTUM / UNITS.MeV)
        self.plots[1].set_ylabel("(f-f_eq) y^2")

    def plot(self, data):
        a = data['a'][-1]

        from particles.NonEqParticle import energy_density, density

        self.plots[0].scatter(a, energy_density(self.particle)
                              / (self.particle.mass * density(self.particle)), s=1)

        age_lines(self.plots[1].get_axes().lines)

        yy = GRID.TEMPLATE * GRID.TEMPLATE / UNITS.MeV**2

        f = self.particle._distribution
        feq = self.particle.equilibrium_distribution()
        self.plots[1].plot(GRID.TEMPLATE / UNITS.MeV, yy*(f-feq))


class EquilibrationMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].set_xlabel("a")
        self.plots[0].set_xscale("log")
        self.plots[0].set_ylabel("|I|")

        self.plots[1].set_xlabel("a")
        self.plots[1].set_ylabel("numerator")

    def plot(self, data):
        a = data['a'][-1]

        from particles.NonEqParticle import numerator

        self.plots[0].scatter(a, numpy.max(numpy.fabs(self.particle.collision_integral)), s=1)
        self.plots[1].scatter(a, numerator(self.particle), s=1)


def age_lines(lines):
    """ Slightly decrease the opacity of plotted lines until they are barely visible.\
        Then, remove them. Saves up on memory and clears the view of the plots. """
    for line in lines:
        alpha = line.get_alpha() or 1.
        if alpha < 0.1:
            line.remove()
        else:
            line.set_alpha((line.get_alpha() or 1.) * 0.8)


def plot_integrand(integrand, name, p0, filename=__file__):
    """ Save a 3D plot of the distribution function integrand into a file. """
    fig = plt.figure(3)
    ax = fig.gca(projection='3d')
    plt.cla()
    X, Y = numpy.meshgrid(GRID.TEMPLATE, GRID.TEMPLATE)
    Z = numpy.array([integrand([x, y]) for x, y in zip(numpy.ravel(X), numpy.ravel(Y))])\
        .reshape(X.shape)

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.1)
    ax.contourf(X, Y, Z, zdir='z', offset=numpy.amin(Z), cmap=cm.coolwarm)
    ax.contourf(X, Y, Z, zdir='x', offset=ax.get_xlim()[0], cmap=cm.coolwarm)
    ax.contourf(X, Y, Z, zdir='y', offset=ax.get_ylim()[1], cmap=cm.coolwarm)

    ax.set_xlabel('p1')
    ax.set_ylabel('p2')
    ax.set_title('{} p0 = {}'.format(name, p0 / UNITS.MeV))

    plt.savefig(os.path.join(os.path.split(filename)[0], 'logs/plt_{}.png'.format(p0 / UNITS.MeV)))


def plot_points(points, name):
    """ Draw a scatter plot for a number of `points` tuples `(x, y)` """
    plt.figure(4)
    plt.title(name)
    plt.scatter(*zip(*points))
    plt.show()
