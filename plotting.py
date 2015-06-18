# -*- coding: utf-8 -*-

import os
import itertools
import math

import numpy
import matplotlib.pyplot as plt

from common import UNITS, GRID
from common.utils import ring_deque, getboolenv
from particles import STATISTICS


class Plotting(object):
    particles = None

    def __init__(self, show=True):
        """ Initialize plots, setup basic styling """

        self.show_plots = getboolenv("SHOW_PLOTS", show)

        self.params_figure, self.plots = plt.subplots(2, 2, num=1)
        self.plots = list(itertools.chain(*self.plots))
        self.params_figure.subplots_adjust(hspace=0.5, wspace=0.5)

        self.plot_map = ['T', 'a', 'aT', 'N_eff']
        self.divider_map = [UNITS.MeV, 1, UNITS.MeV, 1]

        self.plots[0].set_title("Temperature")
        self.plots[0].set_xlabel("time, s")
        self.plots[0].set_xscale("log")
        self.plots[0].set_yscale("log")
        self.plots[0].set_ylabel("T, MeV")
        self.plots[0].set_yscale("log")

        self.plots[1].set_title("Scale factor")
        self.plots[1].set_xlabel("time, s")
        self.plots[1].set_xscale("log")
        self.plots[1].set_yscale("log")
        self.plots[1].set_ylabel("a, 1")
        self.plots[1].set_ylim(0, 1)

        self.plots[2].set_title("T * a")
        self.plots[2].set_xlabel("time, s")
        self.plots[2].set_xscale("log")
        self.plots[2].set_ylabel("T * a, MeV")
        self.plots[2].set_ylim(1, 1.1)

        self.plots[3].set_title("N_eff")
        self.plots[3].set_xlabel("time, s")
        self.plots[3].invert_xaxis()
        self.plots[3].set_xscale("log")
        self.plots[3].set_ylabel("N_eff")

        self.lines = []
        self.plots_data = []
        self.times = ring_deque([], 1000)
        for plot in self.plots:
            self.lines.append(plot.plot([], [], 'b-')[0])
            self.plots_data.append(ring_deque([], 1000))

        if self.show_plots:
            self.params_figure.show()

    def save(self, filename):
        """ Save cosmological and monitored particles plots to the file in the same folder as \
            `filename` """

        folder = os.path.split(filename)[0]
        self.params_figure.savefig(os.path.join(folder, 'plots.svg'))
        if self.particles:
            self.particles_figure.savefig(os.path.join(folder, 'particles.svg'))

    def monitor(self, map):
        """ Setup the detailed distribution function and energy density plots for specific \
            particle species """

        self.particles_figure, self.particles_plots = plt.subplots(len(map), 2, num=2)
        self.particles_figure.subplots_adjust(hspace=0.5, wspace=0.5)

        self.particles = self.particles if self.particles else []

        for i, (particle, monitor) in enumerate(map):
            self.particles.append((particle, monitor(particle, self.particles_plots[i])))

        if self.show_plots:
            self.particles_figure.show()

    def plot(self, data):
        """ Plot cosmological parameters and monitored particles distribution functions """

        last_t = data['t'].iloc[-1] / UNITS.s
        self.times.append(last_t)

        for i, plot in enumerate(self.plots):
            _, xmax = plot.get_xlim()
            ymin, ymax = plot.get_ylim()

            if last_t >= xmax:
                plot.set_xlim(self.times[0], last_t * 1.1)

            last_data = data[self.plot_map[i]].iloc[-1] / self.divider_map[i]
            self.plots_data[i].append(last_data)

            if last_data >= ymax:
                plot.set_ylim(self.plots_data[i].min, 1.1 * last_data)
            if last_data <= ymin:
                plot.set_ylim(last_data / 1.1, self.plots_data[i].max)

            self.lines[i].set_data(self.times, self.plots_data[i])

        self.params_figure.canvas.draw()

        if self.particles:
            for i, (_, monitor) in enumerate(self.particles):
                monitor.plot(data)

            self.particles_figure.canvas.draw()


class ParticleMonitor(object):
    def __init__(self, particle, plots):
        raise NotImplementedError()

    def plot(self, data):
        raise NotImplementedError()


class RadiationParticleMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].set_xlabel("T, MeV")
        self.plots[0].invert_xaxis()
        self.plots[0].set_xscale("log")
        self.plots[0].set_ylabel("rho/rho_eq")

        self.plots[1].set_xlabel("y, MeV")
        self.plots[1].set_xlim(self.particle.grid.MIN_MOMENTUM / UNITS.MeV,
                               self.particle.grid.MAX_MOMENTUM / UNITS.MeV)
        self.plots[1].set_ylabel("f/f_eq")

    def comparison_distributions(self, data):
        T = self.particle.params.T
        aT = self.particle.aT

        rhoeq = self.particle.energy_density() / (
            self.particle.dof * math.pi**2 / 30. * (aT / self.particle.params.a)**4
            * (7./8. if self.particle.statistics == STATISTICS.FERMION else 1.)
        )
        feq = self.particle.equilibrium_distribution(aT=aT)

        return (T, rhoeq), feq

    def plot(self, data):
        (T, rhoeq), feq = self.comparison_distributions(data)

        if not self.particle.in_equilibrium:
            ratio = numpy.vectorize(self.particle.distribution)(self.particle.grid.TEMPLATE) / feq
        else:
            ratio = numpy.ones(self.particle.grid.TEMPLATE.shape)
            rhoeq = 1.

        self.plots[0].scatter(T / UNITS.MeV, rhoeq, s=1)

        age_lines(self.plots[1].get_axes().lines)
        self.plots[1].plot(self.particle.grid.TEMPLATE / UNITS.MeV, ratio)


class EquilibriumRadiationParticleMonitor(RadiationParticleMonitor):
    def comparison_distributions(self, data):
        T = data['T'].iloc[-1]
        aT = data['aT'].iloc[-1]

        rhoeq = self.particle.energy_density() / (
            self.particle.dof * math.pi**2 / 30 * T**4
            * (7./8. if self.particle.statistics == STATISTICS.FERMION else 1.)
        )
        feq = self.particle.equilibrium_distribution(aT=aT)

        return (T, rhoeq), feq


class EffectiveTemperatureRadiationPartileMonitor(RadiationParticleMonitor):
    def comparison_distributions(self, data):
        rho = self.particle.energy_density()
        const = (
            self.particle.dof * math.pi**2 / 30.
            * (7./8. if self.particle.statistics == STATISTICS.FERMION else 1.)
        )

        T = self.particle.params.T
        T_eff = (rho / const)**0.25
        aT = T_eff * self.particle.params.a

        rhoeq = rho / const / (self.particle.aT / self.particle.params.a)**4
        feq = self.particle.equilibrium_distribution(aT=aT)

        self.plots[1].set_title("T ~ {:3e}".format(T_eff / UNITS.MeV))

        return (T, rhoeq), feq


class MassiveParticleMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].invert_xaxis()
        self.plots[0].set_xlabel("T, MeV")
        self.plots[0].set_xscale("log")
        self.plots[0].set_ylabel("rho/(n M)")

        self.plots[1].set_xlabel("y, MeV")
        self.plots[1].set_xlim(self.particle.grid.MIN_MOMENTUM / UNITS.MeV,
                               self.particle.grid.MAX_MOMENTUM / UNITS.MeV)
        self.plots[1].set_ylabel("(f-f_eq) y^2")

    def plot(self, data):
        T = data['T'].iloc[-1]

        from particles.NonEqParticle import energy_density, density

        self.plots[0].scatter(T / UNITS.MeV, energy_density(self.particle)
                              / (self.particle.mass * density(self.particle)), s=1)

        age_lines(self.plots[1].get_axes().lines)

        yy = self.particle.grid.TEMPLATE * self.particle.grid.TEMPLATE / UNITS.MeV**2

        f = self.particle._distribution
        feq = self.particle.equilibrium_distribution()
        self.plots[1].plot(self.particle.grid.TEMPLATE / UNITS.MeV, yy*(f-feq))


class EquilibrationMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].set_xlabel("a")
        self.plots[0].set_xscale("log")
        self.plots[0].set_ylabel("max|I|")

        self.plots[1].set_xlabel("a")
        self.plots[1].set_xscale("log")
        self.plots[1].set_ylabel("numerator, MeV^-1")

    def plot(self, data):
        a = data['a'].iloc[-1]

        from particles.NonEqParticle import numerator

        self.plots[0].scatter(a, numpy.max(numpy.fabs(self.particle.collision_integral))
                              * UNITS.MeV, s=1)
        self.plots[1].scatter(a, numerator(self.particle) * UNITS.MeV, s=1)


class AbundanceMonitor(ParticleMonitor):
    def __init__(self, particle, plots):
        self.particle, self.plots = particle, plots

        self.plots[0].set_title(particle.name)
        self.plots[0].set_xlabel("T, MeV")
        self.plots[0].invert_xaxis()
        self.plots[0].set_xscale("log")
        self.plots[0].set_yscale("log")
        self.plots[0].set_ylabel("rho fraction")

        self.plots[1].set_xlabel("T, MeV")
        self.plots[1].invert_xaxis()
        self.plots[1].set_xscale("log")
        self.plots[1].set_yscale("log")
        self.plots[1].set_ylabel("n a^3")

    def plot(self, data):
        T = data['T'].iloc[-1]

        total_rho = data['rho'].iloc[-1]
        rho = self.particle.energy_density()
        self.plots[0].scatter(T / UNITS.MeV, rho / total_rho, s=1)

        density = self.particle.density() * self.particle.params.a**3 / UNITS.MeV**3
        self.plots[1].scatter(T / UNITS.MeV, density, s=1)


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

    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D

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

    plt.savefig(os.path.join(os.path.split(filename)[0], 'logs/plt_{}.svg'.format(p0 / UNITS.MeV)))


def plot_points(points, name):
    """ Draw a scatter plot for a number of `points` tuples `(x, y)` """
    plt.figure(4)
    plt.title(name)
    plt.scatter(*zip(*points))
    plt.show()
