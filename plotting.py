import os
import itertools

import numpy
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
# from multiprocessing import Process, Pipe

from collections import deque

from common import UNITS, GRID, PARAMS


class ring_deque(deque):
    max = 0
    min = 0

    def __init__(self, data, length):
        self.length = length
        super(ring_deque, self).__init__(data)

    def append_more(self, data):

        if len(self) > self.length:
            self.popleft()

        self.max = max(data, self.max)
        self.min = min(data, self.min)
        super(ring_deque, self).append(data)

    def append(self, data):
        self.max = data
        self.min = data

        self.append = self.append_more

        super(ring_deque, self).append(data)


class Plotting:
    particles = None

    def __init__(self):

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
        self.plots[0].set_ylim(0, 10)

        self.plots[1].set_title("Scale factor")
        self.plots[1].set_xlabel("time, s")
        self.plots[1].set_ylabel("a, 1")
        self.plots[0].set_ylim(1, 10)

        self.plots[2].set_title("T * a")
        self.plots[2].set_xlabel("time, s")
        self.plots[2].set_ylabel("T * a, MeV")
        self.plots[0].set_ylim(1, 10)

        self.plots[3].set_title("Total energy density")
        self.plots[3].set_yscale("log")
        self.plots[3].set_xlabel("time, s")
        self.plots[3].set_ylabel("rho, eV**4")
        self.plots[0].set_ylim(1, 10)

        self.lines = []
        self.plots_data = []
        self.times = ring_deque([], 1000)
        for plot in self.plots:
            self.lines.append(plot.plot([], [], 'b-')[0])
            self.plots_data.append(ring_deque([], 1000))

        self.params_figure.show()

        # self.pipe, self.plot_pipe = Pipe(False)
        # self.plot_process = Process(target=self,
        #                             name="Plotting process")
        # self.plot_process.daemon = True
        # self.plot_process.start()

    def save(self, filename):
        folder = os.path.split(filename)[0]
        plt.figure(1)
        plt.savefig(os.path.join(folder, 'plots.png'))
        if self.particles:
            plt.figure(2)
            plt.savefig(os.path.join(folder, 'particles.png'))

    def monitor(self, particles=[]):
        self.particles = particles
        if self.particles:
            self.particles_figure, self.particles_plots = plt.subplots(len(particles)*2, 1, num=2)
            self.particles_figure.subplots_adjust(hspace=0.5, wspace=0.5)

            for i, particle in enumerate(self.particles):
                self.particles_plots[i*2].set_title(particle.name)
                self.particles_plots[i*2].set_xlabel("y, MeV")
                self.particles_plots[i*2].set_ylabel("f")

                self.particles_plots[i*2 + 1].set_title(particle.name)
                self.particles_plots[i*2 + 1].set_xlabel("y, MeV")
                self.particles_plots[i*2 + 1].set_ylabel("f/f_eq")

            self.particles_figure.show()

    def plot(self, data, full=False):

        last_t = data['t'][-1] / UNITS.s
        self.times.append(last_t)

        for i, plot in enumerate(self.plots):
            xmin, xmax = plot.get_xlim()
            ymin, ymax = plot.get_ylim()

            if last_t >= xmax:
                plot.set_xlim(self.times[0], last_t * 1.5)

            last_data = data[self.plot_map[i]][-1] / self.divider_map[i]
            self.plots_data[i].append(last_data)

            if last_data >= ymax:
                plot.set_ylim(self.plots_data[i].min, 1.5 * last_data)
            if last_data <= ymin:
                plot.set_ylim(last_data / 1.5, self.plots_data[i].max)

            self.lines[i].set_data(self.times, self.plots_data[i])

        plt.figure(1)
        plt.draw()

        if self.particles:
            for i, particle in enumerate(self.particles):
                if not particle.in_equilibrium:
                    self.particles_plots[i*2].plot(GRID.TEMPLATE / UNITS.MeV,
                                                   particle._distribution)

                    effective_temperature = (
                        240. * particle.energy_density() / 7. / numpy.pi**2 / particle.dof
                    ) ** 0.25

                    feq = particle.distribution_function(
                        numpy.vectorize(particle.energy_normalized)(GRID.TEMPLATE)
                        / PARAMS.a / effective_temperature
                    )
                    self.particles_plots[i*2 + 1].set_title(
                        particle.name + " (T' = {})".format(effective_temperature / UNITS.MeV)
                    )
                    for line in self.particles_plots[i*2 + 1].get_axes().lines:
                        alpha = line.get_alpha() or 1.
                        if alpha < 0.1:
                            line.remove()
                        else:
                            line.set_alpha((line.get_alpha() or 1.) * 0.85)

                    self.particles_plots[i*2 + 1].plot(
                        GRID.TEMPLATE / UNITS.MeV,
                        particle._distribution / feq
                    )

            plt.figure(2)
            plt.draw()

    # def terminate(self):
    #     plt.close('all')

    # def poll_draw(self):
    #     print 'poll_draw'
    #     while 1:
    #         print 'Callback'
    #         if not self.pipe.poll(None):
    #             print 'not'
    #             break

    #         command = self.pipe.recv()
    #         print command
    #         if command is None:
    #             self.terminate()
    #             return False

    #         else:
    #             self.plot(command[0], full=command[1])

    #     return True

    # def __call__(self):
    #     print 'Called'
    #     # self.gid = gobject.timeout_add(100, self.poll_draw())
    #     self.poll_draw()


def plot_integrand(integrand, particle, p0):
    fig = plt.figure(3)
    ax = fig.gca(projection='3d')
    plt.cla()
    X, Y = numpy.meshgrid(GRID.TEMPLATE, GRID.TEMPLATE)
    Z = numpy.array([integrand([x, y]) for x, y in zip(numpy.ravel(X), numpy.ravel(Y))])\
        .reshape(X.shape)

    ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=0.1)
    ax.contourf(X, Y, Z, zdir='z', offset=numpy.amin(Z), cmap=cm.coolwarm)
    ax.contourf(X, Y, Z, zdir='x', offset=ax.get_xlim()[0], cmap=cm.coolwarm)
    ax.contourf(X, Y, Z, zdir='y', offset=ax.get_ylim()[1], cmap=cm.coolwarm)

    ax.set_xlabel('p1')
    ax.set_ylabel('p2')
    ax.set_title('{} p0 = {}'.format(particle.name, p0 / UNITS.MeV))

    plt.savefig(os.path.join(os.path.split(__file__)[0], 'logs/plt_{}.png'.format(p0 / UNITS.MeV)))
