import os
import itertools
import matplotlib.pyplot as plt
from collections import deque

from common import UNITS, GRID


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
                    self.particles_plots[i*2].plot(GRID.TEMPLATE / UNITS.MeV, particle._distribution)
                    feq = particle.distribution_function_vectorized(
                        particle.energy_normalized_vectorized(GRID.TEMPLATE)
                        / particle.T
                    )
                    self.particles_plots[i*2 + 1].plot(
                        GRID.TEMPLATE / UNITS.MeV,
                        particle._distribution / feq
                    )

            plt.figure(2)
            plt.draw()
