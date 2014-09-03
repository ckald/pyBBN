import numpy
import matplotlib.pyplot as plt
import numericalunits as nu

from common import UNITS, GRID


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
        self.plots[1][1].set_ylabel("rho, MeV**4")

        self.line00, = self.plots[0][0].plot([1], [1], 'b.')
        self.line01, = self.plots[0][1].plot([1], [1], 'b.')
        self.line10, = self.plots[1][0].plot([1], [1], 'b.')
        self.line11, = self.plots[1][1].plot([1], [1], 'b.')

        self.params_figure.show()

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

        self.line00.set_xdata(ts / UNITS.s)
        self.line00.set_ydata(Ts / nu.MeV)
        self.line01.set_xdata(ts / UNITS.s)
        self.line01.set_ydata(As)
        self.line10.set_xdata(ts / UNITS.s)
        self.line10.set_ydata(Ts * As / nu.MeV)
        self.line11.set_xdata(ts / UNITS.s)
        self.line11.set_ydata(rhos / nu.MeV**4)
        # line00, = self.plots[0][0].plot(ts / UNITS.s, Ts / nu.MeV, style, animated=True)
        # line01, = self.plots[0][1].plot(ts / UNITS.s, As, style, animated=True)
        # line10, = self.plots[1][0].plot(ts / UNITS.s, Ts * As / nu.MeV, style, animated=True)
        # line11, = self.plots[1][1].plot(ts / UNITS.s, rhos / nu.MeV**4, style, animated=True)

        # print(line00, line01, line10, line11)

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
