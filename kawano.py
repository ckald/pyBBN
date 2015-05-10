# -*- coding: utf-8 -*-

import itertools
import numpy
import pandas
from collections import namedtuple
from common import GRID, CONST
from common.integrators import integrate_1D
from library import SM


q = (SM.particles.hadrons.neutron['mass'] - SM.particles.hadrons.proton['mass'])
m_e = SM.particles.leptons.electron['mass']

a = None

Particles = namedtuple("Particles", "electron neutrino")
particles = None


default_bounds = (GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,)


def init_kawano(electron=None, neutrino=None):
    global particles
    particles = Particles(electron=electron, neutrino=neutrino)


def _rate1(y):
    """ n + ν_e ⟶  e + p """
    E_e = q*a + y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * (1. - particles.electron.distribution(y_e)) * particles.neutrino.distribution(y))


def _rate1b(y):
    """ e + p ⟶  n + ν_e """
    E_e = q*a + y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * particles.electron.distribution(y_e) * (1. - particles.neutrino.distribution(y)))


def _rate2(y):
    """ n ⟶  e + ν_e' + p """
    E_e = q*a - y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * (1. - particles.electron.distribution(y_e))
            * (1. - particles.neutrino.distribution(y)))


def _rate2b(y):
    """ e + ν_e' + p ⟶  n """
    E_e = q*a - y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y**2 * y_e * E_e
            * particles.electron.distribution(y_e) * particles.neutrino.distribution(y))


def _rate3(y):
    """ n + e' ⟶  ν_e' + p """
    E_e = -q*a + y
    y_e = numpy.sqrt(E_e**2 - (m_e*a)**2)
    return (y_e**2 * y_e * E_e
            * particles.electron.distribution(y_e) * (1. - particles.neutrino.distribution(y)))


def _rate3b(y):
    """ ν_e' + p ⟶  n + e' """
    E_e = numpy.sqrt(y**2 + (m_e*a)**2)
    y_n = q*a + E_e

    return (y**2 * y_n**2
            * (1. - particles.electron.distribution(y)) * particles.neutrino.distribution(y_n))


def baryonic_rates(_a):
    global a
    a = _a

    return [
        CONST.rate_normalization / a**5 * integrate_1D(numpy.vectorize(integrand), bounds=bounds)[0]
        if bounds[0] < bounds[1] else 0.
        for integrand, bounds in [
            (_rate1, default_bounds),
            (_rate1b, default_bounds),
            (_rate2, (GRID.MIN_MOMENTUM, q*a - (m_e*a))),
            (_rate2b, (GRID.MIN_MOMENTUM, q*a - (m_e*a))),
            (_rate3, (q*a + (m_e*a), GRID.MAX_MOMENTUM)),
            (_rate3b, default_bounds)]
    ]


Plotting = namedtuple('Plotting', 'figure plots')
parameters_plots = None
rates_plots = None

heading = ("t[s]", "x",
           "Tg[10^9K]", "dTg/dt[10^9K/s]",
           "rho_tot[g cm^-3]", "H[s^-1]",
           "n nue->p e", "p e->n nue",
           "n->p e nue", "p e nue->n",
           "n e->p nue", "p nue->n e")


def import_data(filepath):
    with open(filepath) as f:
        line = f.readline()
        try:
            line = line.split()
            if int(line[0]):
                f.readline()
        except:
            pass

        data = pandas.DataFrame(
            ({heading[i]: float(value) for i, value in enumerate(line.strip("\n").split("\t"))}
             for line in f), columns=heading)
    return data


def plot_kawano(data, label=None):
    global parameters_plots, rates_plots

    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    plt.rcParams['toolbar'] = 'None'

    if not parameters_plots:
        figure, plots = plt.subplots(3, 2, num="KAWANO parameters")
        plots = list(itertools.chain(*plots))
        figure.subplots_adjust(hspace=0.7, wspace=0.5)

        #     t[s]         x    Tg[10^9K]   dTg/dt[10^9K/s] rho_tot[g cm^-3]     H[s^-1]

        parameters_plots = Plotting(figure=figure, plots=plots)

        for plot in plots:
            plot.set_xlabel("time, s")
            plot.set_xscale("log")
            plot.set_yscale("log")

        plots[0].set_title("Scale factor")
        plots[0].set_ylabel("a, 1")
        plots[0].set_yscale("log")

        plots[1].set_title("Temperature")
        plots[1].set_ylabel("T, 10^9 K")

        plots[2].set_title("Temperature derivative")
        plots[2].set_ylabel("dT/dt, 10^9 K/s")
        plots[2].set_yscale('linear')

        plots[3].set_title("Total energy density")
        plots[3].set_ylabel("rho, g/cm^3")

        plots[4].set_title("Hubble rate")
        plots[4].set_ylabel("H, s^-1")

        plots[5].set_title("Nuclear rates")
        plots[5].set_ylabel("rate")

    if not rates_plots:
        figure, plots = plt.subplots(3, 2, num="KAWANO rates")
        plots = list(itertools.chain(*plots))
        figure.subplots_adjust(hspace=0.7, wspace=0.5)

        # n nue->p e  p e->n nue  n->p e nue  p e nue->n  n e->p nue  p nue->n e

        rates_plots = Plotting(figure=figure, plots=plots)

        for i, plot in enumerate(plots, 6):
            plot.set_title(heading[i])
            plot.set_xlabel("time, s")
            plot.set_xscale("log")
            plot.set_ylabel("Rate")
            plot.set_yscale("log")

    time_series = data[heading[0]]
    parameters_plots.plots[0].plot(time_series, data[heading[1]])
    parameters_plots.plots[1].plot(time_series, data[heading[2]])
    parameters_plots.plots[2].plot(time_series, data[heading[3]])
    parameters_plots.plots[3].plot(time_series, data[heading[4]])
    parameters_plots.plots[4].plot(time_series, data[heading[5]])

    rates = data.ix[:, 6:12]
    for i, rate in enumerate(rates):
        parameters_plots.plots[5].plot(time_series, rates[rate])
        rates_plots.plots[i].plot(time_series, rates[rate], label=label)
