import os
import csv
from itertools import izip

import numpy
import matplotlib

from plotting import plt
from common import UNITS


def articles_comparison_plots(universe, particles):

    """ Ready for copy-paste docstring:

        ### Plots for comparison with articles

        ### JCAP10(2012)014, Figure 9
        <img src="figure_9.svg" width=100% />

        ### JCAP10(2012)014, Figure 10
        <img src="figure_10.svg" width=100% />
        <img src="figure_10_full.svg" width=100% />
    """

    plt.ion()

    """
    ### JCAP10(2012)014, Figure 9
    <img src="figure_9.svg" width=100% />
    """

    plt.figure(9)
    plt.title('Figure 9')
    plt.xlabel('MeV/T')
    plt.ylabel(u'aT')
    plt.xscale('log')
    plt.xlim(0.5, UNITS.MeV/universe.params.T)
    plt.xticks([1, 2, 3, 5, 10, 20])
    plt.axes().get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.plot(UNITS.MeV / numpy.array(universe.data['T']),
             numpy.array(universe.data['aT']) / UNITS.MeV)
    plt.show()
    plt.savefig(os.path.join(universe.folder, 'figure_9.svg'))

    """
    ### JCAP10(2012)014, Figure 10
    <img src="figure_10.svg" width=100% />
    <img src="figure_10_full.svg" width=100% />
    """

    plt.figure(10)
    plt.title('Figure 10')
    plt.xlabel('Conformal momentum y = pa')
    plt.ylabel('f/f_eq')
    plt.xlim(0, 20)

    # Distribution functions arrays
    distributions_file = open(os.path.join(universe.folder, 'distributions.txt'), "w")

    for neutrino in particles:
        f = neutrino._distribution
        feq = neutrino.equilibrium_distribution()
        plt.plot(neutrino.grid.TEMPLATE/UNITS.MeV, f/feq, label=neutrino.name)

        numpy.savetxt(distributions_file, (f, feq, f/feq), header=str(neutrino),
                      footer='-'*80, fmt="%1.5e")

    plt.legend()
    plt.draw()
    plt.show()
    plt.savefig(os.path.join(universe.folder, 'figure_10_full.svg'))

    plt.xlim(0, 10)
    plt.ylim(0.99, 1.06)
    plt.draw()
    plt.show()
    plt.savefig(os.path.join(universe.folder, 'figure_10.svg'))

    distributions_file.close()


def cosmic_neutrino_temperature(universe):
    data = universe.data
    with open(os.path.join(universe.folder, 'sm_temp.dat'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['# x, MeV', 'aT, MeV'])
        for i, row in data.iterrows():
            writer.writerow([row['x'] / UNITS.MeV, row['aT'] / UNITS.MeV])


def spectrum(universe, particle):
    with open(os.path.join(universe.folder, 'spectrum_{}.dat'.format(particle.name)), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['# y, MeV', 'f'])
        for y, f in izip(particle.grid.TEMPLATE, particle._distribution):
            writer.writerow([y / UNITS.MeV, f])


def spectrum_distortion(universe, particle):
    with open(os.path.join(universe.folder, 'distortion_{}.dat'.format(particle.name)), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['# y, MeV', 'y^2 (f-f_eq)'])
        for y, f in izip(particle.grid.TEMPLATE, particle._distribution):
            writer.writerow([
                y / UNITS.MeV,
                y**2 * (f - particle.equilibrium_distribution(y))
                / UNITS.MeV**2
            ])


def energy_density_deviation(universe, particle):
    with open(os.path.join(universe.folder,
              'energy_deviation_{}.dat'.format(particle.name)), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['# a', 'rho/rho_eq - 1'])
        for a, rho in izip(universe.data['a'], particle.data['energy_density']):
            writer.writerow([
                a,
                rho / (7 * particle.dof * numpy.pi**2 / 240 * (universe.params.T)**4) - 1
            ])
