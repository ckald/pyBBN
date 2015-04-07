"""
## Kinetics-thermodynamics consistency test

<img src="plots.png" width=100% />
<img src="particles.png" width=100% />

This test checks that active neutrinos do not get any non-equilibrium corrections at temperatures\
$\sim 50 \div 10 MeV$

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os
import numpy
import matplotlib

from plotting import plt
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params, GRID


params = Params(T_initial=50. * UNITS.MeV,
                T_final=10 * UNITS.MeV,
                dy=0.025)

universe = Universe(params=params,
                    logfile='tests/standard_model_bbn/log.txt')

photon = Particle(params=params, **SMP.photon)
electron = Particle(params=params, **SMP.leptons.electron)
neutrino_e = Particle(params=params, **SMP.leptons.neutrino_e)
neutrino_mu = Particle(params=params, **SMP.leptons.neutrino_mu)
neutrino_tau = Particle(params=params, **SMP.leptons.neutrino_tau)

universe.particles += [
    photon,
    electron,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
]

universe.interactions += \
    SMI.neutrino_interactions(leptons=[electron], neutrinos=[neutrino_e, neutrino_mu, neutrino_tau])

universe.graphics.monitor(particles=[
    neutrino_e,
    neutrino_mu,
    neutrino_tau
])


universe.evolve()

universe.graphics.save(__file__)


""" ## Plots for comparison with articles """

folder = os.path.split(__file__)[0]
plt.ion()

""" ### JCAP10(2012)014, Figure 9
    <img src="figure_9.png" width=100% /> """

plt.figure(9)
plt.title('Figure 9')
plt.xlabel('MeV/T')
plt.ylabel(u'aT')
plt.xscale('log')
plt.xlim(0.5, UNITS.MeV/universe.params.T_final)
plt.xticks([1, 2, 3, 5, 10, 20])
plt.axes().get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.plot(UNITS.MeV / numpy.array(universe.data['T']), numpy.array(universe.data['aT']) / UNITS.MeV)
plt.show()
plt.savefig(os.path.join(folder, 'figure_9.png'))

""" ### JCAP10(2012)014, Figure 10
    <img src="figure_10.png" width=100% />
    <img src="figure_10_full.png" width=100% /> """

plt.figure(10)
plt.title('Figure 10')
plt.xlabel('Conformal momentum y = pa')
plt.ylabel('f/f_eq')
plt.xlim(0, 20)

f_e = neutrino_e._distribution
feq_e = neutrino_e.equilibrium_distribution()
plt.plot(GRID.TEMPLATE/UNITS.MeV, f_e/feq_e, label="nu_e")

f_mu = neutrino_mu._distribution
feq_mu = neutrino_mu.equilibrium_distribution()
plt.plot(GRID.TEMPLATE/UNITS.MeV, f_mu/feq_mu, label="nu_mu")

f_tau = neutrino_tau._distribution
feq_tau = neutrino_tau.equilibrium_distribution()
plt.plot(GRID.TEMPLATE/UNITS.MeV, f_tau/feq_tau, label="nu_tau")

plt.legend()
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10_full.png'))

plt.xlim(0, 10)
plt.ylim(0.99, 1.06)
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10.png'))

# Distribution functions arrays
distributions_file = open(os.path.join(folder, 'distributions.txt'), "w")
numpy.savetxt(distributions_file, (f_e, feq_e, f_e/feq_e), header=str(neutrino_e),
              footer='-'*80, fmt="%1.5e")
numpy.savetxt(distributions_file, (f_mu, feq_mu, f_mu/feq_mu), header=str(neutrino_mu),
              footer='-'*80, fmt="%1.5e")
numpy.savetxt(distributions_file, (f_tau, feq_tau, f_tau/feq_tau), header=str(neutrino_tau),
              footer='-'*80, fmt="%1.5e")

distributions_file.close()

# Just to be sure everything is okay
import ipdb
ipdb.set_trace()

raw_input("...")
