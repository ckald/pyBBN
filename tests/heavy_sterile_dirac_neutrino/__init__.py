# -*- coding: utf-8 -*-
"""
## Heavy sterile dirac neutrino

$$ M = 33.9 MeV $$

$$ \theta_\tau \approx 7.6 10^{-4} \sim \tau_N \approx 0.3 sec $$

http://arxiv.org/pdf/hep-ph/0002223v2.pdf

"""

import os
import numpy
import matplotlib
from collections import defaultdict

from plotting import plt
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, GRID


folder = os.path.split(__file__)[0]

params = Params(T_initial=100. * UNITS.MeV,
                T_final=1. * UNITS.MeV,
                dy=0.025)

universe = Universe(params=params, logfile=os.path.join(folder, 'log.txt'))

photon = Particle(params, **SMP.photon)
electron = Particle(params, **SMP.leptons.electron)
muon = Particle(params, **SMP.leptons.muon)
neutrino_e = Particle(params, **SMP.leptons.neutrino_e)
neutrino_mu = Particle(params, **SMP.leptons.neutrino_mu)
neutrino_tau = Particle(params, **SMP.leptons.neutrino_tau)
sterile = Particle(params, **NuP.sterile_neutrino(33.9 * UNITS.MeV))

sterile.decoupling_temperature = params.T_initial

universe.particles += [
    photon,
    electron,
    muon,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
    sterile,
]

thetas = defaultdict(float, {
    'tau': 7.6 * 1e-4,
})

universe.interactions += (
    SMI.neutrino_interactions(
        leptons=[electron],
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau]
    ) + NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron, muon]
    )
)

universe.graphics.monitor(particles=[
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
    sterile
])


universe.evolve()

universe.graphics.save(__file__)


""" ## Plots for comparison with articles"""

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

f_sterile = sterile._distribution
feq_sterile = sterile.equilibrium_distribution()
plt.plot(GRID.TEMPLATE/UNITS.MeV, f_sterile/feq_sterile, label="sterile")

plt.legend()
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10_full.png'))

plt.xlim(0, 10)
plt.ylim(0.99, 1.06)
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10.png'))


""" ### JCAP10(2012)014, Figure 12
    <img src="figure_12.png" width=100% /> """

plt.figure(12)
plt.title("Figure 12")
plt.xlabel("Scale factor a, 1")
plt.ylabel("rho/(n M)")
plt.xlim(params.a_initial, 1)
plt.ylim(0.99, 4.01)

import itertools
from particles.NonEqParticle import energy_density, density

for distribution, a in itertools.izip(sterile.data['distribution'], universe.data['a']):
    sterile._distribution = distribution
    plt.scatter(a, energy_density(sterile) / (sterile.mass * density(sterile)), s=1)

plt.savefig(os.path.join(folder, 'figure_12.png'))

# Distribution functions arrays
distributions_file = open(os.path.join(folder, 'distributions.txt'), "w")
numpy.savetxt(distributions_file, (f_e, feq_e, f_e/feq_e), header=str(neutrino_e),
              footer='-'*80, fmt="%1.5e")
numpy.savetxt(distributions_file, (f_mu, feq_mu, f_mu/feq_mu), header=str(neutrino_mu),
              footer='-'*80, fmt="%1.5e")
numpy.savetxt(distributions_file, (f_tau, feq_tau, f_tau/feq_tau), header=str(neutrino_tau),
              footer='-'*80, fmt="%1.5e")
numpy.savetxt(distributions_file, (f_sterile, feq_sterile, f_sterile/feq_sterile),
              header=str(sterile), footer='-'*80, fmt="%1.5e")

distributions_file.close()

# Just to be sure everything is okay
import ipdb
ipdb.set_trace()

raw_input("...")
