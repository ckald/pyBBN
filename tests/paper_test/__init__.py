# -*- coding: utf-8 -*-
"""
## Heavy sterile dirac neutrino

$$ M = 33.9 MeV $$

$$ \theta_\tau \approx 7.6 10^{-3} \sim \tau_N \approx 0.3 sec $$

http://arxiv.org/pdf/hep-ph/0002223v2.pdf

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

"""

import argparse
import os
import numpy
import matplotlib
from collections import defaultdict

from plotting import plt, RadiationParticleMonitor, MassiveParticleMonitor, EquilibrationMonitor
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, GRID, utils


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', required=True)
parser.add_argument('--theta', required=True)
args = parser.parse_args()

mass = float(args.mass)
theta = float(args.theta)

folder = utils.ensure_dir(os.path.split(__file__)[0],
                          "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta))


params = Params(T_initial=100. * UNITS.MeV,
                T_final=0.0008 * UNITS.MeV,
                dy=0.025)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)
sterile = Particle(**NuP.sterile_neutrino(mass))

sterile.decoupling_temperature = params.T_initial
neutrino_e.decoupling_temperature = 10 * UNITS.MeV
neutrino_mu.decoupling_temperature = 10 * UNITS.MeV
neutrino_tau.decoupling_temperature = 10 * UNITS.MeV

universe.add_particles([
    photon,
    electron,
    muon,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
    sterile,
])

thetas = defaultdict(float, {
    'electron': theta,
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

universe.init_kawano(electron=electron, neutrino=neutrino_e)

universe.graphics.monitor([
    (neutrino_e, RadiationParticleMonitor),
    (neutrino_mu, RadiationParticleMonitor),
    (neutrino_tau, RadiationParticleMonitor),
    (sterile, MassiveParticleMonitor),
    (sterile, EquilibrationMonitor)
])

universe.evolve()


""" ## Plots for comparison with articles"""

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
plt.xlim(0.5, UNITS.MeV/universe.params.T_final)
plt.xticks([1, 2, 3, 5, 10, 20])
plt.axes().get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.plot(UNITS.MeV / numpy.array(universe.data['T']), numpy.array(universe.data['aT']) / UNITS.MeV)
plt.show()
plt.savefig(os.path.join(folder, 'figure_9.svg'))

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
distributions_file = open(os.path.join(folder, 'distributions.txt'), "w")

for neutrino in [neutrino_e, neutrino_mu, neutrino_tau, sterile]:
    f = neutrino._distribution
    feq = neutrino.equilibrium_distribution()
    plt.plot(GRID.TEMPLATE/UNITS.MeV, f/feq, label=neutrino.name)

    numpy.savetxt(distributions_file, (f, feq, f/feq), header=str(neutrino),
                  footer='-'*80, fmt="%1.5e")

plt.legend()
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10_full.svg'))

plt.xlim(0, 10)
plt.ylim(0.99, 1.06)
plt.draw()
plt.show()
plt.savefig(os.path.join(folder, 'figure_10.svg'))

distributions_file.close()
