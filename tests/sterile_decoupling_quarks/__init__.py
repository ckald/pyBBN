"""
## Sterile neutrinos decoupling above $\Lambda_{QCD}$

<img src="plots.png" width=100% />
<img src="particles.png" width=100% />

This tests simulates the decoupling of sterile neutrinos in the quark-gluon plasma.

[Log file](log.txt)
[Distribution functions](distributions.txt)


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


params = Params(T_initial=1. * UNITS.GeV,
                T_final=100 * UNITS.MeV)

universe = Universe(params=params,
                    logfile='tests/sterile_decoupling_quarks/log.txt')

photon = Particle(params=params, **SMP.photon)

electron = Particle(params=params, **SMP.leptons.electron)
muon = Particle(params=params, **SMP.leptons.muon)
tau = Particle(params=params, **SMP.leptons.tau)

neutrino_e = Particle(params=params, **SMP.leptons.neutrino_e)
neutrino_mu = Particle(params=params, **SMP.leptons.neutrino_mu)
neutrino_tau = Particle(params=params, **SMP.leptons.neutrino_tau)

up = Particle(params=params, **SMP.quarks.up)
down = Particle(params=params, **SMP.quarks.down)
# charm = Particle(params=params, **SMP.quarks.charm)
strange = Particle(params=params, **SMP.quarks.strange)
# top = Particle(params=params, **SMP.quarks.top)
# bottom = Particle(params=params, **SMP.quarks.bottom)

sterile = Particle(params=params,
                   **NuP.sterile_neutrino(400 * UNITS.MeV))

universe.particles += [
    photon,

    electron,
    muon,
    tau,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    up,
    down,
    # charm,
    strange,
    # top,
    # bottom,

    sterile
]

thetas = defaultdict(float, {
    'electron': 1e-4,
})

universe.interactions += (
    SMI.neutrino_interactions(
        leptons=[electron, muon, tau],
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau]
    ) +
    NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron, muon, tau]
    )
    + NuI.sterile_quark_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron, muon, tau],
        quarks=[up, down, strange]
        # quarks=[up, down, charm, strange, top, bottom]
    )
)

universe.graphics.monitor(particles=[
    neutrino_e,
    sterile
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
