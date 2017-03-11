"""
## Sterile neutrinos decoupling above $\Lambda_{QCD}$

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

This tests simulates the decoupling of sterile neutrinos in the quark-gluon plasma.

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os
import numpy
import matplotlib

from collections import defaultdict

from plotting import plt, RadiationParticleMonitor, MassiveParticleMonitor
from particles import Particle
from library.SM import particles as SMP  # , interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, GRID

folder = os.path.split(__file__)[0]

T_initial = 1.5 * UNITS.GeV
T_final = 100 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.025)

universe = Universe(params=params, logfile=os.path.join(folder, 'log.txt'))

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
tau = Particle(**SMP.leptons.tau)

neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

up = Particle(**SMP.quarks.up)
down = Particle(**SMP.quarks.down)
# charm = Particle(**SMP.quarks.charm)
strange = Particle(**SMP.quarks.strange)
# top = Particle(**SMP.quarks.top)
# bottom = Particle(**SMP.quarks.bottom)

sterile = Particle(**NuP.sterile_neutrino(300 * UNITS.MeV))
sterile.decoupling_temperature = T_initial

universe.add_particles([
    photon,

    electron,
    muon,
    tau,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    up,
    down,
    strange,

    sterile,
])

thetas = defaultdict(float, {
    'electron': 1e-4,
})

universe.interactions += (
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
    )
)

universe.evolve(T_final)

""" ## Plots for comparison with articles """

"""
### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% /> """

"""
### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% /> """
