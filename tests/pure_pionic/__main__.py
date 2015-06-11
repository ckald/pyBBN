# -*- coding: utf-8 -*-
"""

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

"""

import argparse
import os
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, utils, HeuristicGrid


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', required=True)
parser.add_argument('--theta', required=True)
parser.add_argument('--tau', required=True)
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)
lifetime = float(args.tau) * UNITS.s

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "mass={:e}_theta={:e}_tau={:e}".format(mass / UNITS.MeV, theta, lifetime / UNITS.s)
    + args.comment
)


params = Params(T_initial=200. * UNITS.MeV,
                T_final=0.0008 * UNITS.MeV,
                dy=0.025)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
tau = Particle(**SMP.leptons.tau)

neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

neutral_pion = Particle(**SMP.hadrons.neutral_pion)
charged_pion = Particle(**SMP.hadrons.charged_pion)

sterile = Particle(**NuP.dirac_sterile_neutrino(mass))


grid = HeuristicGrid(mass, lifetime)

sterile.decoupling_temperature = params.T_initial
for neutrino in [neutrino_e, neutrino_mu, neutrino_tau]:
    neutrino.decoupling_temperature = 10 * UNITS.MeV
    neutrino.set_grid(grid)


universe.add_particles([
    photon,

    electron,
    muon,
    tau,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    neutral_pion,
    charged_pion,

    sterile,
])

thetas = defaultdict(float, {
    'electron': theta,
})

universe.interactions += (
    SMI.neutrino_interactions(
        leptons=[electron, muon],
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau]
    )
    # + NuI.sterile_leptons_interactions(
    #     thetas=thetas, sterile=sterile,
    #     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    #     leptons=[electron, muon, tau]
    # )
    + NuI.sterile_hadrons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron, muon, tau],
        hadrons=[neutral_pion, charged_pion]
    )
)

universe.init_kawano(electron=electron, neutrino=neutrino_e)
universe.init_oscillations(SMP.leptons.oscillations_map(), (neutrino_e, neutrino_mu, neutrino_tau))

if universe.graphics:
    from plotting import (RadiationParticleMonitor, EffectiveTemperatureRadiationPartileMonitor,
                          MassiveParticleMonitor, AbundanceMonitor)
    universe.graphics.monitor([
        (neutrino_e, EffectiveTemperatureRadiationPartileMonitor),
        (neutrino_mu, RadiationParticleMonitor),
        (neutrino_tau, RadiationParticleMonitor),
        (sterile, MassiveParticleMonitor),
        (sterile, AbundanceMonitor)
    ])

universe.evolve()

"""
### Plots for comparison with articles

### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% />

### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% />
"""

if universe.graphics:
    from tests.plots import articles_comparison_plots
    articles_comparison_plots(universe, [neutrino_e, neutrino_mu, neutrino_tau, sterile])
