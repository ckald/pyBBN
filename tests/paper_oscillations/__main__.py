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
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, utils, HeuristicGrid
from semianalytic import step_monitor


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', required=True)
parser.add_argument('--theta', required=True)
parser.add_argument('--tau', required=True)
parser.add_argument('--Tdec', default=100)
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)
lifetime = float(args.tau) * UNITS.s
T_dec = float(args.Tdec) * UNITS.MeV

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "mass={:e}_tau={:e}_theta={:e}".format(mass / UNITS.MeV, lifetime / UNITS.s, theta)
    + args.comment
)


T_kawano = 12 * UNITS.MeV
T_initial = max(T_dec, T_kawano)
T_interactions_freeze_out = 0.05 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.05)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)
sterile = Particle(**NuP.dirac_sterile_neutrino(mass))

grid = HeuristicGrid(mass, lifetime)

sterile.decoupling_temperature = T_dec
for neutrino in [neutrino_e, neutrino_mu, neutrino_tau]:
    neutrino.decoupling_temperature = 5 * UNITS.MeV
    neutrino.set_grid(grid)

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
universe.init_oscillations(SMP.leptons.oscillations_map(), (neutrino_e, neutrino_mu, neutrino_tau))
universe.step_monitor = step_monitor


if universe.graphics:
    from plotting import RadiationParticleMonitor, MassiveParticleMonitor, AbundanceMonitor
    universe.graphics.monitor([
        (neutrino_e, RadiationParticleMonitor),
        (neutrino_mu, RadiationParticleMonitor),
        (neutrino_tau, RadiationParticleMonitor),
        (sterile, MassiveParticleMonitor),
        (sterile, AbundanceMonitor)
    ])

universe.evolve(T_kawano, export=False)
universe.params.dy = 0.0125
universe.evolve(T_interactions_freeze_out, export=False)
universe.interactions = tuple()
universe.evolve(T_final)

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
