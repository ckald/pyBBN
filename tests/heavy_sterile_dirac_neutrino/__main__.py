# -*- coding: utf-8 -*-
"""
## Heavy sterile dirac neutrino

$$ M = 33.9 MeV $$

$$ \theta_\tau \approx 4.86 10^{-2} \sim \tau_N \approx 0.3 sec $$

http://arxiv.org/pdf/hep-ph/0002223v2.pdf

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

"""

import os
import argparse
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', default='33.9')
parser.add_argument('--theta', default='0.0486')
parser.add_argument('--tau', default='0.3')
parser.add_argument('--Tdec', default='5')
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)
lifetime = float(args.tau) * UNITS.s
T_dec = float(args.Tdec) * UNITS.MeV


folder = os.path.join(os.path.split(__file__)[0], "output", args.tau)

T_initial = max(50. * UNITS.MeV, T_dec)
T_interaction_freezeout = 0.05 * UNITS.MeV
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

sterile.decoupling_temperature = T_initial
neutrino_e.decoupling_temperature = 5 * UNITS.MeV
neutrino_mu.decoupling_temperature = 5 * UNITS.MeV
neutrino_tau.decoupling_temperature = 5 * UNITS.MeV

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
    'tau': theta,
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

universe.evolve(T_interaction_freezeout, export=False)
universe.interactions = tuple()
universe.params.dy = 0.0125
universe.evolve(T_final)


"""
### Plots for comparison with articles

### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% />

### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% />
"""
