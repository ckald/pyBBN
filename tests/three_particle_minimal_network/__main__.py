# -*- coding: utf-8 -*-
"""
## Minimal network with two-body decays

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

Study a network of Boltzmann equation to get control over instabilities related to the addition of
decays. For simplicity investigate a minimal network which exhibits the and has the following
features:

  * ms > mπ so the decay channel into pions opens
  * neglect back-reaction → standard cosmological evolution
  * minimal set of interactions for equilibration of sterile neutrinos and all other species
  * no mixing

Take into account sterile neutrinos (single generation) νs, pions π0 (may be treated as photons)
and a single generation of active neutrinos να (possibly also treated as photons)

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os

from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params


folder = os.path.join(os.path.split(__file__)[0], 'output')

T_initial = 200 * UNITS.MeV
T_final = 50 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.025)

universe = Universe(params=params, logfile=os.path.join(folder, 'log.txt'))

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron)

neutrino_e = Particle(**SMP.leptons.neutrino_e)

neutral_pion = Particle(**SMP.hadrons.neutral_pion)

sterile = Particle(**NuP.sterile_neutrino(300 * UNITS.MeV))
sterile.decoupling_temperature = T_initial

universe.add_particles([
    photon,

    electron,

    neutrino_e,

    neutral_pion,

    sterile,
])

thetas = defaultdict(float, {
    'electron': 1e-3,
})

universe.interactions += (
    SMI.neutrino_interactions(
        leptons=[electron],
        neutrinos=[neutrino_e]
    ) +
    NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[electron]
    )
    + NuI.sterile_hadrons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[electron],
        hadrons=[neutral_pion]
    )
)

universe.evolve(T_final)


""" ## Plots for comparison with articles """

""" ### JCAP10(2012)014, Figure 9
    <img src="figure_9.svg" width=100% /> """

""" ### JCAP10(2012)014, Figure 10
    <img src="figure_10.svg" width=100% />
    <img src="figure_10_full.svg" width=100% /> """
