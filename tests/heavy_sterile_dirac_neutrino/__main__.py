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
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params


folder = os.path.split(__file__)[0]

T_initial = 100. * UNITS.MeV
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
sterile = Particle(**NuP.dirac_sterile_neutrino(33.9 * UNITS.MeV))

sterile.decoupling_temperature = T_initial
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
    'tau': 4.86e-2,
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

if universe.graphics:
    from plotting import RadiationParticleMonitor, DensityAndEnergyMonitor
    universe.graphics.monitor([
        (neutrino_e, RadiationParticleMonitor),
        (neutrino_mu, RadiationParticleMonitor),
        (neutrino_tau, RadiationParticleMonitor),
        (sterile, DensityAndEnergyMonitor)
    ])

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

    import os
    import csv
    from itertools import izip
    density_data, energy_data = universe.graphics.particles[3].data

    with open(os.path.join(universe.folder, 'normalized_density_plot.dat'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for x, y in izip(*density_data):
            writer.write([x, y])

    with open(os.path.join(universe.folder, 'normalized_energy_density_plot.dat'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for x, y in izip(*energy_data):
            writer.write([x, y])
