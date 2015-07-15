"""
## Standard Model BBN test

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

This test checks that in the universe filled with photons, electrons and neutrinos:

  * $a * T$ is not conserved by a factor around `1.401` and precise details of this process
  * neutrino non-equilibrium corrections reproduce the results of the Dolgov-Hansen-Semikoz papers

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os

from plotting import RadiationParticleMonitor
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params


folder = os.path.split(__file__)[0]

T_kawano = 12 * UNITS.MeV
T_simple = 0.05 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_kawano,
                dy=0.0125)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu_tau = Particle(**SMP.leptons.neutrino_mu)
neutrino_mu_tau.dof = 4

neutron = Particle(**SMP.hadrons.neutron)
proton = Particle(**SMP.hadrons.proton)

universe.add_particles([
    photon,
    electron,
    neutrino_e,
    neutrino_mu_tau,

    neutron,
    proton
])

universe.init_oscillations(SMP.leptons.oscillations_map(), (neutrino_e, neutrino_mu_tau))

universe.interactions += (
    SMI.neutrino_interactions(leptons=[electron],
                              neutrinos=[neutrino_e, neutrino_mu_tau])
)

universe.init_kawano(electron=electron, neutrino=neutrino_e)

if universe.graphics:
    universe.graphics.monitor([
        (neutrino_e, RadiationParticleMonitor),
        (neutrino_mu_tau, RadiationParticleMonitor)
    ])


universe.evolve(T_simple, export=False)
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

from tests.plots import articles_comparison_plots
articles_comparison_plots(universe, [neutrino_e, neutrino_mu_tau])
