"""
## Massive $\nu_\tau$ ($20 MeV$) test

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

This test checks that in the universe filled with photons, electrons and neutrinos:

  * $a * T$ is not conserved by a factor around `1.477` and precise details of this process
  * neutrino non-equilibrium corrections reproduce the results of the Dolgov-Hansen-Semikoz papers

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os
import argparse

from plotting import plt, RadiationParticleMonitor, MassiveParticleMonitor
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.MassiveNuTau import particles as MNP, interactions as MNI
from evolution import Universe
from common import UNITS, Params, GRID


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', default='20')
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV


folder = os.path.join(os.path.split(__file__)[0], 'mass={}'.format(args.mass))

T_initial = 20. * UNITS.MeV
T_final = 0.015 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.1)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**MNP.leptons.neutrino_tau)
neutrino_tau.mass = mass

neutrino_e.decoupling_temperature = T_initial
neutrino_mu.decoupling_temperature = T_initial
neutrino_tau.decoupling_temperature = T_initial


universe.add_particles([
    photon,
    electron,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
])

universe.interactions += (
    SMI.neutrino_interactions(leptons=[electron], neutrinos=[neutrino_e, neutrino_mu])
    + MNI.neutrino_scattering(neutrino_e, neutrino_tau)
    + MNI.neutrino_scattering(neutrino_mu, neutrino_tau)
)

universe.evolve(T_final)


"""
### Plots for comparison with articles

### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% />

### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% />
"""
