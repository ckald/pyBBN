"""
## Kinetics-thermodynamics consistency test

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

This test checks that active neutrinos do not get any non-equilibrium corrections at temperatures\
$\sim 1000 \div 100 MeV$

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params


folder = os.path.join(os.path.split(__file__)[0], 'output')

T_initial = 1000. * UNITS.MeV
T_final = 100 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.1)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)
tau = Particle(**SMP.leptons.tau)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

universe.add_particles([
    photon,
    electron,
    muon,
    tau,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
])

neutrinos = [neutrino_e, neutrino_mu, neutrino_tau]
for neutrino in neutrinos:
    neutrino.decoupling_temperature = T_initial

universe.interactions += \
    SMI.neutrino_interactions(leptons=[electron, muon, tau], neutrinos=neutrinos)


universe.evolve(T_final)


""" ## Plots for comparison with articles """

"""
### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% /> """

"""
### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% /> """
