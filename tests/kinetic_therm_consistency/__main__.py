"""
## Kinetics-thermodynamics consistency test

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

This test checks that active neutrinos do not get any non-equilibrium corrections at temperatures\
$\sim 50 \div 10 MeV$

[Log file](log.txt)
[Distribution functions](distributions.txt)


"""

import os
import numpy
import matplotlib

from plotting import plt
from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params


folder = os.path.split(__file__)[0]

T_initial = 50. * UNITS.MeV
T_final = 10 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.025)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

universe.add_particles([
    photon,
    electron,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
])

neutrinos = [neutrino_e, neutrino_mu, neutrino_tau]
for neutrino in neutrinos:
    neutrino.decoupling_temperature = T_initial

universe.interactions += \
    SMI.neutrino_interactions(leptons=[electron], neutrinos=neutrinos)


universe.evolve(T_final)


""" ## Plots for comparison with articles """

"""
### JCAP10(2012)014, Figure 9
<img src="figure_9.svg" width=100% /> """

"""
### JCAP10(2012)014, Figure 10
<img src="figure_10.svg" width=100% />
<img src="figure_10_full.svg" width=100% /> """
