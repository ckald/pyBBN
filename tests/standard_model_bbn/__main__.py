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

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params


folder = os.path.join(os.path.split(__file__)[0], 'output')

T_interaction_freezeout = 0.05 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=10. * UNITS.MeV,
                dy=0.003125)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_mu.dof = 4

universe.add_particles([
    photon,
    electron,
    neutrino_e,
    neutrino_mu
])

universe.interactions += (
    SMI.neutrino_interactions(leptons=[electron],
                              neutrinos=[neutrino_e, neutrino_mu])
)

universe.init_kawano(electron=electron, neutrino=neutrino_e)


def step_monitor(universe):
    # explanation of what is inside the file + first row which is a grid on y
    if universe.step == 1:
        for particle in [neutrino_e, neutrino_mu]:
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".distribution.txt"), 'a') as f:
                f.write('# First line is a grid of y; Starting from second line: first number is temperature, next set of numbers is corresponding f/feq on the grid' + '\n')
                f.write('## T     ' + '\t'.join([
                    '{:e}'.format(x)
                    for x in
                    particle.grid.TEMPLATE / UNITS.MeV
                ]) + '\n')

    # Output the distribution function distortion to file every 10 steps, first column is temperature
    if universe.step % 10 == 0:
        for particle in [neutrino_e, neutrino_mu]:
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".distribution.txt"), 'a') as f:
                f.write ('{:e}'.format(universe.params.T/UNITS.MeV) + '\t')
                f.write('\t'.join([
                    '{:e}'.format(x)
                    for x in
                    (particle._distribution
                     / particle.equilibrium_distribution(aT=universe.params.aT/UNITS.MeV))
                ]) + '\n')


universe.step_monitor = step_monitor

universe.evolve(T_interaction_freezeout, export=False)
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
