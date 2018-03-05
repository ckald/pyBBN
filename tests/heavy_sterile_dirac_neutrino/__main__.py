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
parser.add_argument('--tau', default='0.3')
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
lifetime = float(args.tau) * UNITS.s
# theta = 0.5 * numpy.arcsin(numpy.sqrt(4*5.7*10**(-4)/float(args.tau)))
theta = 0.031
T_dec = 50 * UNITS.MeV
print('theta=', theta, ' Tdec=', T_dec / UNITS.MeV)


folder = os.path.join(os.path.split(__file__)[0], "output", args.tau)

T_initial = T_dec
T_washout = 0.1 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125 * 4)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
muon = Particle(**SMP.leptons.muon)

from common import LogSpacedGrid
active_grid = LogSpacedGrid(MOMENTUM_SAMPLES=201, MAX_MOMENTUM=1.5 * mass / (T_washout / UNITS.MeV))
sterile_grid = LogSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=mass)

neutrino_e = Particle(**SMP.leptons.neutrino_e, grid=active_grid)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, grid=active_grid)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, grid=active_grid)
sterile = Particle(**NuP.dirac_sterile_neutrino(mass), grid=sterile_grid)

sterile.decoupling_temperature = T_initial

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


def step_monitor(universe):
    # explanation of what is inside the file
    if universe.step == 1:
        with open(os.path.join(folder, "sterile_densities.txt"), 'a') as f:
            f.write('## a     rho_sterile     n_sterile' + '\n')
        with open(os.path.join(folder, "sterile_distribution.txt"), 'a') as f:
            f.write('# First line is a grid of y; Starting from second line: first number is a, second is temperature, next is set of numbers is corresponding to f(y) on the grid' + '\n')
            f.write('## a     T     ' + '\t'.join([
                '{:e}'.format(x)
                for x in
                sterile.grid.TEMPLATE / UNITS.MeV
            ]) + '\n')
    # Output the density and energy density of sterile neutrino
    if universe.step % 10 == 0:
        with open(os.path.join(folder,"sterile_densities.txt"), 'a') as f:
            f.write('{:e}'.format(universe.params.a) + '\t' + '{:e}'.format(universe.params.T/UNITS.MeV) + '\t'+'{:e}'.format(universe.params.aT/UNITS.MeV) + '\t'+'{:e}'.format(sterile.energy_density/(UNITS.MeV)**4) + '\t')
            f.write('{:e}'.format(sterile.density/(UNITS.MeV)**3) + '\n')
        with open(os.path.join(folder, "sterile_distribution.txt"), 'a') as f:
            f.write('{:e}'.format(universe.params.a) + '\t'+'{:e}'.format(universe.params.T/UNITS.MeV) + '\t')
            f.write('\t'.join([
                '{:e}'.format(x)
                for x in sterile._distribution
            ]) + '\n')


universe.step_monitor = step_monitor

universe.evolve(5 * UNITS.MeV, export=False)
universe.params.dy = 0.003125
universe.params.infer()

universe.evolve(T_washout, export=False)
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

