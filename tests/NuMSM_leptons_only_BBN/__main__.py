mass = 20
theta = 1e-1

import numpy as np
import argparse
import os
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, utils, LogSpacedGrid

mass *= UNITS.MeV

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}_flav={}".format(mass / UNITS.MeV, theta, "all")
)

T_initial = 50 * UNITS.MeV
T_freeze = 0.01 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.0005)

universe = Universe(params=params, folder=folder)

from common import LinearSpacedGrid
linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=101, MAX_MOMENTUM=100*UNITS.MeV)
linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron, grid=linear_grid)

neutrino_e = Particle(**SMP.leptons.neutrino_e, grid=linear_grid)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, grid=linear_grid)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, grid=linear_grid)

neutrino_e.decoupling_temperature = 5. * UNITS.MeV
neutrino_mu.decoupling_temperature = 5. * UNITS.MeV
neutrino_tau.decoupling_temperature = 5. * UNITS.MeV

sterile = Particle(**NuP.dirac_sterile_neutrino(mass), grid=linear_grid_s)
sterile.decoupling_temperature = T_initial

universe.add_particles([
    photon,

    electron,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    sterile,
])

thetas = defaultdict(float, {
    'electron': theta,
    'muon': theta,
    'tau': theta
})

interaction = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron]
)

universe.interactions += (
    NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron]
    ) +
    SMI.neutrino_interactions(
        leptons=[electron],
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    )
)


universe.init_kawano(electron=electron, neutrino=neutrino_e)

def step_monitor(universe):
    if universe.step == 1:
        for particle in [neutrino_e, neutrino_mu, neutrino_tau]:
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".distribution.txt"), 'a') as f:
                f.write('# First line is a grid of y; Starting from second line: first number is a, second is temperature, next is set of numbers is corresponding to f(y) on the grid' + '\n')
                f.write('## a     T     ' + '\t'.join([
                    '{:e}'.format(x)
                    for x in particle.grid.TEMPLATE / UNITS.MeV
                ]) + '\n')
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".collision_integrals.txt"), 'a') as f:
                f.write('# First line is a grid of y; Starting from second line each line is a set of numbers is corresponding to Icoll(f) on the grid y with temperature equal to the T in .distribution.txt' + '\n')
                f.write('##     ' + '\t'.join([
                    '{:e}'.format(x)
                    for x in particle.grid.TEMPLATE / UNITS.MeV
                ]) + '\n')
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".rho.txt"), 'a') as f:
                f.write('## a     T     aT    rho_nu' + '\n')
        with open(os.path.join(folder, "rho_nu.txt"), 'a') as f:
                    f.write('# a    rho_nu' + '\n')

    if universe.step % 10 == 0:
        for particle in [neutrino_e, neutrino_mu, neutrino_tau]:
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".distribution.txt"), 'a') as f:
                f.write(
                    '{:e}\t{:e}\t'.format(universe.params.a, universe.params.T / UNITS.MeV)
                    + '\t'.join([
                        '{:e}'.format(x)
                        for x in particle._distribution
                    ]) + '\n'
                )
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".rho.txt"), 'a') as f:
                f.write('{:e}\t{:e}\t{:e}\t{:e}\n'.format(
                    universe.params.a, universe.params.T / UNITS.MeV,
                    universe.params.aT / UNITS.MeV, particle.energy_density / UNITS.MeV**4
                ))
        with open(os.path.join(folder, "rho_nu.txt"), 'a') as f:
            f.write('{:e}\t{:e}\n'.format(
                universe.params.a,
                sum(p.energy_density for p in [neutrino_e, neutrino_mu, neutrino_tau]) / UNITS.MeV**4
            ))

universe.step_monitor = step_monitor

universe.evolve(T_final)

