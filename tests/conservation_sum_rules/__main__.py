"""
Compute and save collision integrals of active neutrinos
in the case of elastic processes only.
"""

import os
import sys
from collections import Counter

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import UNITS, Params
from common import LinearSpacedGrid


folder = os.path.join(os.path.split(__file__)[0], 'output')

T_final = 0.1 * UNITS.MeV
params = Params(T=10 * UNITS.MeV, dy=0.003125)

universe = Universe(params=params, folder=folder)

linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=201, MAX_MOMENTUM=20 * UNITS.MeV)

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)
neutrino_e = Particle(**SMP.leptons.neutrino_e, grid=linear_grid)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, grid=linear_grid)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, grid=linear_grid)

neutrino_e.decoupling_temperature = 10 * UNITS.MeV
neutrino_mu.decoupling_temperature = 10 * UNITS.MeV
neutrino_tau.decoupling_temperature = 10 * UNITS.MeV

universe.add_particles([photon, electron, neutrino_e, neutrino_mu, neutrino_tau])

interactions_pre = (
    SMI.neutrino_interactions(leptons=[electron], neutrinos=[neutrino_e, neutrino_mu, neutrino_tau])
)

integrals = [[] for n in range(len(interactions_pre))]

for num, interaction in enumerate(interactions_pre):
    for integral in interaction.integrals:
        left_species = Counter(item.specie for item in integral.reaction if item.side == -1)
        right_species = Counter(item.specie for item in integral.reaction if item.side == 1)

        if left_species == right_species:
            integrals[num].append(integral)

for num, interaction in enumerate(interactions_pre):
    interaction.integrals = integrals[num]


universe.interactions += interactions_pre

def step_monitor(universe):

    # Output the collision integrals to file every step, first row is momentum grid
    # Output the number and energy densities of active neutrinos
    if universe.step % 1 == 0:
        for particle in [neutrino_e, neutrino_mu, neutrino_tau]:
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".collision_integrals.txt"), 'a') as f1:
                if universe.step == 1:
                    f1.write('\t'.join(['{:e}'.format(x) for x in particle.grid.TEMPLATE / UNITS.MeV]) + '\n')
                f1.write('\t'.join(['{:e}'.format(x) for x in particle.collision_integral]) + '\n')
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".number_and_energy_density.txt"), 'a') as f2:
                f2.write('{:e}'.format(particle.density / UNITS.MeV**3) + '\t' + '{:e}'.format(universe.total_energy_density() / UNITS.MeV**4) + '\n')

universe.step_monitor = step_monitor

universe.evolve(T_final)
