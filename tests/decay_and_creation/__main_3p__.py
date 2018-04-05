# -*- coding: utf-8 -*-

import argparse
import os
from collections import defaultdict, Counter

os.environ['SPLIT_COLLISION_INTEGRAL'] = ''

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from interactions.four_particle.cpp.integral import CollisionIntegralKind
from evolution import Universe
from common import UNITS, Params, utils, LinearSpacedGrid, LogSpacedGrid
from scipy.integrate import simps
import numpy as np

mass = 150 * UNITS.MeV
theta = 1e-3

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}_3p".format(mass / UNITS.MeV, theta)
)

T_initial = 20. * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125)

universe = Universe(params=params, folder=folder)

from common import LinearSpacedGrid
linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=1001, MAX_MOMENTUM=500*UNITS.MeV)
log_grid = LogSpacedGrid(MOMENTUM_SAMPLES=1001, MAX_MOMENTUM=200*UNITS.MeV)
linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)

photon = Particle(**SMP.photon)
neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': log_grid, 'thermal_dyn': False})
neutral_pion = Particle(**SMP.hadrons.neutral_pion, **{'grid': linear_grid, 'thermal_dyn': False})

sterile = Particle(**NuP.dirac_sterile_neutrino(mass), **{'grid': linear_grid_s})
sterile.decoupling_temperature = T_initial

universe.add_particles([
    photon,

    neutrino_e,

    neutral_pion,

    sterile,
])


thetas = defaultdict(float, {
    'electron': theta,
})

interactions_Ff = NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],
        leptons = [],
        mesons = [neutral_pion],
        kind=CollisionIntegralKind.F_f_vacuum_decay
)

interactions_F1 = NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],
        leptons = [],
        mesons = [neutral_pion],
        kind=CollisionIntegralKind.F_1_vacuum_decay
)

def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)


def Count(reaction):
    val = Counter(item.specie.name for item in reaction if
            item.specie.name in ['Muon', 'Sterile neutrino (Dirac)', 'Electron'])
    return val['Muon'] * val['Sterile neutrino (Dirac)'] * val['Electron']

for inter in interactions_Ff:
    inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [1]]

for inter in interactions_F1:
    inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [-1]
                        and integral.reaction[0].specie.name in ["Neutral pion", "Electron neutrino"]]


universe.interactions += (interactions_Ff + interactions_F1)


def step_monitor(universe):
    if len(universe.data) > 1:
        for particle in universe.particles:
            data = particle.data['params']
            if particle.name == 'Sterile neutrino (Dirac)':
                com_den_HNL_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_HNL_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_HNL = com_den_HNL_1st - com_den_HNL_2nd
            if particle.name == 'Neutral pion':
                com_den_pion_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_pion_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_pion = com_den_pion_1st - com_den_pion_2nd
            if particle.name == 'Electron neutrino':
                com_den_nu_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_nu_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_nu = com_den_nu_1st - com_den_nu_2nd
        print(delta_den_HNL, delta_den_pion * 2, delta_den_nu * 2)
        with open(os.path.join(folder, particle.name.replace(' ', '_') + "Densities_3p.txt"), 'a') as f1:
            f1.write('{T:e}\t{a:e}\t{nHNL:e}\t{nmu:e}\t{nnu:e}\n'
                .format(T=particle.params.T / UNITS.MeV,
                        a=particle.params.a,
                        nHNL=com_den_HNL_1st / UNITS.MeV**3,
                        nmu=com_den_pion_1st / UNITS.MeV**3,
                        nnu=com_den_nu_1st / UNITS.MeV**3))

universe.step_monitor = step_monitor

universe.evolve(T_final)
