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
from common import UNITS, Params, utils, LogSpacedGrid
from scipy.integrate import simps, cumtrapz
import numpy as np

mass = 150 * UNITS.MeV
theta = 1e-3

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)
)

T_initial = 5. * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125)

universe = Universe(params=params, folder=folder)

from common import LinearSpacedGrid
linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=101, MAX_MOMENTUM=100*UNITS.MeV)
linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron, **{'grid': linear_grid})
muon = Particle(**SMP.leptons.muon, **{'grid': linear_grid, 'thermal': False})#, 'decoupling_temperature': T_initial})

neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': linear_grid})
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, **{'grid': linear_grid})
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, **{'grid': linear_grid})

for neutrino in [neutrino_e, neutrino_mu, neutrino_tau]:
    neutrino.decoupling_temperature = 0 * UNITS.MeV


sterile = Particle(**NuP.dirac_sterile_neutrino(mass), **{'grid': linear_grid_s})
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
    'electron': theta#,
    #'muon': theta,
    #'tau': theta
})

interactions_Ff = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    kind=CollisionIntegralKind.F_f_vacuum_decay
)

interactions_F1 = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    kind=CollisionIntegralKind.F_1_vacuum_decay
)

def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)

def Count(reaction):
    val = Counter(item.specie.name for item in reaction if
            item.specie.name in ['Muon', 'Sterile neutrino (Dirac)'])
    return val['Muon'] * val['Sterile neutrino (Dirac)']

for inter in interactions_Ff:
    inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
                        and reaction_type(integral.reaction) in [2]]

for inter in interactions_F1:
    inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
                        and reaction_type(integral.reaction) in [-2]]

universe.interactions += (interactions_Ff + interactions_F1)

def step_monitor(universe):
    if len(universe.data) > 1:
        for particle in universe.particles:
            data = particle.data['params']
            if particle.mass > 0:
                if particle.name == 'Sterile neutrino (Dirac)':
                    com_den_HNL_1st = data['density'][-1] * particle.params.a**3
                    com_den_HNL_2nd = data['density'][-2] * particle.params.a**3
                    delta_den_HNL = com_den_HNL_2nd - com_den_HNL_1st
                if particle.name == 'Muon':
                    com_den_muon_1st = data['density'][-1] * particle.params.a**3
                    com_den_muon_2nd = data['density'][-2] * particle.params.a**3
                    delta_den_muon = com_den_muon_2nd - com_den_muon_1st
        # print(delta_den_HNL, delta_den_muon)
        with open(os.path.join(folder, particle.name.replace(' ', '_') + "Densities.txt"), 'a') as f1:
            f1.write('{T:e}\t{a:e}\t{nHNL:e}\t{nmu:e}\n'
                .format(T=particle.params.T / UNITS.MeV,
                        a=particle.params.a,
                        nHNL=com_den_HNL_1st / UNITS.MeV**3,
                        nmu=com_den_muon_1st / UNITS.MeV**3))

universe.step_monitor = step_monitor

#for aa in interactions_F1:
#   print(aa)

universe.evolve(T_final)

