# -*- coding: utf-8 -*-

import argparse
import os
from collections import defaultdict, Counter, namedtuple
from interactions.four_particle import FourParticleIntegral

os.environ['SPLIT_COLLISION_INTEGRAL'] = ''

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI, SterileM
from interactions.four_particle.cpp.integral import CollisionIntegralKind
from evolution import Universe
from common import UNITS, Params, utils, LinearSpacedGrid
from scipy.integrate import simps
import numpy as np

mass = 150 * UNITS.MeV
theta = 1e-3

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)
)

T_initial = 5 * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125)

universe = Universe(params=params, folder=folder)

linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=301, MAX_MOMENTUM=600*UNITS.MeV)
linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron, **{'grid': linear_grid})
muon = Particle(**SMP.leptons.muon, **{'grid': linear_grid, 'thermal_dyn': False})

neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': linear_grid, 'thermal_dyn': False})
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, **{'grid': linear_grid, 'thermal_dyn': False})
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, **{'grid': linear_grid})

for neutrino in [neutrino_e, neutrino_tau]:
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
    'electron': theta,
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


interactions_F1 = utils.interaction_filter(
                ['Electron'],
                interactions_F1
                )

interactions_Ff = utils.interaction_filter(
                ['Electron'],
                interactions_Ff
                )

def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)

## Test for N <--> electron + muon + muon neutrino
# def Count(reaction):
#     val = Counter(item.specie.name for item in reaction if
#             item.specie.name in ['Muon', 'Sterile neutrino (Dirac)', 'Electron'])
#     return val['Muon'] * val['Sterile neutrino (Dirac)'] * val['Electron']

# for inter in interactions_Ff:
#     inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
#                         and reaction_type(integral.reaction) in [2]]

# for inter in interactions_F1:
#     inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
#                         and reaction_type(integral.reaction) in [-2]
#                         and integral.reaction[0].specie.name in ["Muon", "Muon neutrino"]]


## Test for N <--> 2 muon neutrinos + 1 electron neutrino
# def Count(reaction):
#     val = Counter(item.specie.name for item in reaction if
#             item.specie.name in ['Muon neutrino', 'Sterile neutrino (Dirac)', 'Electron neutrino'])
#     return val['Muon neutrino'] * val['Sterile neutrino (Dirac)'] * val['Electron neutrino']

# for inter in interactions_Ff:
#     inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
#                         and reaction_type(integral.reaction) in [2]]

# for inter in interactions_F1:
#     inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
#                         and reaction_type(integral.reaction) in [-2]
#                         and integral.reaction[0].specie.name in ["Muon neutrino", "Electron neutrino"]]

## Test for N <--> 3 electron neutrinos
def Count(reaction):
    val = Counter(item.specie.name for item in reaction if
            item.specie.name in ['Sterile neutrino (Dirac)', 'Electron neutrino'])
    return val['Sterile neutrino (Dirac)'] * val['Electron neutrino']

for inter in interactions_Ff:
    inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
                        and reaction_type(integral.reaction) in [2]
                        and integral.reaction[2].specie.name == 'Electron neutrino']

for inter in interactions_F1:
    inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
                        and reaction_type(integral.reaction) in [-2]
                        and integral.reaction[0].specie.name in ["Electron neutrino"]
                        and integral.reaction[1].specie.name == 'Electron neutrino']


# Add integral manually
# IntegralItem = namedtuple('IntegralItem', ['side', 'specie', 'antiparticle', 'index', 'crossed'])
# reaction = (IntegralItem(side=-1, specie=neutrino_e, antiparticle=False, index=0, crossed=False),
#             IntegralItem(side=-1, specie=neutrino_e, antiparticle=False, index=1, crossed=False),
#             IntegralItem(side=-1, specie=neutrino_e, antiparticle=True, index=2, crossed=False),
#             IntegralItem(side=1, specie=sterile, antiparticle=False, index=3, crossed=False))

# Ms = (SterileM(theta=theta, K1=4., order=(0, 1, 2, 3)),)

# interactions_F1[0].integrals.append(FourParticleIntegral(
#             particle=neutrino_e,
#             reaction=reaction,
#             washout_temperature=0,
#             Ms=Ms,
#             kind=CollisionIntegralKind.F_1_vacuum_decay
#         ))



universe.interactions += (interactions_Ff + interactions_F1)


def step_monitor(universe):
    if len(universe.data) > 1:
        for particle in universe.particles:
            data = particle.data['params']
            if particle.name == 'Sterile neutrino (Dirac)':
                com_den_HNL_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_HNL_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_HNL = com_den_HNL_1st - com_den_HNL_2nd
            # if particle.name == 'Muon neutrino':
            #     com_den_muon_1st = data['density'][-1] * universe.data['a'][-1]**3
            #     com_den_muon_2nd = data['density'][-2] * universe.data['a'][-2]**3
            #     delta_den_muon = com_den_muon_1st - com_den_muon_2nd
            if particle.name == 'Electron neutrino':
                com_den_nu_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_nu_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_nu = com_den_nu_1st - com_den_nu_2nd
        print(delta_den_HNL, delta_den_nu)
        with open(os.path.join(folder, particle.name.replace(' ', '_') + "Densities.txt"), 'a') as f1:
            f1.write('{T:e}\t{a:e}\t{nHNL:e}\t{nmu:e}\t{nnu:e}\n'
                .format(T=particle.params.T / UNITS.MeV,
                        a=particle.params.a,
                        nHNL=com_den_HNL_1st / UNITS.MeV**3,
                        # nmu=com_den_muon_1st / UNITS.MeV**3,
                        nnu=com_den_nu_1st / UNITS.MeV**3))


universe.step_monitor = step_monitor

universe.evolve(T_final)