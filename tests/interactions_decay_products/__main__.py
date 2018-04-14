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
import numpy as np

mass = 150 * UNITS.MeV
theta = 5e-3

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)
)

T_initial = 5. * UNITS.MeV
T_final = 0.01 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.0001)

universe = Universe(params=params, folder=folder)

from common import LinearSpacedGrid
linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=201, MAX_MOMENTUM=200*UNITS.MeV)
linear_grid_nu = LinearSpacedGrid(MOMENTUM_SAMPLES=401, MAX_MOMENTUM=300*UNITS.MeV)
linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)

photon = Particle(**SMP.photon)

electron = Particle(**SMP.leptons.electron, **{'grid': linear_grid})
muon = Particle(**SMP.leptons.muon, **{'grid': linear_grid, 'thermal_dyn': False})

charged_pion = Particle(**SMP.hadrons.charged_pion, **{'grid': linear_grid, 'thermal_dyn': False})
neutral_pion = Particle(**SMP.hadrons.neutral_pion, **{'grid': linear_grid, 'thermal_dyn': False})

neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': linear_grid_nu, 'thermal_dyn': False})
neutrino_mu = Particle(**SMP.leptons.neutrino_mu, **{'grid': linear_grid_nu, 'thermal_dyn': False})
neutrino_tau = Particle(**SMP.leptons.neutrino_tau, **{'grid': linear_grid_nu})

neutrino_tau.decoupling_temperature = 0 * UNITS.MeV


sterile = Particle(**NuP.dirac_sterile_neutrino(mass), **{'grid': linear_grid_s})
sterile.decoupling_temperature = T_initial

universe.add_particles([
    photon,

    electron,
    muon,

    charged_pion,
    neutral_pion,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    sterile,
])

thetas = defaultdict(float, {
    'electron': theta
})

interactions_SM = SMI.neutrino_interactions(
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron]
)


## N --> mu + e + nu_mu and mu --> e + nu_mu + nu_e
interactions_primary_Ff = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    kind=CollisionIntegralKind.F_f_vacuum_decay
)

interactions_primary_F1 = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    kind=CollisionIntegralKind.F_1_vacuum_decay
)

interactions_secondary_Ff = NuI.interactions_decay_products(
    interactions_primary=[interactions_primary_Ff],
    interactions_SM=interactions_SM,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[electron, muon],
    kind=CollisionIntegralKind.F_f_vacuum_decay
)

interactions_secondary_F1 = NuI.interactions_decay_products(
    interactions_primary=[interactions_primary_F1],
    interactions_SM=interactions_SM,
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

for inter in interactions_primary_Ff:
    inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
                        and reaction_type(integral.reaction) in [2]]

for inter in interactions_primary_F1:
    inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
                        and reaction_type(integral.reaction) in [-2]]

for inter in interactions_secondary_Ff:
    inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [2]]

for inter in interactions_secondary_F1:
    inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [-2]]


## N --> pi + e and pi --> mu + nu_mu and mu --> e + nu_e + nu_mu
# interactions_primary_Ff = NuI.sterile_hadrons_interactions(
#     thetas=thetas, sterile=sterile,
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     leptons=[electron, muon],
#     mesons=[charged_pion],
#     kind=CollisionIntegralKind.F_f_vacuum_decay
# )

# interactions_primary_F1 = NuI.sterile_hadrons_interactions(
#     thetas=thetas, sterile=sterile,
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     leptons=[electron, muon],
#     mesons=[charged_pion],
#     kind=CollisionIntegralKind.F_1_vacuum_decay
# )

# interactions_secondary_Ff = NuI.interactions_decay_products(
#     interactions_primary=[interactions_primary_Ff],
#     interactions_SM=interactions_SM,
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     leptons=[electron, muon],
#     mesons=[charged_pion],
#     kind=CollisionIntegralKind.F_f_vacuum_decay
# )

# interactions_secondary_F1 = NuI.interactions_decay_products(
#     interactions_primary=[interactions_primary_F1],
#     interactions_SM=interactions_SM,
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     leptons=[electron, muon],
#     mesons=[charged_pion],
#     kind=CollisionIntegralKind.F_1_vacuum_decay
# )

# def reaction_type(reaction):
#     return sum(reactant.side for reactant in reaction)

# def Count(reaction):
#     val = Counter(item.specie.name for item in reaction if
#             item.specie.name in ['Electron', 'Sterile neutrino (Dirac)'])
#     return val['Electron'] * val['Sterile neutrino (Dirac)']

# for inter in interactions_primary_Ff:
#     inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
#                         and reaction_type(integral.reaction) in [1]]

# for inter in interactions_primary_F1:
#     inter.integrals = [integral for integral in inter.integrals if Count(integral.reaction)
#                         and reaction_type(integral.reaction) in [-1]]

# for inter in interactions_secondary_Ff:
#     inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [1, 2]]

# for inter in interactions_secondary_F1:
#     inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [-1, -2]]


## N --> pi + nu_e and pi --> gamma + gamma
# interactions_primary_Ff = NuI.sterile_hadrons_interactions(
#     thetas=thetas, sterile=sterile,
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     leptons=[],
#     mesons=[neutral_pion],
#     kind=CollisionIntegralKind.F_f_vacuum_decay
# )

# interactions_primary_F1 = NuI.sterile_hadrons_interactions(
#     thetas=thetas, sterile=sterile,
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     leptons=[],
#     mesons=[neutral_pion],
#     kind=CollisionIntegralKind.F_1_vacuum_decay
# )

# interactions_secondary_Ff = NuI.interactions_decay_products(
#     interactions_primary=[interactions_primary_Ff],
#     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
#     photon=[photon],
#     mesons=[neutral_pion],
#     kind=CollisionIntegralKind.F_f_vacuum_decay
# )

# def reaction_type(reaction):
#     return sum(reactant.side for reactant in reaction)

# for inter in interactions_primary_Ff:
#     inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [1]]

# for inter in interactions_primary_F1:
#     inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [-1]]

# for inter in interactions_secondary_Ff:
#     inter.integrals = [integral for integral in inter.integrals if reaction_type(integral.reaction) in [1]]


interactions_primary_F1 = utils.interaction_filter(
                ['Electron'],
                interactions_primary_F1
                )

interactions_primary_Ff = utils.interaction_filter(
                ['Electron'],
                interactions_primary_Ff
                )

interactions_secondary_F1 = utils.interaction_filter(
                ['Electron'],
                interactions_secondary_F1
                )

interactions_secondary_Ff = utils.interaction_filter(
                ['Electron'],
                interactions_secondary_Ff
                )

universe.interactions += (interactions_primary_Ff + interactions_primary_F1 + interactions_secondary_Ff + interactions_secondary_F1)


## N --> mu + e + nu_mu and mu --> e + nu_mu + nu_e
def step_monitor(universe):
    if len(universe.data) > 1:
        for particle in universe.particles:
            data = particle.data['params']
            if particle.name == 'Sterile neutrino (Dirac)':
                com_den_HNL_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_HNL_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_HNL = com_den_HNL_1st - com_den_HNL_2nd
            if particle.name == 'Muon':
                com_den_muon_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_muon_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_muon = com_den_muon_1st - com_den_muon_2nd
            if particle.name == 'Muon neutrino':
                com_den_numu_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_numu_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_numu = com_den_numu_1st - com_den_numu_2nd
            if particle.name == 'Electron neutrino':
                com_den_nuel_1st = data['density'][-1] * universe.data['a'][-1]**3
                com_den_nuel_2nd = data['density'][-2] * universe.data['a'][-2]**3
                delta_den_nuel = com_den_nuel_1st - com_den_nuel_2nd
      #  print(delta_den_HNL, delta_den_muon, delta_den_numu, delta_den_nuel)
        with open(os.path.join(folder, particle.name.replace(' ', '_') + "Densities_mu_dec.txt"), 'a') as f1:
            f1.write('{T:e}\t{a:e}\t{nHNL:e}\t{nmu:e}\t{nnumu:e}\t{nnuel:e}\n'
                .format(T=particle.params.T / UNITS.MeV,
                        a=particle.params.a,
                        nHNL=com_den_HNL_1st / UNITS.MeV**3,
                        nmu=com_den_muon_1st / UNITS.MeV**3,
                        nnumu=com_den_numu_1st / UNITS.MeV**3,
                        nnuel=com_den_nuel_1st / UNITS.MeV**3))


## N --> pi + e and pi --> mu + nu_mu and mu --> e + nu_e + nu_mu
# def step_monitor(universe):
#     if len(universe.data) > 1:
#         for particle in universe.particles:
#             data = particle.data['params']
#             if particle.name == 'Sterile neutrino (Dirac)':
#                 com_den_HNL_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_HNL_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_HNL = com_den_HNL_1st - com_den_HNL_2nd
#             if particle.name == 'Charged pion':
#                 com_den_pion_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_pion_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_pion = com_den_pion_1st - com_den_pion_2nd
#             if particle.name == 'Muon':
#                 com_den_muon_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_muon_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_muon = com_den_muon_1st - com_den_muon_2nd
#             if particle.name == 'Muon neutrino':
#                 com_den_numu_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_numu_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_numu = com_den_numu_1st - com_den_numu_2nd
#             if particle.name == 'Electron neutrino':
#                 com_den_nuel_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_nuel_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_nuel = com_den_nuel_1st - com_den_nuel_2nd
#         print("\n\n", delta_den_HNL, delta_den_pion, delta_den_muon, delta_den_numu, delta_den_nuel, "\n\n")
#         with open(os.path.join(folder, particle.name.replace(' ', '_') + "Densities_pi_mu_dec.txt"), 'a') as f1:
#             f1.write('{T:e}\t{a:e}\t{nHNL:e}\t{npion:e}\t{nmu:e}\t{nnumu:e}\t{nnuel:e}\n'
#                 .format(T=particle.params.T / UNITS.MeV,
#                         a=particle.params.a,
#                         nHNL=com_den_HNL_1st / UNITS.MeV**3,
#                         npion=com_den_pion_1st / UNITS.MeV**3,
#                         nmu=com_den_muon_1st / UNITS.MeV**3,
#                         nnumu=com_den_numu_1st / UNITS.MeV**3,
#                         nnuel=com_den_nuel_1st / UNITS.MeV**3))


## N --> pi + nu_e and pi --> gamma + gamma
# def step_monitor(universe):
#     if len(universe.data) > 1:
#         for particle in universe.particles:
#             data = particle.data['params']
#             if particle.name == 'Sterile neutrino (Dirac)':
#                 com_den_HNL_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_HNL_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_HNL = com_den_HNL_1st - com_den_HNL_2nd
#             if particle.name == 'Neutral pion':
#                 com_den_pion_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_pion_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_pion = com_den_pion_1st - com_den_pion_2nd
#             if particle.name == 'Electron neutrino':
#                 com_den_nuel_1st = data['density'][-1] * universe.data['a'][-1]**3
#                 com_den_nuel_2nd = data['density'][-2] * universe.data['a'][-2]**3
#                 delta_den_nuel = com_den_nuel_1st - com_den_nuel_2nd
#         # print("\n\n", delta_den_HNL, delta_den_pion, delta_den_nuel, universe.data['T'][-1], "\n\n")
#         with open(os.path.join(folder, particle.name.replace(' ', '_') + "Densities_pi_gamma.txt"), 'a') as f1:
#             f1.write('{T:e}\t{a:e}\t{nHNL:e}\t{npion:e}\t{nnuel:e}\n'
#                 .format(T=particle.params.T / UNITS.MeV,
#                         a=particle.params.a,
#                         nHNL=com_den_HNL_1st / UNITS.MeV**3,
#                         npion=com_den_pion_1st / UNITS.MeV**3,
#                         nnuel=com_den_nuel_1st / UNITS.MeV**3))


universe.step_monitor = step_monitor
universe.evolve(T_final)