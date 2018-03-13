# -*- coding: utf-8 -*-
"""

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

"""
import numpy
import sys
import os
sys.path.insert(0, "/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[:-2]))
sys.settrace
import argparse
from collections import defaultdict

os.environ['SPLIT_COLLISION_INTEGRAL'] = ''

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import CONST, UNITS, Params, utils, LinearSpacedGrid, HeuristicGrid
from interactions.four_particle.cpp.integral import CollisionIntegralKind

parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', default=300)
parser.add_argument('--theta', default=0.001)
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    'output',
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)
    + args.comment
)

T_initial = 50. * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.0003125)

universe = Universe(params=params, folder=folder)

linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=20*UNITS.MeV)
linear_grid = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=50*UNITS.MeV)
photon = Particle(**SMP.photon)

#electron = Particle(**SMP.leptons.electron)
#muon = Particle(**SMP.leptons.muon)
#tau = Particle(**SMP.leptons.tau)

neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': linear_grid})
#neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
#neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

neutral_pion = Particle(**SMP.hadrons.neutral_pion)
#charged_pion = Particle(**SMP.hadrons.charged_pion)

#neutral_rho = Particle(**SMP.hadrons.neutral_rho)
#charged_rho = Particle(**SMP.hadrons.charged_rho)

#eta = Particle(**SMP.hadrons.eta)
#eta_prime = Particle(**SMP.hadrons.eta_prime)
#omega = Particle(**SMP.hadrons.omega)
#phi = Particle(**SMP.hadrons.phi)

sterile = Particle(**NuP.dirac_sterile_neutrino(mass), **{'grid': linear_grid_s})


sterile.decoupling_temperature = T_initial
for neutrino in [neutrino_e]:#, neutrino_mu, neutrino_tau]:
    neutrino.decoupling_temperature = 1 * UNITS.MeV
#    neutrino.set_grid(grid)


universe.add_particles([
    photon,

#    electron,
#    muon,
#    tau,

    neutrino_e,
#    neutrino_mu,
#    neutrino_tau,

    neutral_pion,
#    charged_pion,

#    neutral_rho,
#    charged_rho,

#    eta,
#    eta_prime,
#    omega,
#    phi,

    sterile,
])

thetas = defaultdict(float, {
    'electron': theta,
})

# kind = CollisionIntegralKind.{}
# {} =
# Full: I_coll = A + f_1 * B
# F_1: I_coll = A
# F_f: I_coll = f_1 * B
# Full_vacuum_decay: I_coll = I_coll = A_vacuum_decay + f_1 * B_vacuum_decay
# F_1_vacuum_decay: I_coll = A_vacuum_decay
# F_f_vacuum_decay: I_coll = f_1 * B_vacuum_decay
kind = CollisionIntegralKind.F_f_vacuum_decay

def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)

interaction = NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],
        leptons = [],
        mesons = [neutral_pion],
        kind=kind
)[0]

#for inter in interaction:
#    inter.integrals = [integral for integral in inter.integrals
#                             if reaction_type(integral.reaction) in [1]]

interaction.integrals = [integral for integral in interaction.integrals
                            if reaction_type(integral.reaction) in [1]]

universe.interactions += (interaction, )

"""
universe.interactions += (
#    SMI.neutrino_interactions(
#        leptons=[electron, muon],
#        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau]
#    )
    # + NuI.sterile_leptons_interactions(
    #     thetas=thetas, sterile=sterile,
    #     neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    #     leptons=[electron, muon, tau]
    # )
    NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],# neutrino_mu, neutrino_tau],
        leptons = [],#[electron],# muon, tau],
        mesons = [neutral_pion#, charged_pion#, neutral_rho, charged_rho, eta,
#        eta_prime, omega, phi
        ],
    )
)
"""

theo_value = (CONST.G_F * theta * neutral_pion.decay_constant)**2 \
            * sterile.mass**3 * (1 - (neutral_pion.mass / sterile.mass)**2)**2 \
            / (16 * numpy.pi) / UNITS.MeV

def step_monitor(universe):
    for particle in universe.particles:
        if particle.mass > 0 and not particle.in_equilibrium:
            boost = particle.conformal_energy(particle.grid.TEMPLATE) / particle.conformal_mass
            integrand = (particle.collision_integral * particle.params.H * boost / particle.old_distribution)
            decay_rate = -integrand
            print('Average lifetime: {t:.6f} s\nΓ_code / Γ_theory:\n{rates}\n N/S: {NS:.6f}\n\n'
                .format(t=1 / decay_rate.mean() / UNITS.s,
                        rates=decay_rate / UNITS.MeV / theo_value,
                        NS=particle.density * universe.params.a**3 / universe.params.S))

            with open(os.path.join(folder, particle.name.replace(' ', '_') + "Decay_rate_three_particle.txt"), 'a') as f1:
                f1.write('{T:e}\t{a:e}\t{rates}\n'
                    .format(T=particle.params.T / UNITS.MeV,
                            a=particle.params.a,
                            rates='\t'.join(['{:e}'.format(x / UNITS.MeV) for x in decay_rate])))

universe.step_monitor = step_monitor

universe.evolve(T_final)