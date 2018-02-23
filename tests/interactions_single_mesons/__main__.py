# -*- coding: utf-8 -*-
"""

<img src="plots.svg" width=100% />
<img src="particles.svg" width=100% />

"""
import sys
import os
sys.path.insert(0, "/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[:-2]))
sys.settrace
import argparse
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, utils, LinearSpacedGrid, HeuristicGrid


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', default=300)#required=True)
parser.add_argument('--theta', default=0.0001)#required=True)
#parser.add_argument('--tau', required=True)
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)
#lifetime = float(args.tau) * UNITS.s

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    'output',
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)#, lifetime / UNITS.s)
    + args.comment
)


T_initial = 50. * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.0003125)

universe = Universe(params=params, folder=folder)

linear_grid_s = LinearSpacedGrid(MOMENTUM_SAMPLES=51, MAX_MOMENTUM=1*UNITS.MeV)
heuristic_grid = HeuristicGrid(mass, 0.35 * UNITS.s)
print(heuristic_grid.TEMPLATE*0.02)
photon = Particle(**SMP.photon)

#electron = Particle(**SMP.leptons.electron)
#muon = Particle(**SMP.leptons.muon)
#tau = Particle(**SMP.leptons.tau)

neutrino_e = Particle(**SMP.leptons.neutrino_e, **{'grid': heuristic_grid})
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


def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)

interaction = NuI.sterile_hadrons_interactions(
        thetas = thetas, sterile=sterile,
        neutrinos = [neutrino_e],
        leptons = [],
        mesons = [neutral_pion],
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


#universe.init_kawano(electron=electron, neutrino=neutrino_e)
#universe.init_oscillations(SMP.leptons.oscillations_map(), (neutrino_e, neutrino_mu, neutrino_tau))


def step_monitor(universe):
    import numpy
    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium:
            momenta = particle.grid.TEMPLATE 
            density = particle.density
            density_c = particle.density * particle.params.a**3 
            integrand = (particle.collision_integral * particle.params.H * particle.conformal_energy(particle.grid.TEMPLATE) / particle.mass /particle.params.a / particle._distribution)
            decay_rate = -integrand
            print((decay_rate/UNITS.MeV / 3.92472e-21))
            with open(os.path.join(folder, particle.name.replace(' ', '_') + ".decay_rate.txt"), 'a') as f1:
                f1.write('{:e}'.format(particle.params.T / UNITS.MeV) + '\t' + '{:e}'.format(particle.params.a) + '\t' + '\t'.join(['{:e}'.format(x) for x in decay_rate / UNITS.MeV]) + '\n')


universe.step_monitor = step_monitor

universe.evolve(T_final)

