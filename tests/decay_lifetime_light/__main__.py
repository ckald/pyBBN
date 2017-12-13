# -*- coding: utf-8 -*-


import argparse
import os
from collections import defaultdict

from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from library.NuMSM import particles as NuP, interactions as NuI
from evolution import Universe
from common import UNITS, Params, utils, LogSpacedGrid


parser = argparse.ArgumentParser(description='Run simulation for given mass and mixing angle')
parser.add_argument('--mass', default=33.9)
parser.add_argument('--theta', default=0.031)
parser.add_argument('--comment', default='')
args = parser.parse_args()

mass = float(args.mass) * UNITS.MeV
theta = float(args.theta)

folder = utils.ensure_dir(
    os.path.split(__file__)[0],
    "output",
    "mass={:e}_theta={:e}".format(mass / UNITS.MeV, theta)
    + args.comment
)

T_initial = 50. * UNITS.MeV
T_final = 0.0008 * UNITS.MeV
params = Params(T=T_initial,
                dy=0.003125)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)

# electron = Particle(**SMP.leptons.electron)

neutrino_e = Particle(**SMP.leptons.neutrino_e)
neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
neutrino_tau = Particle(**SMP.leptons.neutrino_tau)

for neutrino in [neutrino_e, neutrino_mu, neutrino_tau]:
    neutrino.decoupling_temperature = 0 * UNITS.MeV


sterile = Particle(**NuP.dirac_sterile_neutrino(mass))
sterile.decoupling_temperature = T_initial


universe.add_particles([
    photon,

    # electron,

    neutrino_e,
    neutrino_mu,
    neutrino_tau,

    sterile,
])

thetas = defaultdict(float, {
    'tau': theta,
})

interaction = NuI.sterile_leptons_interactions(
    thetas=thetas, sterile=sterile,
    neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
    leptons=[]
    # leptons=[electron]
)[0]


def reaction_type(reaction):
    return sum(reactant.side for reactant in reaction)


interaction.integrals = [integral for integral in interaction.integrals
                         if reaction_type(integral.reaction) in [-2, 2]]

universe.interactions += (
    interaction,
)


def step_monitor(universe):
    import numpy

    for particle in universe.particles:
        data = particle.data['params']
        if particle.mass > 0 and not particle.in_equilibrium and len(particle.data['params']) > 3:
            decay_rate = -(
                (numpy.log(data['density'][-1]) - numpy.log(data['density'][-2]))
                / (data['t'][-1] - data['t'][-2])
                + 3 * universe.params.H
            )

            print("{}: Γ ={: .3e} MeV, τ ={: .3e} s, Y = n/S ={: .3e}".format(
                particle.symbol,
                decay_rate / UNITS.MeV,
                1 / decay_rate / UNITS.s,
                particle.density * universe.params.a**3 / universe.params.S
            ))


universe.step_monitor = step_monitor


universe.evolve(T_final)
